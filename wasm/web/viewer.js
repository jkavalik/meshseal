// meshseal web app — 3D before/after preview (three.js).
//
// Loaded as an ES module; exposes window.MeshPreview for the classic
// app.js. The input mesh is shown as soon as a file is picked ("Before"),
// and the repaired mesh is added when the repair finishes ("After"). The
// toggle flips which is visible; the camera stays put.
//
// Defect highlighting: each mesh is rendered as TWO overlapping meshes
// sharing one geometry — a FrontSide material in the state colour and a
// BackSide material in red. For a correctly-oriented, watertight surface
// you only ever see front faces (state colour) from outside. You see the
// red BackSide material exactly where something is wrong: through a hole
// (the far interior wall faces away from you) or across an inverted-normal
// patch (its front points inward, so its back faces the camera). So red =
// "hole or flipped normal" — and a clean repair shows little or no red.

import * as THREE from 'three';
import { STLLoader } from './vendor/STLLoader.js';
import { OrbitControls } from './vendor/OrbitControls.js';

const COL_BEFORE  = 0xf0b429; // amber — unrepaired input (front faces)
const COL_AFTER   = 0x2bd47d; // green — repaired result (front faces)
const COL_PROBLEM = 0xff6b6b; // red   — back faces (holes / inverted normals)

let renderer, scene, camera, controls, container;
let groupBefore = null, groupAfter = null;
let resizeObs = null, rafId = 0;
let mode = 'before';

function init() {
  container = document.getElementById('viewer');
  // preserveDrawingBuffer lets the rendered frame be read back / captured.
  renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
  container.appendChild(renderer.domElement);

  scene = new THREE.Scene();
  scene.background = new THREE.Color(0x0f1216);

  camera = new THREE.PerspectiveCamera(45, 1, 0.01, 1e7);

  scene.add(new THREE.HemisphereLight(0xffffff, 0x384048, 1.1));
  const key = new THREE.DirectionalLight(0xffffff, 1.5);
  key.position.set(1, 1.5, 1);
  scene.add(key);
  const fill = new THREE.DirectionalLight(0xffffff, 0.5);
  fill.position.set(-1, -0.6, -1);
  scene.add(fill);

  controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;

  resizeObs = new ResizeObserver(resize);
  resizeObs.observe(container);

  const loop = () => { rafId = requestAnimationFrame(loop); controls.update(); renderer.render(scene, camera); };
  loop();
}

function resize() {
  if (!renderer || !container) return;
  const w = container.clientWidth, h = container.clientHeight;
  if (!w || !h) return;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  renderOnce();
}

// Draw a frame immediately, independent of the rAF loop — guarantees a
// first paint and instant toggle repaint even if rAF is throttled
// (background tab).
function renderOnce() { if (renderer) renderer.render(scene, camera); }

const COL_OPEN_EDGE = 0x19e3ff; // cyan — open (boundary) edges / hole rims

// Find open (boundary) edges and return them as cyan LineSegments, or null
// if the mesh is closed. STLLoader gives a non-indexed soup with duplicated
// vertices, so we weld by quantized position (tol ~ bbox·1e-5), then an
// undirected edge touched by exactly one triangle is a boundary edge.
function buildBoundaryEdges(geom) {
  const pos = geom.attributes.position;
  const n = pos.count;                       // 3 vertices per triangle
  geom.computeBoundingBox();
  const bb = geom.boundingBox;
  const q = (bb.min.distanceTo(bb.max) || 1) * 1e-5;
  const vkey = new Array(n);
  for (let i = 0; i < n; i++) {
    vkey[i] = Math.round(pos.getX(i) / q) + '_' +
              Math.round(pos.getY(i) / q) + '_' +
              Math.round(pos.getZ(i) / q);
  }
  const edges = new Map();                   // undirected key -> {count, i0, i1}
  for (let t = 0; t < n; t += 3) {
    for (let e = 0; e < 3; e++) {
      const i0 = t + e, i1 = t + ((e + 1) % 3);
      const ka = vkey[i0], kb = vkey[i1];
      if (ka === kb) continue;               // skip zero-length (degenerate) edges
      const k = ka < kb ? ka + '|' + kb : kb + '|' + ka;
      const rec = edges.get(k);
      if (rec) rec.count++; else edges.set(k, { count: 1, i0, i1 });
    }
  }
  const verts = [];
  for (const rec of edges.values()) {
    if (rec.count !== 1) continue;
    verts.push(pos.getX(rec.i0), pos.getY(rec.i0), pos.getZ(rec.i0),
               pos.getX(rec.i1), pos.getY(rec.i1), pos.getZ(rec.i1));
  }
  if (!verts.length) return null;
  const lg = new THREE.BufferGeometry();
  lg.setAttribute('position', new THREE.Float32BufferAttribute(verts, 3));
  const lines = new THREE.LineSegments(lg,
    new THREE.LineBasicMaterial({ color: COL_OPEN_EDGE }));
  lines.renderOrder = 2;
  return lines;
}

// One geometry → a Group of {front-face mesh in `frontColor`, back-face
// mesh in red, + cyan open-edge lines}. The two meshes never draw the same
// triangle (a triangle is either front- or back-facing to the camera), so
// there's no z-fighting.
function makeGroup(geom, frontColor) {
  geom.computeVertexNormals();
  const g = new THREE.Group();
  g.add(new THREE.Mesh(geom, new THREE.MeshStandardMaterial({
    color: frontColor, side: THREE.FrontSide, flatShading: true, roughness: 0.6, metalness: 0.0,
  })));
  g.add(new THREE.Mesh(geom, new THREE.MeshStandardMaterial({
    color: COL_PROBLEM, side: THREE.BackSide, flatShading: true, roughness: 0.6, metalness: 0.0,
  })));
  const edges = buildBoundaryEdges(geom);    // hole rims / open boundaries
  if (edges) g.add(edges);
  return g;
}

function disposeGroup(g) {
  if (!g) return;
  scene.remove(g);
  const seenGeom = new Set();   // front+back meshes share one geometry; lines have their own
  g.traverse((o) => {
    if (o.geometry && !seenGeom.has(o.geometry)) { o.geometry.dispose(); seenGeom.add(o.geometry); }
    if (o.material) o.material.dispose();
  });
}

function fitCamera(geom) {
  geom.computeBoundingSphere();
  const s = geom.boundingSphere;
  const r = s.radius || 1;
  const c = s.center;
  controls.target.copy(c);
  const d = r * 2.6;
  camera.position.set(c.x + d * 0.4, c.y + d * 0.35, c.z + d);
  camera.near = Math.max(r / 1000, 1e-4);
  camera.far = r * 1000;
  camera.updateProjectionMatrix();
  controls.update();
}

function applyMode() {
  if (groupBefore) groupBefore.visible = (mode === 'before');
  if (groupAfter)  groupAfter.visible  = (mode === 'after');
}

function triCount(g) {
  if (!g) return 0;
  const m = g.children[0];
  return (m && m.geometry.attributes.position) ? m.geometry.attributes.position.count / 3 : 0;
}

// Build a non-indexed BufferGeometry from a flat triangle-soup position
// array (9 floats per face) — what the engine returns for 3MF (Route B).
function geomFromPositions(arr) {
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(arr, 3));
  return g;
}

// Shared loaders: take a BufferGeometry (from STL bytes or engine positions)
// and install it as the before / after group.
function setBeforeGeom(geom) {
  if (!renderer) init();
  disposeGroup(groupBefore); disposeGroup(groupAfter); groupBefore = groupAfter = null;
  groupBefore = makeGroup(geom, COL_BEFORE);
  scene.add(groupBefore);
  mode = 'before';
  applyMode();
  fitCamera(geom);
  resize();
  renderOnce();
}
function setAfterGeom(geom) {
  if (!renderer) init();
  disposeGroup(groupAfter); groupAfter = null;
  groupAfter = makeGroup(geom, COL_AFTER);
  scene.add(groupAfter);
  mode = 'after';
  applyMode();
  renderOnce();
}

window.MeshPreview = {
  // Show the input mesh immediately (called at file-pick). Resets any
  // prior preview and fits the camera to this geometry.
  loadBefore(buf) { setBeforeGeom(new STLLoader().parse(buf)); },
  // Add the repaired mesh and switch to it (called when repair finishes).
  loadAfter(buf) { setAfterGeom(new STLLoader().parse(buf)); },
  // Engine-geometry variants (3MF / Route B): `arr` is a Float32Array of
  // triangle-soup positions.
  loadBeforePositions(arr) { setBeforeGeom(geomFromPositions(arr)); },
  loadAfterPositions(arr) { setAfterGeom(geomFromPositions(arr)); },
  setMode(m) { mode = m; applyMode(); renderOnce(); },
  hasAfter() { return !!groupAfter; },
  reset() {
    disposeGroup(groupBefore); disposeGroup(groupAfter);
    groupBefore = groupAfter = null;
    mode = 'before';
    renderOnce();
  },
  info() { return { beforeTris: triCount(groupBefore), afterTris: triCount(groupAfter), mode }; },
};
