// meshseal web app — main-thread controller.
//
// Owns the UI; offloads the actual repair to a Web Worker (worker.js) so a
// multi-second repair never freezes the page. Talks to the worker with three
// message types: 'progress', 'done', 'error'.

'use strict';

const $ = (id) => document.getElementById(id);

const drop = $('drop');
const fileInput = $('file');
const filerow = $('filerow');
const fnameEl = $('fname');
const fmetaEl = $('fmeta');
const repairBtn = $('repair');
const bigwarn = $('bigwarn');
const progressWrap = $('progressWrap');
const barfill = $('barfill');
const stagetext = $('stagetext');
const cancelBtn = $('cancelBtn');
const resultEl = $('result');
const verdictEl = $('verdict');
const statsEl = $('stats');
const notesEl = $('notes');
const downloadBtn = $('download');
const againBtn = $('again');
const errorEl = $('error');
const previewWrap = $('previewWrap');
const previewBar = $('previewBar');
const previewNote = $('previewNote');
const viewerEl = $('viewer');
const btnBefore = $('btnBefore');
const btnAfter = $('btnAfter');
let previewActive = false;   // true once an input mesh is shown in the viewer

// 3D preview is size-gated: rendering both before+after geometries for a
// very large mesh is memory-heavy, so above this triangle count we skip it
// with a note rather than risk an unresponsive tab.
const PREVIEW_MAX = 400000;

let currentFile = null;     // the chosen File
let worker = null;          // active Web Worker (one per repair)
let outBlobUrl = null;      // object URL for the repaired download
let outName = 'repaired.stl';
let outMime = 'model/stl';

function fileFormat(name) { return /\.3mf$/i.test(name) ? '3mf' : 'stl'; }

// Binary STL is 84-byte header + 50 bytes/triangle. Use it as a rough
// face-count estimate for the size warning (ASCII is larger per face, so
// this over-estimates faces for ASCII — fine, the warning is conservative).
const BYTES_PER_TRI = 50;
const FACE_WARN = 1_500_000;   // warn above ~1.5M faces (~75 MB binary)

function fmt(n) { return n.toLocaleString('en-US'); }

function estFaces(bytes) { return Math.max(0, Math.round((bytes - 84) / BYTES_PER_TRI)); }

// Stop any in-flight repair worker. Also the cancel mechanism. Terminating
// here (not just on done/error) prevents a stale worker from posting a
// result against a newer file after the user moves on.
function terminateWorker() { if (worker) { worker.terminate(); worker = null; } }

function resetUI() {
  terminateWorker();
  progressWrap.style.display = 'none';
  resultEl.style.display = 'none';
  errorEl.style.display = 'none';
  bigwarn.style.display = 'none';
  barfill.style.width = '0%';
  previewWrap.classList.add('hidden');
  previewActive = false;
  btnAfter.disabled = true;
  btnBefore.classList.add('active');
  btnAfter.classList.remove('active');
  if (window.MeshPreview) window.MeshPreview.reset();
  if (outBlobUrl) { URL.revokeObjectURL(outBlobUrl); outBlobUrl = null; }
}

function pickFile(file) {
  if (!file) return;
  currentFile = file;
  resetUI();
  fnameEl.textContent = file.name;
  if (fileFormat(file.name) === '3mf') {
    // 3MF is zip-compressed XML — file size is no proxy for triangle count,
    // so don't estimate or warn (the input sanity caps still apply at parse).
    fmetaEl.textContent = `${(file.size / 1024).toFixed(0)} KB · 3MF model`;
  } else {
    const faces = estFaces(file.size);
    fmetaEl.textContent = `${(file.size / 1024).toFixed(0)} KB · ~${fmt(faces)} triangles`;
    if (faces > FACE_WARN) {
      bigwarn.textContent =
        `⚠ This is a large mesh (~${fmt(faces)} triangles). Repair may be slow or ` +
        `run out of browser memory. If it fails, try the command-line tool.`;
      bigwarn.style.display = 'block';
    }
  }
  filerow.style.display = 'flex';
  repairBtn.disabled = false;
  previewInput(file);   // show the input mesh right away (STL, under cap)
}

// ---- drag & drop + file picker -------------------------------------------
drop.addEventListener('click', () => fileInput.click());
fileInput.addEventListener('change', (e) => pickFile(e.target.files[0]));
['dragenter', 'dragover'].forEach((t) =>
  drop.addEventListener(t, (e) => { e.preventDefault(); drop.classList.add('over'); }));
['dragleave', 'drop'].forEach((t) =>
  drop.addEventListener(t, (e) => { e.preventDefault(); drop.classList.remove('over'); }));
drop.addEventListener('drop', (e) => {
  const f = e.dataTransfer.files && e.dataTransfer.files[0];
  if (f) pickFile(f);
});

// ---- repair --------------------------------------------------------------
repairBtn.addEventListener('click', startRepair);
againBtn.addEventListener('click', () => {
  resetUI();
  filerow.style.display = 'none';
  fileInput.value = '';
  currentFile = null;
});

let progressCount = 0;
const EXPECTED_STAGES = 24;   // ~ pipeline stage-boundary count; soft cap

async function startRepair() {
  if (!currentFile) return;
  // Partial reset: clear prior result/error but KEEP the input mesh preview
  // visible through the repair (don't call resetUI, which tears it down).
  resultEl.style.display = 'none';
  errorEl.style.display = 'none';
  if (outBlobUrl) { URL.revokeObjectURL(outBlobUrl); outBlobUrl = null; }
  repairBtn.disabled = true;
  progressWrap.style.display = 'block';
  stagetext.textContent = 'Loading engine…';
  progressCount = 0;
  barfill.style.width = '0%';

  const buf = await currentFile.arrayBuffer();
  const bytes = new Uint8Array(buf);
  const fmtType = fileFormat(currentFile.name);
  outName = currentFile.name.replace(/\.(stl|3mf)$/i, '') + '_fixed.' + fmtType;
  outMime = fmtType === '3mf' ? 'model/3mf' : 'model/stl';

  // Fresh worker per run — terminating it is also our cancel mechanism.
  terminateWorker();
  worker = new Worker('worker.js');
  worker.onmessage = onWorkerMessage;
  worker.onerror = (e) => showError(`Worker error: ${e.message || e}`);

  // Transfer the input buffer (zero-copy) into the worker. previewMax tells
  // the 3MF path whether to also return preview geometry (Route B).
  worker.postMessage({ type: 'repair', bytes, format: fmtType, previewMax: PREVIEW_MAX }, [bytes.buffer]);
}

function onWorkerMessage(ev) {
  const m = ev.data;
  if (m.type === 'progress') {
    progressCount++;
    const pct = m.stage === 'done' ? 100
      : Math.min(95, Math.round((progressCount / EXPECTED_STAGES) * 100));
    barfill.style.width = pct + '%';
    stagetext.textContent = m.stage === 'done'
      ? 'Finishing…'
      : `Working… ${m.stage}  (${fmt(m.faceCount)} faces, ${Math.round(m.elapsedMs)} ms)`;
  } else if (m.type === 'done') {
    barfill.style.width = '100%';
    showResult(m);
    maybePreview(m);
    terminateWorker();
  } else if (m.type === 'error') {
    showError(m.message);
    terminateWorker();
  }
}

cancelBtn.addEventListener('click', () => {
  terminateWorker();
  progressWrap.style.display = 'none';
  repairBtn.disabled = false;   // back to the pre-repair state (file still selected)
});

function showResult(m) {
  progressWrap.style.display = 'none';
  const s = m.summary;

  // Watertight is the property that decides whether the mesh is sliceable,
  // so it drives the headline. partial_failure / low confidence mean the
  // pipeline leaned on a fallback — worth disclosing as a caveat, but not a
  // reason to call a genuinely watertight+manifold result "failed".
  let cls, label, caveat = '';
  if (s.outputFaces === 0) {
    cls = 'bad'; label = '✗ Could not repair (empty result)';
  } else if (!s.watertight) {
    cls = 'partial'; label = '⚠ Partial repair — still has open or non-manifold edges';
  } else {
    cls = 'ok'; label = '✓ Watertight & manifold — ready to slice';
    if (s.partialFailure || s.confidence < 0.99) {
      caveat = `Confidence ${(s.confidence * 100).toFixed(0)}% — some regions needed a ` +
        `fallback repair, so the result is closed but may differ from the original there.`;
    }
  }
  verdictEl.className = 'verdict ' + cls;
  verdictEl.textContent = label;
  const caveatEl = $('caveat');
  caveatEl.textContent = caveat;
  caveatEl.classList.toggle('hidden', !caveat);

  const rows = [
    ['Watertight', s.watertight ? 'yes' : 'no'],
    ['Solid volume', s.isVolume ? 'yes' : 'no'],
    ['Components', fmt(s.components)],
  ];
  if (s.volumes > 1) rows.push(['Volumes (3MF, repaired separately)', fmt(s.volumes)]);
  rows.push(
    ['Triangles (in → out)', `${fmt(s.inputFaces)} → ${fmt(s.outputFaces)}`],
    ['Self-intersections', s.selfIntersections < 0 ? 'not measured' : fmt(s.selfIntersections)],
    ['Confidence', (s.confidence * 100).toFixed(0) + '%'],
    ['Repair time', Math.round(s.ms) + ' ms'],
  );
  statsEl.innerHTML = rows.map(
    ([k, v]) => `<tr><th>${k}</th><td class="num">${v}</td></tr>`).join('');

  notesEl.textContent = (s.notes && s.notes.length)
    ? 'Stages: ' + (s.stages || []).join(' · ') + '\n' + s.notes.join('\n')
    : '';
  $('noteslabel').classList.toggle('hidden', !notesEl.textContent);

  // Prepare the download.
  if (outBlobUrl) URL.revokeObjectURL(outBlobUrl);
  const blob = new Blob([m.bytes], { type: outMime });
  outBlobUrl = URL.createObjectURL(blob);
  downloadBtn.disabled = (s.outputFaces === 0);

  resultEl.style.display = 'block';
  repairBtn.disabled = false;
}

// ---- 3D before/after preview --------------------------------------------
function showViewer() {
  previewNote.classList.add('hidden');
  previewBar.classList.remove('hidden');
  viewerEl.classList.remove('hidden');
  previewWrap.classList.remove('hidden');
}
function showPreviewNote(text) {
  previewBar.classList.add('hidden');
  viewerEl.classList.add('hidden');
  previewNote.textContent = text;
  previewNote.classList.remove('hidden');
  previewWrap.classList.remove('hidden');
}

// Show the input ("Before") mesh immediately when a file is picked — STL
// only, and only when the file-size triangle estimate is under the cap.
// 3MF / large meshes get a note after repair instead.
async function previewInput(file) {
  if (!window.MeshPreview) return;
  if (fileFormat(file.name) !== 'stl') return;
  if (estFaces(file.size) > PREVIEW_MAX) return;
  try {
    const buf = await file.arrayBuffer();
    if (file !== currentFile) return;   // a newer file was picked while reading
    window.MeshPreview.loadBefore(buf);
    btnAfter.disabled = true;
    setPreviewMode('before');
    showViewer();
    previewActive = true;
  } catch (e) { /* leave preview hidden; repair still works */ }
}

// When the repair finishes, add the repaired ("After") mesh.
//   STL  — parse the input + output STL bytes in JS (STLLoader).
//   3MF  — use the before/after geometry the engine returned (Route B): no
//          3MF parser in the browser. The 3MF preview appears post-repair
//          (both meshes at once), since the engine produces the geometry.
async function maybePreview(m) {
  const s = m.summary;
  if (!window.MeshPreview) return;
  if (s.outputFaces === 0) return;
  const fmtType = fileFormat(currentFile.name);
  const tooBig = Math.max(s.inputFaces, s.outputFaces) > PREVIEW_MAX;
  try {
    if (fmtType === '3mf') {
      if (tooBig || !m.beforePositions || !m.afterPositions) {
        showPreviewNote(`Mesh too large for the in-browser 3D preview (> ${fmt(PREVIEW_MAX)} triangles). The repair is done — use the download.`);
        return;
      }
      window.MeshPreview.loadBeforePositions(m.beforePositions);
      window.MeshPreview.loadAfterPositions(m.afterPositions);
    } else {
      if (tooBig) {
        showPreviewNote(`Mesh too large for the in-browser 3D preview (> ${fmt(PREVIEW_MAX)} triangles). The repair is done — use the download.`);
        return;
      }
      if (!previewActive) {   // input wasn't shown at pick (race / earlier skip)
        window.MeshPreview.loadBefore(await currentFile.arrayBuffer());
      }
      window.MeshPreview.loadAfter(m.bytes.buffer);   // repaired STL bytes
    }
    btnAfter.disabled = false;
    setPreviewMode('after');
    showViewer();
    previewActive = true;
  } catch (e) {
    showPreviewNote('3D preview unavailable for this file.');
  }
}

function setPreviewMode(mode) {
  if (mode === 'after' && btnAfter.disabled) return;   // no repaired mesh yet
  if (window.MeshPreview) window.MeshPreview.setMode(mode);
  btnBefore.classList.toggle('active', mode === 'before');
  btnAfter.classList.toggle('active', mode === 'after');
}
btnBefore.addEventListener('click', () => setPreviewMode('before'));
btnAfter.addEventListener('click', () => setPreviewMode('after'));

downloadBtn.addEventListener('click', () => {
  if (!outBlobUrl) return;
  const a = document.createElement('a');
  a.href = outBlobUrl; a.download = outName;
  document.body.appendChild(a); a.click(); a.remove();
});

function showError(msg) {
  progressWrap.style.display = 'none';
  errorEl.textContent = '✗ ' + msg;
  errorEl.style.display = 'block';
  repairBtn.disabled = false;
}

// ---- PWA: service worker + install affordance ----------------------------
if ('serviceWorker' in navigator) {
  window.addEventListener('load', () => {
    navigator.serviceWorker.register('sw.js').catch(() => { /* non-fatal */ });
  });
}

let deferredInstall = null;
const installBtn = $('install');
window.addEventListener('beforeinstallprompt', (e) => {
  // Chrome/Edge: stash the event and surface our own install button.
  e.preventDefault();
  deferredInstall = e;
  if (installBtn) installBtn.classList.remove('hidden');
});
if (installBtn) {
  installBtn.addEventListener('click', async () => {
    if (!deferredInstall) return;
    deferredInstall.prompt();
    await deferredInstall.userChoice.catch(() => {});
    deferredInstall = null;
    installBtn.classList.add('hidden');
  });
}
window.addEventListener('appinstalled', () => {
  deferredInstall = null;
  if (installBtn) installBtn.classList.add('hidden');
});
