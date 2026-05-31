// Minimal Node harness for the meshseal WASM Phase-0 spike.
//
//   node harness.cjs <input.stl> [output.stl]
//
// Loads the Emscripten module, repairs the STL entirely in the WASM heap,
// prints the RepairResult summary, and writes the repaired STL out.

const fs = require('fs');
const path = require('path');

const createMeshseal = require(path.join(__dirname, '..', 'build-wasm', 'wasm', 'meshseal.js'));

async function main() {
  const inPath = process.argv[2];
  const outPath = process.argv[3] || '/tmp/wasm_repaired.stl';
  if (!inPath) {
    console.error('usage: node harness.cjs <input.stl> [output.stl]');
    process.exit(2);
  }

  const inputBytes = fs.readFileSync(inPath);
  const Module = await createMeshseal();

  const stageLog = [];
  const onProgress = (e) => { stageLog.push(e.stage); return true; };
  const is3mf = inPath.toLowerCase().endsWith('.3mf');

  const t0 = process.hrtime.bigint();
  const res = is3mf
    ? Module.repair3mf(new Uint8Array(inputBytes), onProgress, 2_000_000)
    : Module.repairStl(new Uint8Array(inputBytes), onProgress);
  const ms = Number(process.hrtime.bigint() - t0) / 1e6;
  if (is3mf && res.ok && res.beforePositions)
    console.log(`preview geom: before ${res.beforePositions.length/9} tris, after ${res.afterPositions.length/9} tris`);
  if (stageLog.length) console.log('progress:', stageLog.join(' -> '));

  if (!res.ok) {
    console.error('FAILED:', res.error);
    process.exit(1);
  }

  const out = Buffer.from(res.bytes);
  fs.writeFileSync(outPath, out);

  console.log(JSON.stringify({
    input: path.basename(inPath),
    inputBytes: inputBytes.length,
    inputFaces: res.inputFaces,
    outputFaces: res.outputFaces,
    watertight: res.watertight,
    isVolume: res.isVolume,
    components: res.components,
    selfIntersections: res.selfIntersections,
    confidence: res.confidence,
    partialFailure: res.partialFailure,
    outBytes: out.length,
    ms: Math.round(ms),
  }, null, 2));
  console.log('stages:', (res.stages || []).join(' '));
  if (res.notes && res.notes.length) console.log('notes :', res.notes.join(' | '));
}

main().catch(e => { console.error('harness error:', e); process.exit(1); });
