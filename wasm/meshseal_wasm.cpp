// WebAssembly entry shim for meshseal.
//
// Exposes a single embind function, repairStl(Uint8Array) -> object, that
// runs the full repair pipeline entirely in the browser / Node WASM heap:
//
//     read_stl_bytes  ->  repair()  ->  write_stl_bytes
//
// The return object carries the repaired STL bytes plus the RepairResult
// summary fields the web UI needs (watertight, component count, self-
// intersection count, confidence, notes). On any thrown exception the
// object instead carries { ok:false, error }.
//
// This file is compiled ONLY by the Emscripten toolchain (guarded by
// MESHSEAL_BUILD_WASM in CMake); it is never part of the native library,
// CLI, or test builds. It reaches into the library's private src/ header
// for the buffer-I/O entry points (read_stl_bytes / write_stl_bytes),
// which is fine because the shim is built inside the meshseal tree, not
// as an external consumer.

#include <emscripten/bind.h>
#include <emscripten/val.h>

#include <meshseal/meshseal.h>
#include "stl_io.h"   // meshseal::read_stl_bytes / write_stl_bytes
#include "3mf_io.h"   // meshseal::read_3mf_volumes_bytes / write_3mf*_bytes

#include <algorithm>
#include <cstdint>
#include <exception>
#include <string>
#include <vector>

using emscripten::val;

namespace {

// Copy a JS Uint8Array (or any array-like of numbers) into a std::vector.
// convertJSArrayToNumberVector is embind's purpose-built bridge for exactly
// this — it reads the typed array directly and needs no exported heap view.
std::vector<uint8_t> bytes_from_js(const val& input) {
    return emscripten::convertJSArrayToNumberVector<uint8_t>(input);
}

// Build a JS Uint8Array that OWNS a copy of `data` (detached from the wasm
// heap, so it survives after `data` is freed on return).
val js_bytes_copy(const std::vector<uint8_t>& data) {
    val heap_view = val(emscripten::typed_memory_view(data.size(), data.data()));
    // slice() with no args copies into a fresh, JS-owned ArrayBuffer.
    return heap_view.call<val>("slice");
}

// Wire a JS progress callback (or null/undefined) to opts.on_progress.
// An optional prefix is prepended to each stage name (used to tag the
// volume index during multi-volume 3MF repair). A `false` return from the
// JS callback requests cooperative cancellation.
void attach_progress(meshseal::RepairOptions& opts, const val& cb,
                     const std::string& prefix = std::string()) {
    if (cb.isUndefined() || cb.isNull()) return;
    opts.on_progress = [cb, prefix](const meshseal::ProgressEvent& e) -> bool {
        val ev = val::object();
        ev.set("stage", prefix.empty() ? e.stage_name : (prefix + e.stage_name));
        ev.set("elapsedMs", e.elapsed_ms);
        ev.set("faceCount", static_cast<double>(e.face_count));
        ev.set("vertexCount", static_cast<double>(e.vertex_count));
        ev.set("depth", e.recursion_depth);
        val r = cb(ev);
        if (r.isUndefined() || r.isNull()) return true;
        return r.as<bool>();
    };
}

void set_string_array(val& obj, const char* key,
                      const std::vector<std::string>& v) {
    val arr = val::array();
    for (std::size_t i = 0; i < v.size(); ++i)
        arr.set(static_cast<unsigned>(i), v[i]);
    obj.set(key, arr);
}

// Build a JS Float32Array owning a copy of `data` (survives the C++ vector).
val js_floats_copy(const std::vector<float>& data) {
    val view = val(emscripten::typed_memory_view(data.size(), data.data()));
    return view.call<val>("slice");
}

// Flatten a Mesh into a non-indexed triangle-soup position array (9 floats
// per face), appended to `out`. This is the geometry the 3D preview wants —
// returning it from the engine means the browser needs no 3MF parser
// (Route B). float32 matches what a renderer uses anyway.
void append_positions(std::vector<float>& out, const meshseal::Mesh& m) {
    out.reserve(out.size() + m.faces.size() * 9);
    for (const auto& f : m.faces)
        for (int k = 0; k < 3; ++k) {
            const auto& v = m.vertices[f[k]];
            out.push_back(static_cast<float>(v[0]));
            out.push_back(static_cast<float>(v[1]));
            out.push_back(static_cast<float>(v[2]));
        }
}

} // namespace

// repairStl(Uint8Array stlBytes, progressCb) -> {
//     ok, bytes, watertight, isVolume, components, selfIntersections,
//     confidence, partialFailure, inputFaces, outputFaces, notes[], stages[]
// }
//
// progressCb is a JS function (or null/undefined to disable). It is called
// synchronously at every pipeline stage boundary with an object
//   { stage, elapsedMs, faceCount, vertexCount, depth }
// and may return false to request cooperative cancellation (repair() then
// returns early with partialFailure=true). From a Web Worker the callback
// typically postMessage()s the event to the main thread and returns true;
// cancellation in the MVP is handled by terminating the worker instead.
val repair_stl(val input, val progressCb) {
    val result = val::object();
    try {
        std::vector<uint8_t> buf = bytes_from_js(input);

        meshseal::Mesh in = meshseal::read_stl_bytes(buf.data(), buf.size());
        const std::size_t in_faces = in.faces.size();

        meshseal::RepairOptions opts;
        attach_progress(opts, progressCb);
        meshseal::RepairResult res = meshseal::repair(in, opts);

        std::vector<uint8_t> out = meshseal::write_stl_bytes(res.mesh);

        result.set("ok", true);
        result.set("bytes", js_bytes_copy(out));
        result.set("watertight", res.watertight);
        result.set("isVolume", res.is_volume);
        result.set("components", res.component_count);
        result.set("selfIntersections", res.self_intersections);
        result.set("confidence", res.confidence);
        result.set("partialFailure", res.partial_failure);
        result.set("inputFaces", static_cast<double>(in_faces));
        result.set("outputFaces", static_cast<double>(res.mesh.faces.size()));
        result.set("volumes", 1.0);
        set_string_array(result, "notes", res.notes);
        set_string_array(result, "stages", res.stages_applied);
    } catch (const std::exception& e) {
        result.set("ok", false);
        result.set("error", std::string(e.what()));
    } catch (...) {
        result.set("ok", false);
        result.set("error", std::string("unknown error"));
    }
    return result;
}

// repair3mf(Uint8Array threeMfBytes, progressCb) -> same shape as repairStl,
// plus a `volumes` count.
//
// Mirrors the CLI's 3MF -> 3MF path: read all volumes (PrusaSlicer-style
// multi-volume metadata is honoured), repair each one INDEPENDENTLY, and
// write a multi-<object> 3MF preserving the per-volume separation. A
// single-volume 3MF takes the single-mesh path. The summary fields are
// aggregated across volumes (watertight = all; components = sum;
// confidence = min; self-intersections summed, or "not measured" if any
// volume was over the measurement cap; partialFailure = any).
//
// previewMaxFaces (3rd arg, 0 = off): when > 0 and the input is under the
// cap, the result also carries flat before/after position arrays
// (`beforePositions`/`afterPositions`, 9 floats per triangle) for the 3D
// preview — Route B, so the browser needs no 3MF parser.
val repair_3mf(val input, val progressCb, val previewMaxFaces) {
    val result = val::object();
    try {
        std::vector<uint8_t> buf = bytes_from_js(input);
        auto volumes = meshseal::read_3mf_volumes_bytes(buf.data(), buf.size());
        if (volumes.empty())
            throw std::runtime_error("3MF contains no mesh data");

        std::size_t in_faces = 0;
        for (const auto& v : volumes) in_faces += v.mesh.faces.size();

        // Route-B preview geometry (see previewMaxFaces note above).
        const double pmax = (previewMaxFaces.isUndefined() || previewMaxFaces.isNull())
            ? 0.0 : previewMaxFaces.as<double>();
        const bool want_geom = pmax > 0.0 && static_cast<double>(in_faces) <= pmax;
        std::vector<float> beforePos, afterPos;

        bool all_watertight = true, all_volume = true, any_partial = false;
        int total_components = 0;
        long long si_sum = 0;
        bool si_measured = true;
        double min_conf = 1.0;
        std::size_t out_faces = 0;
        std::vector<std::string> notes;
        std::vector<std::string> stages;
        auto merge_stages = [&](const std::vector<std::string>& s) {
            for (const auto& st : s)
                if (std::find(stages.begin(), stages.end(), st) == stages.end())
                    stages.push_back(st);
        };
        auto fold = [&](const meshseal::RepairResult& sub) {
            all_watertight = all_watertight && sub.watertight;
            all_volume     = all_volume && sub.is_volume;
            any_partial    = any_partial || sub.partial_failure;
            total_components += sub.component_count;
            if (sub.self_intersections < 0) si_measured = false;
            else si_sum += sub.self_intersections;
            if (sub.confidence < min_conf) min_conf = sub.confidence;
            out_faces += sub.mesh.faces.size();
            merge_stages(sub.stages_applied);
        };

        std::vector<uint8_t> out_bytes;

        if (volumes.size() > 1) {
            std::vector<meshseal::ThreeMfVolume> out_vols;
            out_vols.reserve(volumes.size());
            for (std::size_t i = 0; i < volumes.size(); ++i) {
                meshseal::RepairOptions opts;
                const std::string prefix = "vol " + std::to_string(i + 1) + "/" +
                    std::to_string(volumes.size()) + ": ";
                attach_progress(opts, progressCb, prefix);
                meshseal::RepairResult sub = meshseal::repair(volumes[i].mesh, opts);
                fold(sub);
                if (want_geom) {
                    append_positions(beforePos, volumes[i].mesh);
                    append_positions(afterPos, sub.mesh);   // before sub.mesh is moved
                }
                std::string tag = "[vol " + std::to_string(i + 1);
                if (!volumes[i].name.empty()) tag += " " + volumes[i].name;
                tag += "] ";
                for (const auto& n : sub.notes) notes.push_back(tag + n);
                if (!sub.mesh.faces.empty()) {
                    meshseal::ThreeMfVolume tv;
                    tv.mesh = std::move(sub.mesh);
                    tv.name = volumes[i].name;
                    tv.extruder = volumes[i].extruder;
                    out_vols.push_back(std::move(tv));
                }
            }
            if (out_vols.empty())
                throw std::runtime_error("repair produced an empty result for every volume");
            out_bytes = meshseal::write_3mf_volumes_bytes(out_vols, true);
            result.set("volumes", static_cast<double>(volumes.size()));
        } else {
            meshseal::RepairOptions opts;
            attach_progress(opts, progressCb);
            meshseal::RepairResult sub = meshseal::repair(volumes.front().mesh, opts);
            fold(sub);
            if (want_geom) {
                append_positions(beforePos, volumes.front().mesh);
                append_positions(afterPos, sub.mesh);
            }
            notes = sub.notes;
            out_bytes = meshseal::write_3mf_bytes(sub.mesh);
            result.set("volumes", 1.0);
        }

        result.set("ok", true);
        result.set("bytes", js_bytes_copy(out_bytes));
        result.set("watertight", all_watertight);
        result.set("isVolume", all_volume);
        result.set("components", total_components);
        result.set("selfIntersections", si_measured ? static_cast<double>(si_sum) : -1.0);
        result.set("confidence", min_conf);
        result.set("partialFailure", any_partial);
        result.set("inputFaces", static_cast<double>(in_faces));
        result.set("outputFaces", static_cast<double>(out_faces));
        set_string_array(result, "notes", notes);
        set_string_array(result, "stages", stages);
        // Route-B preview geometry (only when under the cap).
        if (want_geom && !beforePos.empty() &&
            static_cast<double>(out_faces) <= pmax) {
            result.set("beforePositions", js_floats_copy(beforePos));
            result.set("afterPositions", js_floats_copy(afterPos));
        }
    } catch (const std::exception& e) {
        result.set("ok", false);
        result.set("error", std::string(e.what()));
    } catch (...) {
        result.set("ok", false);
        result.set("error", std::string("unknown error"));
    }
    return result;
}

EMSCRIPTEN_BINDINGS(meshseal_wasm) {
    emscripten::function("repairStl", &repair_stl);
    emscripten::function("repair3mf", &repair_3mf);
}
