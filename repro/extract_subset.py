"""Reproducible response-blind subsets + full-band ephys from Chen et al. 2026.

Do not run with untrusted HDF5 files. Only HTTP 206 ranges are accepted.
Cell selection: 32 equal-index-spaced contiguous blocks of 32 source rows.
The source cell coordinate arrays are voxel footprints, not neuron centers;
centroids are computed with the pinned authors' pullCenters function.
"""
import concurrent.futures
import hashlib
import json
import pathlib
import sys
import time
import h5py
import numpy as np
import requests
import scipy.signal
from range_access import RangeFile

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent / "gutBrainPipeline"
sys.path.insert(0, str(REPO / "gutBrainCode"))
import gutBrainPipeline as source

REPO_COMMIT = "4c62612ecd540ec786c43fe9ce67cc78259a6ec1"


def fetch_exact(url, begin, end):
    for attempt in range(5):
        # Resolve a fresh Figshare redirect for each range; signed target URLs expire.
        r = requests.get(url, params={"range_start": begin, "attempt": attempt},
                         headers={"Range": f"bytes={begin}-{end}"}, timeout=90, stream=True)
        if r.status_code == 206:
            cr = r.headers.get("Content-Range", "")
            if not cr.startswith(f"bytes {begin}-{end}/"):
                r.close()
                raise ValueError(f"Wrong Content-Range {cr}")
            out = r.content
            if len(out) != end - begin + 1: raise ValueError("Truncated range")
            return out
        r.close()
        if attempt == 4: raise RuntimeError(f"Range request status {r.status_code}")
        time.sleep(attempt + 1)


def extract_cells(file_info, fish):
    remote = RangeFile(file_info["id"], file_info["size"], HERE / "range-cache")
    schema = {}
    with h5py.File(remote, "r") as h:
        n, nt = int(h["n"][()]), int(h["t"][()])
        dimensions = [int(h[k][()]) for k in ("x", "y", "z")]
        for key in ("cell_timeseries", "cell_baseline", "cell_x", "cell_y", "cell_z"):
            d = h[key]
            assert d.chunks is None and d.compression is None
            schema[key] = {"shape": d.shape, "dtype": d.dtype.str, "offset": d.id.get_offset()}
    starts = np.linspace(0, n - 32, 32, dtype=int)
    cell_ids = np.concatenate([np.arange(s, s + 32) for s in starts])
    arrays = {}
    tasks = []
    total_bytes = 0
    for key, spec in schema.items():
        row_bytes = int(np.prod(spec["shape"][1:])) * np.dtype(spec["dtype"]).itemsize
        arrays[key] = np.empty((len(cell_ids), *spec["shape"][1:]), dtype=spec["dtype"])
        for block, start in enumerate(starts):
            begin = spec["offset"] + int(start) * row_bytes
            end = begin + 32 * row_bytes - 1
            tasks.append((key, block, begin, end))
            total_bytes += end - begin + 1
    print(f"{fish}: {n} source cells x {nt} frames; extracting {len(cell_ids)} cells, {total_bytes / 1e6:.2f} MB", flush=True)
    if total_bytes > 60_000_000: raise ValueError("Unexpected cell extraction size")
    def one(task):
        key, block, begin, end = task
        cache = HERE / "subset-cache" / f"{file_info['id']}-{key}-{block}.bin"
        cache.parent.mkdir(exist_ok=True)
        if cache.exists(): data = cache.read_bytes()
        else:
            data = fetch_exact(file_info["download_url"], begin, end)
            cache.write_bytes(data)
        spec = schema[key]
        return key, block, np.frombuffer(data, dtype=spec["dtype"]).reshape((32, *spec["shape"][1:]))
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for i, (key, block, values) in enumerate(pool.map(one, tasks)):
            arrays[key][block * 32:(block + 1) * 32] = values
            if (i + 1) % 32 == 0: print(f"  completed {i+1}/{len(tasks)} cell ranges", flush=True)
    coords = source.pullCenters(arrays["cell_x"], arrays["cell_y"], arrays["cell_z"])
    with np.errstate(divide="ignore", invalid="ignore"):
        traces = (arrays["cell_timeseries"] - arrays["cell_baseline"]) / (arrays["cell_baseline"] - 100)
    return {"cell_ids": cell_ids, "traces": traces, "coords": coords,
            "fluorescence": arrays["cell_timeseries"], "baseline": arrays["cell_baseline"],
            "baseline_min": arrays["cell_baseline"].min(axis=1)}, {
            "source_n_cells": n, "source_n_frames": nt, "dimensions_xyz_voxels": dimensions,
            "selection": "32 blocks of 32 consecutive cell IDs, starts np.linspace(0,n-32,32,dtype=int); independent of neural responses",
            "selected_cell_ids": cell_ids.tolist(), "raw_range_bytes": total_bytes,
            "schema": schema, "cell_source": file_info,
            "coordinate_method": "Authors' pullCenters: unweighted mean of nonnegative footprint voxels on each axis; raw per-fish voxel coordinates, NOT registered anatomical labels",
            "dff_method": "Authors' (cell_timeseries-cell_baseline)/(cell_baseline-100); source fixed camera background100, not HDF5 background estimate"}


def extract_ephys(file_info, fish, expected_frames):
    path = HERE / file_info["name"]
    if path.stat().st_size != file_info["size"]: raise ValueError("Incomplete ephys file")
    digest = hashlib.md5()
    with path.open("rb") as fd:
        for chunk in iter(lambda: fd.read(8 << 20), b""): digest.update(chunk)
    assert digest.hexdigest() == file_info["supplied_md5"]
    # This is the authors' binary layout:26 sample-interleaved float32 channels at6kHz.
    a = source.load(str(path), num_channels=26, memmap=True)
    sr = 6000
    lightsheet = np.flatnonzero(np.diff((a[2] > 3.8).astype(np.int8)) == -1)
    if len(lightsheet) != expected_frames:
        raise ValueError(f"Trigger count {len(lightsheet)} != imaging frames {expected_frames}")
    frame_time = lightsheet / sr
    # Prefer recorded UV pulses when present; blood recordings require the documented
    # digital trial-command fallback. In gut fish3 the command precedes UV by~1.45s.
    command_starts = np.flatnonzero(np.diff((a[17] > .5).astype(np.int8), prepend=0) == 1)
    uv_above = np.flatnonzero(a[10] > 3.8)
    analog_starts = uv_above[np.r_[True, np.diff(uv_above) > 3 * sr]] if len(uv_above) else np.array([], dtype=int)
    analog_ends = uv_above[np.r_[np.diff(uv_above) > 3 * sr, True]] if len(uv_above) else np.array([], dtype=int)
    if len(analog_starts):
        if len(analog_starts) != len(command_starts): raise ValueError("UV pulse/command event-count mismatch")
        if not np.array_equal(a[22, analog_starts], a[22, command_starts]): raise ValueError("UV pulse/command site-code mismatch")
        uv_starts = analog_starts
        timing_method = "Recorded analog ch_UV10 >3.8V; pulse trains split by >3s gaps; onset of each train. Train counts/site codes independently match ch_uvtrial17 commands."
    else:
        uv_starts = command_starts
        timing_method = "Digital ch_uvtrial17 rising edges >0.5; analog ch_UV10 has no >3.8V pulses. This is command timing, not independently verified optical delivery timing."
    uv_positions = a[22, uv_starts]
    reference_position = 3. if fish == "glucose_fish3" else 1.
    observed_noncontrol = uv_positions[(uv_positions != 0) & (uv_positions != reference_position)]
    if not len(observed_noncontrol): raise ValueError("No noncontrol position recorded")
    target_position = float(observed_noncontrol[0])
    events = [{"sample": int(s), "seconds": float(s / sr), "position": float(p),
               "condition": "reference_position" if p == reference_position else "selected_other_position" if p == target_position else "other_unmapped_position",
               "preceding_frame": int(np.searchsorted(lightsheet, s) - 1)} for s, p in zip(uv_starts, uv_positions)]
    positions = sorted(np.unique(uv_positions).tolist())
    if not set(positions).issubset({0., 1., 2., 3., 4., 5., 6., 7.}):
        raise ValueError(f"Unreviewed position codes {positions}")
    # Compatibility shim only: SciPy1.13 moved gaussian to scipy.signal.windows.
    # The actual full-band variance algorithm is executed from the pinned MIT source.
    scipy.signal.gaussian = scipy.signal.windows.gaussian
    motor = np.zeros(len(lightsheet))
    for channel in (0, 1):
        motor += source.windowed_variance(a[channel], fs=sr)[0][lightsheet] / 2
    visual = a[13, lightsheet].copy()
    gain = a[15, lightsheet].copy()
    return {"time": frame_time, "motor": motor, "visual": visual, "visual_gain": gain,
            "visual_open": (visual > 0) & (gain == 0),
            "visual_closed": (visual > 0) & (gain > 0),
            "gut_event_times": uv_starts[uv_positions == target_position] / sr,
            "nutrient_event_times": uv_starts[uv_positions == target_position] / sr,
            "control_event_times": uv_starts[uv_positions == reference_position] / sr,
            "other_event_times": uv_starts[(uv_positions != target_position) & (uv_positions != reference_position)] / sr,
            "all_event_times": uv_starts / sr, "all_event_positions": uv_positions,
            "trial_command_times": command_starts / sr,
            "selected_target_position": np.array(target_position),
            "reference_position": np.array(reference_position)}, {
            "ephys_source": file_info, "ephys_md5_verified": digest.hexdigest(), "sample_rate_hz": sr,
            "source_samples": int(a.shape[1]), "matched_frame_count": len(lightsheet),
            "frame_interval_seconds_median": float(np.median(np.diff(frame_time))),
            "events": events, "uv_event_method": timing_method,
            "trial_command_times": (command_starts / sr).tolist(),
            "recorded_uv_onset_times": (analog_starts / sr).tolist(),
            "recorded_uv_train_durations_seconds": ((analog_ends - analog_starts) / sr).tolist(),
            "selected_target_position": target_position,
            "reference_position": reference_position,
            "position_mapping_basis": "Blood recordings reference code1 follows authors' example notebook control convention. Gut glucose fish3 reference code3 is the initial recorded site; selected target code2 matches the authors' gut example convention but no experiment-specific parameters.pickle is deposited. These are position-code comparisons: anatomical identity is NOT independently verified and NOT claimed identical between fish. All other positions retained as nuisance events; position0 not interpreted as anatomy.",
            "motor_method": "source.windowed_variance on each complete6kHz voltage channel0/1, mean of channels sampled at matched LS falling edges;40ms Gaussian kernels,4ms sigma",
            "visual_method": "source channel13 velocity sampled at LS timestamps, recorded units; channel15 gain saved separately. visual_open=(velocity>0)&(gain==0), visual_closed=(velocity>0)&(gain>0), matching source definitions. Raw visual is NOT moving-time fraction until binarized then averaged.",
            "analog_uv_min_max": [float(a[10].min()), float(a[10].max())]}


def main(fish):
    metadata = json.loads((HERE / "figshare-31967640.json").read_text())
    files = metadata["files"]
    cell_info = next(f for f in files if f["name"] == fish + "_cells0_clean.hdf5")
    ep_info = next(f for f in files if f["name"].startswith(fish + "_") and "chFlt" in f["name"])
    arrays, cell_manifest = extract_cells(cell_info, fish)
    ephys, ep_manifest = extract_ephys(ep_info, fish, arrays["traces"].shape[1])
    arrays.update(ephys)
    target = HERE / (fish + "_subset.npz")
    np.savez_compressed(target, **arrays)
    manifest = {"fish": fish, "dataset_doi": metadata["doi"], "data_license": metadata["license"],
                "repo": "https://github.com/benJames1212/gutBrainPipeline", "commit": REPO_COMMIT,
                "source_code_executed": ["load", "pullCenters", "windowed_variance"],
                "subset_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
                **cell_manifest, **ep_manifest}
    (HERE / (fish + "_manifest.json")).write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"file": str(target), "shape": list(arrays["traces"].shape),
                      "npz_MB": target.stat().st_size / 1e6,
                      "nutrient_events": arrays["nutrient_event_times"].tolist(),
                      "control_events": arrays["control_event_times"].tolist(),
                      "finite_cells": int(np.isfinite(arrays["traces"]).all(axis=1).sum()),
                      "positive_denominator_cells": int((arrays["baseline_min"] > 100).sum())}, indent=2), flush=True)


if __name__ == "__main__": main(sys.argv[1])
