"""Reproduce the two response-independent, position-coded recording audits."""
import json
from pathlib import Path
import numpy as np
from model_audit import analyze_fish

ROOT = Path(__file__).resolve().parent
RECORDINGS = [
    ("blood_glucose_fish1", "fish1-analysis.json", 1, "Blood D-glucose"),
    ("blood_glucose_fish2", "fish2-analysis.json", 1, "Blood D-glucose"),
    ("glucose_fish3", "gut3-analysis.json", 3, "Gut D-glucose"),
]
for recording_id, output_name, reference_code, route in RECORDINGS:
    path = ROOT / f"{recording_id}_subset.npz"
    if not path.exists():
        print(f"Not available: {path.name}")
        continue
    with np.load(path) as z:
        valid = z["baseline_min"] > 100
        result = analyze_fish(
            z["time"], z["traces"][valid], z["gut_event_times"], z["control_event_times"],
            z["motor"], (z["visual"] > 0).astype(float), z["cell_ids"][valid], z["coords"][valid],
            other_event_times=z["other_event_times"], recording_id=recording_id,
            route=f"{route} dataset; UV positions anatomically unmapped",
            selected_position_code=int(z["selected_target_position"]), reference_position_code=reference_code)
        result["counts"]["sampled_cells_before_camera_baseline_qc"] = int(len(valid))
        result["counts"]["excluded_camera_baseline_at_or_below_100"] = int((~valid).sum())
        result["all_position_events"] = [
            {"time":float(time),"position_code":int(position)}
            for time,position in zip(z["all_event_times"],z["all_event_positions"])]
        if recording_id == "glucose_fish3":
            all_times, all_positions = z["all_event_times"],z["all_event_positions"]
            alternate = analyze_fish(
                z["time"],z["traces"][valid],all_times[all_positions == 2],all_times[all_positions == 1],
                z["motor"],(z["visual"]>0).astype(float),z["cell_ids"][valid],z["coords"][valid],
                other_event_times=all_times[(all_positions != 2) & (all_positions != 1)],
                recording_id="glucose_fish3_reference1",route="Gut D-glucose; alternate position-code reference",
                selected_position_code=2,reference_position_code=1)
            alternate_by_id={cell["id"]:cell for cell in alternate["cells"]}
            for cell in result["cells"]:
                other=alternate_by_id[cell["id"]]
                cell["alternative_reference"]={
                    "reference_position":1,"metrics":other["metrics"],
                    "original_no_drift":other["original_no_drift"],"trial_metrics":other["trial_metrics"]}
            result["alternative_reference"]={
                "reference_position":1,"protocol":alternate["protocol"],
                "trial_support":alternate["counts"],"folds":alternate["folds"],
                "site_labels":alternate["site_labels"],
                "rank_policy":"Primary discovery ranks unchanged; alternate comparison is a sensitivity check, never used to re-rank."}
            (ROOT / "gut3-reference1-analysis.json").write_text(
                json.dumps(alternate,separators=(",",":"),allow_nan=False))
    target = ROOT / output_name
    target.write_text(json.dumps(result,separators=(",",":"),allow_nan=False))
    print(target.name,result["protocol"]["name"],result["counts"])
    print("rank / cell / training gain % / held-out position gain % / full r / motor gain %")
    for cell in sorted(result["cells"], key=lambda item:item["discovery_rank"])[:10]:
        print(cell["discovery_rank"],cell["id"],cell["discovery_score"],
              cell["metrics"]["nutrient_error_reduction_pct"],cell["metrics"]["full_r"],
              cell["metrics"]["multimodal_error_reduction_pct"])
