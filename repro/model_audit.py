"""Small, inspectable held-out-trial audit for external zebrafish recordings.

This is an extension, NOT an exact reproduction of Chen et al. (2026).
The zero-padded lag construction follows the MIT-licensed gutBrainPipeline
constructLag/makeInputMat, commit 4c62612ecd540ec786c43fe9ce67cc78259a6ec1.
Unlike that repository's runBothModels, this module constructs lags on the
intact timeline before selecting rows. It explicitly uses D.T @ D penalties.
No cell is selected using a response threshold before analysis. The caller
must provide a response-independent cell sample and document its sampling.
"""
from __future__ import annotations

import numpy as np

HORIZONS = {"uv": 55.0, "nutrient": 55.0, "motor": 10.0, "visual": 5.0}
MODEL_KEYS = {"UVM": ("uv", "motor", "visual"),
              "UG": ("uv", "nutrient"),
              "FULL": ("uv", "nutrient", "motor", "visual")}
RIDGE, SMOOTH = 2.0, 20.0
POST_START, POST_END, EMBARGO = 2.0, 55.0, 55.0


def lag_matrix(values, order):
    """Zero-padded causal lags, equivalent to the cited official helper."""
    values = np.asarray(values, dtype=float)
    out = np.zeros((len(values), order), dtype=float)
    for lag in range(order):
        out[lag:, lag] = values[:len(values) - lag]
    return out


def penalty_matrix(orders, nuisance_columns=1):
    """Block-wise ridge + squared second-difference smoothness; intercept free."""
    p = nuisance_columns + sum(orders)
    penalty = np.eye(p) * RIDGE
    penalty[:nuisance_columns, :nuisance_columns] = 0
    offset = nuisance_columns
    for order in orders:
        d = np.diag(np.full(order, 2.0))
        if order > 1:
            d += np.diag(np.full(order - 1, -1.0), 1)
            d += np.diag(np.full(order - 1, -1.0), -1)
            d[0, 0] = d[-1, -1] = 1.0
        penalty[offset:offset + order, offset:offset + order] += SMOOTH * (d.T @ d)
        offset += order
    return penalty


def _event_rows(time, event, before=0.0, after=POST_END):
    return (time >= event["time"] - before) & (time < event["time"] + after)


def _score_rows(time, event):
    return (time >= event["time"] + POST_START) & (time < event["evaluation_end"])


def _training_mask(time, events, heldout):
    mask = np.ones(len(time), dtype=bool)
    for event in events:
        # Do not fit the UV artifact/copy-forward interval.
        mask &= ~_event_rows(time, event, after=POST_START)
        # Unmapped positions may have different responses. Do not force them
        # into the common UV-response kernel used for the two selected codes.
        if event["kind"] == "unclassified_uv":
            mask &= ~_event_rows(time, event, after=POST_END)
    for event in heldout:
        mask &= ~((time >= event["time"] - EMBARGO) &
                  (time < event["evaluation_end"] + EMBARGO))
    # Initial zero-padding is not an observed stimulus history.
    mask &= time >= time[0] + max(HORIZONS.values())
    return mask


def _corr_columns(actual, predicted):
    a = actual - np.mean(actual, axis=0, keepdims=True)
    b = predicted - np.mean(predicted, axis=0, keepdims=True)
    denom = np.sqrt(np.sum(a*a, axis=0) * np.sum(b*b, axis=0))
    return np.divide(np.sum(a*b, axis=0), denom,
                     out=np.full(actual.shape[1], np.nan), where=denom > 1e-12)


def _fit_models(signals, y, time, train_mask, drift_degree=3):
    """All signal scaling and coefficients are estimated on training rows only."""
    dt = float(np.median(np.diff(time)))
    raw_lagged = {key: lag_matrix(value, max(1, int(np.ceil(HORIZONS[key] / dt))))
                  for key, value in signals.items()}
    # Training excludes UV-artifact frames, so pulse amplitudes themselves are
    # absent there. Scale each entire lagged BLOCK on its training rows, not
    # just the lag-zero impulse column (which would incorrectly have RMS zero).
    scales = {key: max(float(np.sqrt(np.mean(value[train_mask] ** 2))), 1e-10)
              for key, value in raw_lagged.items()}
    lagged = {key: value / scales[key] for key,value in raw_lagged.items()}
    center = float(np.mean(time[train_mask]))
    time_scale = max(float(np.std(time[train_mask])), 1.0)
    normalized_time = (time - center) / time_scale
    nuisance = [np.ones(len(time))] + [normalized_time ** d for d in range(1,drift_degree+1)]
    scales.update(time_center_seconds=center, time_scale_seconds=time_scale,
                  time_polynomial_degree=drift_degree)
    predictions, coefficients, designs = {}, {}, {}
    for name, keys in MODEL_KEYS.items():
        keys = tuple(key for key in keys
                     if (np.any(np.abs(raw_lagged[key][train_mask]) > 1e-12)
                         if key in ("uv","nutrient")
                         else np.max(np.std(raw_lagged[key][train_mask],axis=0)) > 1e-12))
        design = np.column_stack(nuisance + [lagged[key] for key in keys])
        orders = [lagged[key].shape[1] for key in keys]
        x_train = design[train_mask]
        beta = np.linalg.solve(x_train.T @ x_train + penalty_matrix(orders,len(nuisance)),
                               x_train.T @ y[train_mask])
        predictions[name] = design @ beta
        coefficients[name] = beta
        designs[name] = design
    return predictions, coefficients, scales


def _score(y, predictions, rows):
    actual = y[rows]
    errors = {key: np.mean((pred[rows] - actual) ** 2, axis=0)
              for key, pred in predictions.items()}
    corrs = {key: _corr_columns(actual, pred[rows]) for key, pred in predictions.items()}
    def reduction(base):
        return np.divide(100 * (errors[base] - errors["FULL"]), errors[base],
                         out=np.full(actual.shape[1], np.nan), where=errors[base] > 1e-12)
    return {"nutrient_error_reduction_pct": reduction("UVM"),
            "multimodal_error_reduction_pct": reduction("UG"),
            "nutrient_delta_r": corrs["FULL"] - corrs["UVM"],
            "multimodal_delta_r": corrs["FULL"] - corrs["UG"],
            "full_r": corrs["FULL"], "uvm_r": corrs["UVM"], "ug_r": corrs["UG"],
            "full_mse": errors["FULL"], "uvm_mse": errors["UVM"], "ug_mse": errors["UG"]}


def _number(value):
    return round(float(value), 5) if np.isfinite(value) else None


def _array(values, digits=4):
    return [round(float(x), digits) if np.isfinite(x) else None for x in values]


def analyze_fish(time, traces, gut_event_times, control_event_times, motor, visual,
                 cell_ids, coords, *, recording_id="unspecified", route="unspecified",
                 max_display_points=1100, other_event_times=(),
                 selected_position_code=None, reference_position_code=1):
    """Return JSON-ready analysis. All event times use the supplied seconds axis.

    `gut_event_times` is retained for API compatibility; it means targeted
    nutrient uncaging and may refer to vascular delivery when `route` says so.
    `traces` are cells x imaging frames, already aligned/source-preprocessed.
    """
    time = np.asarray(time, dtype=float)
    traces = np.asarray(traces, dtype=float)
    motor, visual = np.asarray(motor, dtype=float), np.asarray(visual, dtype=float)
    coords, cell_ids = np.asarray(coords, dtype=float), np.asarray(cell_ids)
    if traces.shape != (len(cell_ids), len(time)) or coords.shape != (len(cell_ids), 3):
        raise ValueError("Expected traces[cells,time], coords[cells,3], and matching IDs")
    if len(time) < 200 or np.any(np.diff(time) <= 0):
        raise ValueError("Need at least 200 strictly increasing imaging timestamps")
    if motor.shape != time.shape or visual.shape != time.shape:
        raise ValueError("Motor and visual vectors must match imaging timestamps")
    if not all(np.isfinite(v).all() for v in (time, motor, visual)):
        raise ValueError("Time, motor and visual values must be finite")
    dt = float(np.median(np.diff(time)))
    if np.max(np.abs(np.diff(time) - dt)) > dt * 0.05:
        raise ValueError("Lag audit requires a regular imaging grid (<=5% timing jitter)")
    original_n_frames, original_dt = len(time), dt
    # Mean-bin to approximately 1 Hz: bounded design size, lower export weight,
    # and anti-aliasing rather than simply dropping motor/calcium samples.
    bin_width = max(1, int(np.ceil(1.0 / dt)))
    if bin_width > 1:
        n_bins = len(time) // bin_width
        end = n_bins * bin_width
        time = time[:end].reshape(n_bins, bin_width).mean(axis=1)
        traces = traces[:, :end].reshape(len(cell_ids), n_bins, bin_width).mean(axis=2)
        motor = motor[:end].reshape(n_bins, bin_width).mean(axis=1)
        visual = visual[:end].reshape(n_bins, bin_width).mean(axis=1)
        dt = float(np.median(np.diff(time)))
    valid = (np.isfinite(traces).all(axis=1) & np.isfinite(coords).all(axis=1)
             & (np.var(traces, axis=1) > 1e-12))
    n_input = len(cell_ids)
    traces, coords, cell_ids = traces[valid], coords[valid], cell_ids[valid]
    if not len(cell_ids):
        raise ValueError("No finite, nonconstant cells passed transparent input QC")
    y = traces.T
    events = []
    for kind, times in (("nutrient", gut_event_times), ("control", control_event_times),
                        ("unclassified_uv", other_event_times)):
        for j, onset in enumerate(np.sort(np.asarray(times, dtype=float))):
            if not np.isfinite(onset):
                raise ValueError("Event times must be finite")
            events.append({"id": f"{kind}-{j + 1}", "kind": kind, "time": float(onset),
                           "complete": False, "role": "excluded"})
    events.sort(key=lambda item: item["time"])
    for j, event in enumerate(events):
        next_time = events[j + 1]["time"] if j + 1 < len(events) else time[-1]
        event["evaluation_end"] = float(min(event["time"] + POST_END, next_time, time[-1]))
        duration = event["evaluation_end"] - event["time"] - POST_START
        complete = (event["time"] >= time[0] + 10 and duration >= 20
                    and event["kind"] != "unclassified_uv")
        preceded_unknown = any(e["kind"] == "unclassified_uv" and
                               e["time"] < event["evaluation_end"] and
                               e["time"] + POST_END > event["time"] + POST_START
                               for e in events)
        complete = complete and not preceded_unknown
        event.update(complete=bool(complete), role="discovery" if complete else "excluded",
                     scored_seconds=round(max(0.0, duration), 3),
                     overlaps_unclassified_response=bool(preceded_unknown))
    complete = {kind: [e for e in events if e["kind"] == kind and e["complete"]]
                for kind in ("nutrient", "control")}
    if not all(complete.values()):
        raise ValueError("Need at least one complete nutrient and UV-control event")
    signals = {key: np.zeros(len(time)) for key in ("uv", "nutrient")}
    for event in events:
        idx = int(np.argmin(np.abs(time - event["time"])))
        if time[0] <= event["time"] <= time[-1]:
            signals["uv"][idx] += 1
            if event["kind"] == "nutrient":
                signals["nutrient"][idx] += 1
    signals.update(motor=motor, visual=visual)
    heldout = complete["nutrient"][-2:] + complete["control"][-2:]
    train_mask = _training_mask(time, events, heldout)
    discovery_counts = {}
    discovery_frame_counts = {}
    for kind in complete:
        discovery_counts[kind] = sum(
            bool(train_mask[_score_rows(time, e)].all())
            for e in complete[kind] if e not in heldout)
        discovery_frame_counts[kind] = sum(
            int(np.sum(train_mask & _score_rows(time, e)))
            for e in complete[kind] if e not in heldout)
    # Partial earlier responses can train a fixed ridge model. Very small
    # support is disclosed, not used to justify selection on held-out data.
    independent = (min(map(len, complete.values())) >= 3
                   and min(discovery_frame_counts.values()) * dt >= 10
                   and train_mask.sum() >= max(80, int(60 / dt)))
    trial_scores, original_trial_scores, fold_records = [], [], []
    if independent:
        protocol_name = "late-within-condition independent holdout"
        predictions, betas, scales = _fit_models(signals, y, time, train_mask)
        original_predictions, _, _ = _fit_models(signals, y, time, train_mask,drift_degree=0)
        training_scores = _score(y, predictions, train_mask)
        # Freeze rank using training outcomes, never validation metrics.
        discovery_score = training_scores["nutrient_error_reduction_pct"]
        for event in heldout:
            event["role"] = "holdout"
            rows = _score_rows(time, event)
            trial_scores.append((event, _score(y, predictions, rows)))
            original_trial_scores.append((event, _score(y, original_predictions, rows)))
        fold_records.append({"heldout_event_ids": [e["id"] for e in heldout],
                             "train_frames": int(train_mask.sum()),
                             "signal_rms_training_only": {k: float(v) for k, v in scales.items()}})
        test_mask = np.zeros(len(time), bool)
        for event in heldout:
            test_mask |= _score_rows(time, event)
    else:
        protocol_name = "leave-event-out exploratory audit"
        heldout = complete["nutrient"] + complete["control"]
        predictions = {name: np.full_like(y, np.nan) for name in MODEL_KEYS}
        test_mask = np.zeros(len(time), bool)
        for event in heldout:
            event["role"] = "cross_validation"
            fold_train = _training_mask(time, events, [event])
            if fold_train.sum() < max(150, int(120 / dt)):
                event["role"] = "excluded_insufficient_training"
                continue
            fold_pred, _, scales = _fit_models(signals, y, time, fold_train)
            original_pred, _, _ = _fit_models(signals, y, time, fold_train,drift_degree=0)
            rows = _score_rows(time, event)
            for name in MODEL_KEYS:
                predictions[name][rows] = fold_pred[name][rows]
            trial_scores.append((event, _score(y, fold_pred, rows)))
            original_trial_scores.append((event, _score(y, original_pred, rows)))
            test_mask |= rows
            fold_records.append({"heldout_event_ids": [event["id"]],
                                 "train_frames": int(fold_train.sum()),
                                 "signal_rms_training_only": {k: float(v) for k,v in scales.items()}})
        discovery_score = np.zeros(len(cell_ids))
    if not trial_scores:
        raise ValueError("No test windows remain after QC and temporal buffering")
    nutrient_trials = [scores for event, scores in trial_scores if event["kind"] == "nutrient"]
    control_trials = [scores for event, scores in trial_scores if event["kind"] == "control"]
    if not nutrient_trials:
        raise ValueError("No nutrient test trials available")
    def average_scores(items):
        if not items:
            return {}
        result = {}
        for key in items[0]:
            values = np.stack([item[key] for item in items])
            count = np.isfinite(values).sum(axis=0)
            result[key] = np.divide(np.nansum(values,axis=0),count,
                                    out=np.full(values.shape[1],np.nan),where=count>0)
        return result
    # Average TRIAL metrics, not pseudo-independent calcium frames.
    mean_scores = average_scores(nutrient_trials)
    mean_control = average_scores(control_trials)
    original_mean = average_scores([scores for event,scores in original_trial_scores
                                    if event["kind"] == "nutrient"])
    if not independent:
        discovery_score = mean_scores["nutrient_error_reduction_pct"]
    ranking = np.argsort(-np.nan_to_num(discovery_score, nan=-np.inf), kind="stable")
    ranks = np.empty(len(ranking), int)
    ranks[ranking] = np.arange(1, len(ranking) + 1)
    positive_fraction = np.mean(np.stack([s["nutrient_error_reduction_pct"] > 0
                                          for s in nutrient_trials]), axis=0)
    positive_multimodal = np.mean(np.stack([s["multimodal_error_reduction_pct"] > 0
                                            for s in nutrient_trials]), axis=0)
    stride = max(1, int(np.ceil(len(time) / max_display_points)))
    display = np.arange(0, len(time), stride)
    cells = []
    for i, cell_id in enumerate(cell_ids):
        metrics = {key: _number(value[i]) for key, value in mean_scores.items()}
        metrics.update(nutrient_positive_trial_fraction=float(positive_fraction[i]),
                       multimodal_positive_trial_fraction=float(positive_multimodal[i]),
                       control_full_mse=_number(mean_control["full_mse"][i]) if mean_control else None)
        cells.append({"id": int(cell_id), "xyz": _array(coords[i], 3),
                      "discovery_rank": int(ranks[i]), "discovery_score": _number(discovery_score[i]),
                      "metrics": metrics,
                      "original_no_drift": {key: _number(value[i]) for key,value in original_mean.items()},
                      "drift_sensitivity": {
                          "position_gain_change_percentage_points": _number(mean_scores["nutrient_error_reduction_pct"][i] - original_mean["nutrient_error_reduction_pct"][i]),
                          "position_gain_sign_changed": bool((mean_scores["nutrient_error_reduction_pct"][i] > 0) != (original_mean["nutrient_error_reduction_pct"][i] > 0))},
                      "trace": _array(traces[i, display]),
                      "predictions": {name: _array(pred[display, i]) for name, pred in predictions.items()},
                      "trial_metrics": [{"event_id": event["id"], "kind": event["kind"],
                                          **{key: _number(value[i]) for key, value in score.items()}}
                                        for event, score in trial_scores]})
    return {
        "recording_id": recording_id, "route": route,
        "site_labels": {"selected": f"UV position code {selected_position_code}",
                        "reference": f"UV position code {reference_position_code}",
                        "selected_position_code": selected_position_code,
                        "reference_position_code": reference_position_code,
                        "anatomical_mapping_verified": False,
                        "interpretation": "position-associated prediction; not nutrient-specific or anatomical targeting proof"},
        "modalities": {
            "visual": {"varies": bool(np.ptp(visual) > 1e-9),
                       "min": float(np.min(visual)), "max": float(np.max(visual)),
                       "interpretation": "No varying visual stimulus: visual contribution is not identifiable and constant columns are omitted."
                                         if np.ptp(visual) <= 1e-9 else "Recorded visual regressor varies."},
            "motor": {"min": float(np.min(motor)), "max": float(np.max(motor)),
                      "std": float(np.std(motor)),
                      "interpretation": "Recorded motor-channel variance; low-amplitude fluctuations may be recording noise, not verified swimming."}},
        "protocol": {"name": protocol_name, "independent_shortlist_validation": bool(independent),
                     "rank_definition": "training-only nutrient error reduction" if independent else
                                        "exploratory leave-event-out nutrient error reduction; selection not independently validated",
                     "evaluation_window_seconds": [POST_START, POST_END],
                     "evaluation_censoring": "stop at next UV pulse or recording end; require at least20 scored seconds",
                     "temporal_embargo_seconds": EMBARGO, "lag_horizons_seconds": HORIZONS,
                     "lambda_ridge": RIDGE, "lambda_smooth": SMOOTH,
                     "penalty": "2I + 20 D.T@D, independent lag blocks; intercept and common time-polynomial nuisance unpenalized",
                     "common_session_drift": "degree3 polynomial included identically in all models; time centered/scaled on training rows only; degree fixed without response-based tuning",
                     "drift_sensitivity_audit": "original_no_drift retains the same split and models without timepolynomial; this sensitivity check was added after a blocked-design negative-control audit, not prospectively registered",
                     "scaling": "training-row RMS of each intact lagged regressor block, zero-preserving",
                     "analysis_binning": f"mean of {bin_width} consecutive sourceframes; visual value is moving-grating fraction within bin",
                     "error_reduction_definition": "100*(MSE baseline - MSE full)/MSE baseline; mean of per-trial values",
                     "warnings": ["Conditional predictability is not causal or a pure interoceptive signal.",
                                  "Internal nutrient field names refer to the selected position code; anatomical mapping is unverified.",
                                  "Unclassified UV-position response windows (55s) are excluded from training and evaluation.",
                                  "Position blocks can be confounded with session time, nonlinear drift and compound depletion; random interleaving is not demonstrated.",
                                  "A common cubic time nuisance tests slow-drift sensitivity but cannot remove all time/condition confounding.",
                                  "Two held-out trials per condition are a small within-recording check, not a statistically validated cell population.",
                                  "Motor activity can mediate nutrient responses, not only confound them.",
                                  "Late-within-condition holdout is not prospective whole-session forecasting.",
                                  "Source-preprocessed fluorescence is not a validated online preprocessing system.",
                                  "Cells are not independent animal replicates; coordinates are native recording space."],
                     "source_code_commit": "4c62612ecd540ec786c43fe9ce67cc78259a6ec1"},
        "counts": {"input_cells": n_input, "analyzed_cells": len(cell_ids),
                   "excluded_cells_nonfinite_or_constant": int(n_input - len(cell_ids)),
                   "original_frames": original_n_frames, "analysis_frames": len(time),
                   "display_frames": len(display),
                   "complete_nutrient_events": len(complete["nutrient"]),
                   "complete_control_events": len(complete["control"]),
                   "validation_nutrient_events": len(nutrient_trials),
                   "validation_control_events": len(control_trials),
                   "discovery_events_with_full_buffer": discovery_counts,
                   "discovery_response_frames": discovery_frame_counts,
                   "discovery_response_seconds": {k: round(v*dt,3) for k,v in discovery_frame_counts.items()}},
        "sample_interval_seconds": dt, "source_sample_interval_seconds": original_dt,
        "analysis_bin_width": bin_width, "display_stride": stride,
        "time": _array(time[display], 4), "motor": _array(motor[display], 7),
        "visual": _array(visual[display], 4), "events": events,
        "folds": fold_records, "validation_mask": test_mask[display].astype(int).tolist(),
        "cells": cells,
    }
