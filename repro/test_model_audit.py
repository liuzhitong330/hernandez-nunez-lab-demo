"""Deterministic tests; synthetic signals are QA only, never the website data."""
import ast
from pathlib import Path
import unittest
import numpy as np
import model_audit as m


class ModelAuditTests(unittest.TestCase):
    def synthetic(self):
        rng = np.random.default_rng(901)
        time = np.arange(2600.0)
        nutrient = np.arange(200, 1801, 200.0)
        control = np.arange(100, 1701, 200.0)
        motor = rng.uniform(0, 1, len(time))
        visual = (np.sin(time / 17) > 0).astype(float)
        g = np.zeros(len(time)); u = g.copy()
        g[nutrient.astype(int)] = 1
        u[control.astype(int)] = 1
        u += g
        kernel = np.exp(-np.arange(55)/12)
        signal = np.convolve(g, kernel, mode="full")[:len(time)]
        uv = np.convolve(u, kernel, mode="full")[:len(time)]
        traces = np.vstack([signal + 0.1*uv, 0.3*motor + 0.3*signal,
                            0.2*visual + 0.2*signal, 0.2*uv])
        traces += rng.normal(0, .01, traces.shape)
        return time, traces, nutrient, control, motor, visual, np.arange(4), np.arange(12).reshape(4,3)

    def test_lag_matches_official_mit_helper(self):
        path = Path(__file__).parents[1] / "gutBrainPipeline/gutBrainCode/gutBrainPipeline.py"
        tree = ast.parse(path.read_text())
        selected = [node for node in tree.body if isinstance(node, ast.FunctionDef)
                    and node.name in ("shift", "constructLag")]
        scope = {"np": np}
        exec(compile(ast.Module(body=selected, type_ignores=[]), str(path), "exec"), scope)
        values = np.random.default_rng(44).normal(size=50)
        np.testing.assert_array_equal(m.lag_matrix(values, 9), scope["constructLag"](values, 9))

    def test_no_lags_cross_deleted_timeline(self):
        impulse = np.zeros(30); impulse[10] = 1
        lagged = m.lag_matrix(impulse, 6)
        mask = np.ones(30, bool); mask[11:20] = False
        # The event at 10 cannot acquire a lag-1 influence at timestamp 20.
        self.assertEqual(lagged[20, 1], 0)
        self.assertEqual(lagged[mask][11, 1], 0)

    def test_holdout_outcomes_do_not_change_prediction_or_rank(self):
        args = self.synthetic()
        first = m.analyze_fish(*args)
        self.assertTrue(first["protocol"]["independent_shortlist_validation"])
        self.assertGreater(first["folds"][0]["signal_rms_training_only"]["uv"], .001)
        self.assertGreater(first["folds"][0]["signal_rms_training_only"]["nutrient"], .001)
        changed = list(args)
        changed[1] = args[1].copy()
        for event in first["events"]:
            if event["role"] == "holdout":
                rows = (args[0] >= event["time"] - 55) & (args[0] < event["time"] + 110)
                changed[1][:, rows] += 7
        second = m.analyze_fish(*changed)
        for a,b in zip(first["cells"], second["cells"]):
            self.assertEqual(a["discovery_rank"], b["discovery_rank"])
            self.assertEqual(a["predictions"], b["predictions"])
            self.assertEqual(a["discovery_score"], b["discovery_score"])
        self.assertNotEqual(first["cells"][0]["metrics"], second["cells"][0]["metrics"])

    def test_penalty_has_no_between_modality_smoothing(self):
        penalty = m.penalty_matrix([4,3])
        self.assertEqual(penalty[4,5], 0)
        self.assertEqual(penalty[0,0], 0)
        self.assertGreaterEqual(np.min(np.linalg.eigvalsh(penalty)), 0)

    def test_event_windows_do_not_overlap(self):
        args = list(self.synthetic())
        args[2] = np.arange(420., 581, 40)
        args[3] = np.array([300., 340., 380., 620., 660., 700.])
        result = m.analyze_fish(*args, other_event_times=[740.])
        events = result["events"]
        for left,right in zip(events,events[1:]):
            self.assertLessEqual(left["evaluation_end"], right["time"])
        self.assertFalse(events[-1]["complete"])
        self.assertEqual(events[-2]["evaluation_end"],740.)

    def test_high_rate_is_mean_binned(self):
        args = list(self.synthetic())
        args[0] = args[0]*0.4
        args[2] = np.array([120.,200.,280.,360.,440.,520.])
        args[3] = np.array([160.,240.,320.,400.,480.,560.])
        result = m.analyze_fish(*args)
        self.assertEqual(result["analysis_bin_width"],3)
        self.assertEqual(result["counts"]["original_frames"],2600)
        self.assertAlmostEqual(result["sample_interval_seconds"],1.2)

    def test_unclassified_uv_response_is_censored(self):
        events = [{"time":100., "kind":"unclassified_uv", "evaluation_end":140.}]
        t = np.arange(300.)
        mask = m._training_mask(t, events, [])
        self.assertFalse(mask[100:155].any())
        self.assertTrue(mask[155])

    def test_polynomial_only_data_are_captured_without_position_terms(self):
        args=list(self.synthetic())
        u=(args[0]-args[0].mean())/args[0].std()
        args[1]=np.vstack([u,u*u,u**3,u+.4*u*u])
        result=m.analyze_fish(*args)
        for cell in result["cells"]:
            # Exact low-order drift must be explained without an event effect.
            self.assertLess(cell["metrics"]["uvm_mse"],1e-8)
            self.assertLess(cell["metrics"]["full_mse"],1e-8)

    def test_insufficient_trials_is_exploratory(self):
        args = list(self.synthetic())
        args[2],args[3] = args[2][:2],args[3][:2]
        result = m.analyze_fish(*args)
        self.assertFalse(result["protocol"]["independent_shortlist_validation"])
        self.assertIn("exploratory", result["protocol"]["name"])


if __name__ == "__main__":
    unittest.main()
