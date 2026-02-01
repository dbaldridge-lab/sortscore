"""Regression guards for the CLI entry point.
"""

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

import sortscore.run_analysis as ra

def _dummy_experiment(tmp_path: Path):
    """Create a minimal experiment object that satisfies run_analysis.main expectations."""

    class DummyExperiment:
        def __init__(self, outdir: Path):
            self.output_dir = str(outdir)
            self.counts = {1: {1: pd.DataFrame({"variant_seq": ["AAA"]})}}
            self.variant_type = "dna"
            self.wt_seq = "AAG"
            self.min_pos = 1
            self.max_pos = 3
            self.experiment_name = "dummy"
            self.avg_method = "rep-weighted"
            self.mutagenesis_variants = None
            self.position_type = "aa"
            self.biophysical_prop = False

        def load_counts(self):
            return None

        @property
        def num_positions(self):
            return 3

        @property
        def num_aa(self):
            return 3

    return DummyExperiment(tmp_path)


def test_batch_flag_invokes_batch_mode(tmp_path, monkeypatch):
    """Ensure --batch routes to batch workflow and passes suffix through."""
    cfg = tmp_path / "batch.json"
    cfg.write_text(json.dumps({"experiment_configs": ["a.json", "b.json"]}))

    called = {}
    monkeypatch.setattr(
        ra, "run_batch_mode", lambda path, suffix=None: called.setdefault("args", (path, suffix))
    )
    monkeypatch.setattr(
        sys, "argv", ["sortscore", "--batch", "-c", str(cfg), "--suffix", "smoke"]
    )

    ra.main()

    assert called["args"] == (str(cfg), "smoke")


def test_cli_smoke_non_batch_fast(tmp_path, monkeypatch):
    """Guard the non-batch path without touching real data."""
    dummy = _dummy_experiment(tmp_path)
    monkeypatch.setattr(ra.ExperimentConfig, "from_json", staticmethod(lambda path: dummy))
    monkeypatch.setattr(ra, "ensure_output_subdirs", lambda outdir: None)

    # Return dummy score paths; only aa_scores is used here.
    aa_scores_file = tmp_path / "aa_scores.csv"
    #TODO: remove codon column? #17
    aa_scores_file.write_text("aa_seq_diff,annotate_aa,avgscore_rep_weighted\nK.1.K,synonymous,0.0\n")
    monkeypatch.setattr(
        ra, "run_variant_analysis_workflow", lambda experiment, outdir, suffix, logger: (None, aa_scores_file)
    )

    # Keep plotting fast/no-op.
    monkeypatch.setattr(ra, "plot_heatmap", lambda *a, **k: None)

    # Lightweight analysis logger stand-in
    monkeypatch.setattr(
        ra,
        "AnalysisLogger",
        lambda *a, **k: type(
            "L",
            (),
            {"add_error": staticmethod(lambda *a, **k: None), "finish": staticmethod(lambda: "log.txt")},
        ),
    )

    # Use a tiny config; contents are ignored because we stub from_json.
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps({"experiment_name": "dummy"}))
    monkeypatch.setattr(sys, "argv", ["sortscore", "-c", str(cfg), "--suffix", "smoke"])

    ra.main()  # Should complete without raising

