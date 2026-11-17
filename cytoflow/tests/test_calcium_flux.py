#!/usr/bin/env python3.8
# coding: latin-1

# (c) Massachusetts Institute of Technology 2015-2018
# (c) Brian Teague 2018-2022
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Tests for CalciumFluxOp and CalciumFluxView.
"""

import unittest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')

import cytoflow as flow
import cytoflow.utility as util
from cytoflow.operations.calcium_flux import CalciumFluxOp, CalciumFluxView


def _make_experiment(n_events=2000, baseline_end=60.0, spike_fraction=0.3,
                     spike_magnitude=3.0, seed=42):
    """
    Build a synthetic Experiment with a Time channel and a Fluo4_A channel.

    * Events are uniformly distributed in time from 0 to 300 s.
    * Pre-baseline events (time <= baseline_end) have Fluo4_A ~ N(100, 10).
    * Post-baseline events are split:
      - (1 - spike_fraction): Fluo4_A ~ N(100, 10)   (non-responding)
      - spike_fraction:       Fluo4_A ~ N(100 * spike_magnitude, 20)
    """
    rng = np.random.default_rng(seed)

    times = rng.uniform(0.0, 300.0, n_events)
    signal = np.full(n_events, np.nan)

    baseline_mask = times <= baseline_end
    signal[baseline_mask] = rng.normal(100.0, 10.0, baseline_mask.sum())

    post_mask = ~baseline_mask
    n_post = post_mask.sum()
    responder_mask = rng.random(n_post) < spike_fraction

    post_signal = np.where(
        responder_mask,
        rng.normal(100.0 * spike_magnitude, 20.0, n_post),
        rng.normal(100.0, 10.0, n_post),
    )
    signal[post_mask] = post_signal

    df = pd.DataFrame({'Time': times, 'Fluo4_A': signal})

    # Build a minimal Experiment by hand (ImportOp would need an FCS file).
    # Channels must be registered first (with data=None on an empty experiment),
    # then events are added via add_events.
    ex = flow.Experiment()
    ex.add_channel('Time', None)
    ex.add_channel('Fluo4_A', None)
    ex.add_events(df[['Time', 'Fluo4_A']], {})

    return ex


class TestCalciumFluxApply(unittest.TestCase):
    """Tests for CalciumFluxOp.apply()"""

    def setUp(self):
        self.ex = _make_experiment()
        self.op = CalciumFluxOp(
            name='Ca',
            time_channel='Time',
            signal_channel='Fluo4_A',
            baseline_end=60.0,
            response_threshold=1.5,
            bin_width=10.0,
        )

    def test_ratio_channel_created(self):
        """The {name}_ratio channel must be present in the result."""
        ex2 = self.op.apply(self.ex)
        self.assertIn('Ca_ratio', ex2.data.columns,
                      "Ca_ratio channel was not added to the experiment")

    def test_ratio_values_near_one_at_baseline(self):
        """Mean F/F0 during the baseline window should be close to 1."""
        ex2 = self.op.apply(self.ex)
        baseline_mask = ex2.data['Time'] <= self.op.baseline_end
        baseline_mean = ex2.data.loc[baseline_mask, 'Ca_ratio'].mean()
        self.assertAlmostEqual(baseline_mean, 1.0, places=1,
                               msg="Baseline mean ratio should be ~1")

    def test_time_bin_condition_created(self):
        """The {name}_time_bin condition must be present."""
        ex2 = self.op.apply(self.ex)
        self.assertIn('Ca_time_bin', ex2.data.columns,
                      "Ca_time_bin condition not found")

    def test_responding_condition_created(self):
        """The {name}_responding condition must exist with correct categories."""
        ex2 = self.op.apply(self.ex)
        self.assertIn('Ca_responding', ex2.data.columns,
                      "Ca_responding condition not found")
        unique_vals = set(ex2.data['Ca_responding'].unique())
        self.assertTrue(unique_vals <= {"Responding", "Not responding"},
                        "Unexpected values in Ca_responding: {}".format(unique_vals))

    def test_responding_fraction(self):
        """Post-baseline responding fraction should be non-trivial."""
        ex2 = self.op.apply(self.ex)
        post_mask = ex2.data['Time'] > self.op.baseline_end
        post_data = ex2.data[post_mask]
        pct = (post_data['Ca_responding'] == 'Responding').mean()
        # We spiked ~30 % of post-baseline cells at 3× baseline, well above 1.5×
        self.assertGreater(pct, 0.1,
                           "Expected > 10 % responding post-baseline, got {:.1%}"
                           .format(pct))

    def test_statistics_created(self):
        """All three statistics entries must be present."""
        ex2 = self.op.apply(self.ex)
        self.assertIn('Ca_mean_ratio', ex2.statistics,
                      "Ca_mean_ratio statistic not found")
        self.assertIn('Ca_sem_ratio', ex2.statistics,
                      "Ca_sem_ratio statistic not found")
        self.assertIn('Ca_pct_responding', ex2.statistics,
                      "Ca_pct_responding statistic not found")

    def test_statistics_indexed_by_time_bin(self):
        """Statistics must be indexed by Ca_time_bin."""
        ex2 = self.op.apply(self.ex)
        stat = ex2.statistics['Ca_mean_ratio']
        self.assertEqual(stat.index.name, 'Ca_time_bin',
                         "Statistics index should be named Ca_time_bin")

    def test_ratio_denominator(self):
        """When ratio_denominator is set the signal is a ratio before F/F0."""
        import copy
        ex_ratio = copy.deepcopy(self.ex)
        # add a constant denominator channel
        ex_ratio.add_channel('Denom_A',
                              pd.Series(np.full(len(ex_ratio.data), 50.0),
                                        index=ex_ratio.data.index))

        op2 = CalciumFluxOp(
            name='CaRatio',
            time_channel='Time',
            signal_channel='Fluo4_A',
            ratio_denominator='Denom_A',
            baseline_end=60.0,
            response_threshold=1.5,
            bin_width=10.0,
        )
        ex2 = op2.apply(ex_ratio)
        self.assertIn('CaRatio_ratio', ex2.data.columns)

    def test_error_no_name(self):
        """apply() must raise if name is empty."""
        self.op.name = ''
        with self.assertRaises(util.CytoflowOpError):
            self.op.apply(self.ex)

    def test_error_bad_time_channel(self):
        """apply() must raise if time_channel is not in experiment."""
        self.op.time_channel = 'NonExistent'
        with self.assertRaises(util.CytoflowOpError):
            self.op.apply(self.ex)

    def test_error_bad_signal_channel(self):
        """apply() must raise if signal_channel is not in experiment."""
        self.op.signal_channel = 'NonExistent'
        with self.assertRaises(util.CytoflowOpError):
            self.op.apply(self.ex)

    def test_history_updated(self):
        """apply() must append to experiment history."""
        ex2 = self.op.apply(self.ex)
        self.assertGreater(len(ex2.history), len(self.ex.history))


class TestCalciumFluxView(unittest.TestCase):
    """Tests for CalciumFluxView.plot()"""

    def setUp(self):
        self.ex = _make_experiment()
        self.op = CalciumFluxOp(
            name='Ca',
            time_channel='Time',
            signal_channel='Fluo4_A',
            baseline_end=60.0,
            response_threshold=1.5,
            bin_width=10.0,
        )
        self.ex2 = self.op.apply(self.ex)

    def tearDown(self):
        import matplotlib.pyplot as plt
        plt.close('all')

    def test_default_view_runs(self):
        """default_view().plot() must complete without error."""
        self.op.default_view().plot(self.ex2)

    def test_plot_no_error_band(self):
        """plot() with show_error_band=False must run without error."""
        view = CalciumFluxView(op=self.op, show_error_band=False)
        view.plot(self.ex2)

    def test_plot_no_stimulation_line(self):
        """plot() with show_stimulation_line=False must run without error."""
        view = CalciumFluxView(op=self.op, show_stimulation_line=False)
        view.plot(self.ex2)

    def test_plot_raises_without_op(self):
        """plot() must raise if op is None."""
        view = CalciumFluxView()
        with self.assertRaises(util.CytoflowViewError):
            view.plot(self.ex2)

    def test_plot_raises_without_statistics(self):
        """plot() must raise if statistics are absent (operation not applied)."""
        view = CalciumFluxView(op=self.op)
        with self.assertRaises(util.CytoflowViewError):
            view.plot(self.ex)   # pass the unprocessed experiment


if __name__ == '__main__':
    unittest.main()
