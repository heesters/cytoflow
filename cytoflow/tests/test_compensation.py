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

import os, unittest, tempfile, warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import cytoflow as flow
import cytoflow.utility as util
from cytoflow.operations.compensation import (parse_spillover, fcs_spillover,
                                              read_spillover_csv,
                                              spillover_to_dataframe)
from cytoflow.operations.import_op import parse_tube, fcs_marker_names

ARIA = 'instruments/BD - FACS Aria II.fcs'
ARIA_CHANNELS = ['R780-A', 'B515-A', 'V450-A', 'V545-A', 'G560-A',
                 'G610-A', 'G660-A', 'G710-A', 'G780-A']

class TestCompensationSynthetic(unittest.TestCase):
    """Compensate data with a known spillover matrix, and get back the original"""

    def setUp(self):
        rng = np.random.default_rng(0)
        self.channels = ['A', 'B', 'C']

        # true signals
        self.x = rng.lognormal(mean = 5, sigma = 1.5, size = (2000, 3))

        # spillover: row = from, column = to
        self.s = np.array([[1.0, 0.20, 0.05],
                           [0.10, 1.0, 0.30],
                           [0.00, 0.15, 1.0]])

        # the cytometer measures x . S
        raw = self.x @ self.s

        self.ex = flow.Experiment()
        for c in self.channels:
            self.ex.add_channel(c)
            self.ex.metadata[c]['range'] = 1e6
            self.ex.metadata[c]['fcs_name'] = c

        self.ex.add_channel('FSC')
        self.ex.metadata['FSC']['range'] = 1e6
        self.ex.add_events(pd.DataFrame({'A' : raw[:, 0],
                                         'B' : raw[:, 1],
                                         'C' : raw[:, 2],
                                         'FSC' : rng.uniform(size = 2000)}),
                           {})
        self.ex.history.append(flow.ImportOp())

        self.df = pd.DataFrame(self.s,
                               index = self.channels,
                               columns = self.channels)

    def testRecover(self):
        op = flow.CompensationOp(spillover = self.df)
        ex2 = op.apply(self.ex)
        np.testing.assert_allclose(ex2.data[self.channels].values, self.x, rtol = 1e-9)

        # untouched channel
        pd.testing.assert_series_equal(ex2.data['FSC'], self.ex.data['FSC'])

        # original experiment is unchanged
        self.assertFalse(np.allclose(self.ex.data[self.channels].values, self.x))

        # column order is preserved
        self.assertEqual(list(ex2.data.columns), list(self.ex.data.columns))

    def testDictOfDicts(self):
        spill = {f : {t : self.s[i, j] for j, t in enumerate(self.channels)}
                 for i, f in enumerate(self.channels)}
        op = flow.CompensationOp(spillover = spill)
        ex2 = op.apply(self.ex)
        np.testing.assert_allclose(ex2.data[self.channels].values, self.x, rtol = 1e-9)

    def testMissingEntriesDefault(self):
        # diagonal defaults to 1, off-diagonal to 0
        op = flow.CompensationOp(spillover = {'A' : {'B' : 0.20, 'C' : 0.05},
                                              'B' : {'A' : 0.10, 'C' : 0.30},
                                              'C' : {'B' : 0.15}})
        ex2 = op.apply(self.ex)
        np.testing.assert_allclose(ex2.data[self.channels].values, self.x, rtol = 1e-9)

    def testDataFrameConversion(self):
        op = flow.CompensationOp(spillover = self.df)
        self.assertIsInstance(op.spillover['A'], dict)
        self.assertAlmostEqual(op.spillover['A']['B'], 0.20)
        self.assertAlmostEqual(op.spillover['B']['A'], 0.10)
        pd.testing.assert_frame_equal(spillover_to_dataframe(op.spillover), self.df)

    def testMetadata(self):
        ex2 = flow.CompensationOp(spillover = self.df).apply(self.ex)
        self.assertTrue(ex2.metadata['B']['compensated'])
        self.assertEqual(ex2.metadata['B']['spillover'], {'A' : 0.20, 'C' : 0.15})
        self.assertNotIn('compensated', ex2.metadata['FSC'])
        self.assertEqual(ex2.metadata['A']['range'], 1e6)
        self.assertNotIn('compensated', self.ex.metadata['A'])
        self.assertIsInstance(ex2.history[-1], flow.CompensationOp)

    def testDoubleCompensationWarns(self):
        ex2 = flow.CompensationOp(spillover = self.df).apply(self.ex)
        with self.assertWarns(util.CytoflowOpWarning):
            flow.CompensationOp(spillover = self.df).apply(ex2)

    def testErrors(self):
        with self.assertRaises(util.CytoflowOpError):
            flow.CompensationOp().apply(self.ex)

        with self.assertRaises(util.CytoflowOpError):
            flow.CompensationOp(spillover = {'A' : {'D' : 0.1}}).apply(self.ex)

        with self.assertRaises(util.CytoflowOpError):
            flow.CompensationOp(spillover = {'A' : {'A' : 1, 'B' : 1},
                                             'B' : {'A' : 1, 'B' : 1}}).apply(self.ex)

        with self.assertRaises(util.CytoflowOpError):
            flow.CompensationOp().estimate(self.ex)

    def testCSV(self):
        # FlowJo-style: header row of names, rows of values only
        fh, filename = tempfile.mkstemp(suffix = '.csv')
        os.close(fh)
        try:
            with open(filename, 'w') as f:
                f.write("A,B,C\n")
                for row in self.s:
                    f.write(",".join(str(v) for v in row) + "\n")
            pd.testing.assert_frame_equal(read_spillover_csv(filename), self.df)

            # with row names, "Comp-" prefixes and percentages
            with open(filename, 'w') as f:
                f.write(",Comp-A,Comp-B,Comp-C\n")
                for c, row in zip(self.channels, self.s):
                    f.write("Comp-" + c + "," + ",".join(str(v * 100) for v in row) + "\n")

            op = flow.CompensationOp(spillover_file = filename)
            op.estimate(self.ex)
            pd.testing.assert_frame_equal(spillover_to_dataframe(op.spillover), self.df)

            ex2 = op.apply(self.ex)
            np.testing.assert_allclose(ex2.data[self.channels].values, self.x, rtol = 1e-9)
        finally:
            os.unlink(filename)

    def testPlot(self):
        flow.CompensationOp(spillover = self.df).default_view().plot(self.ex)
        plt.close('all')

        with self.assertRaises(util.CytoflowViewError):
            flow.CompensationOp().default_view().plot(self.ex)


class TestParseSpillover(unittest.TestCase):

    def testRowMajor(self):
        df = parse_spillover("2,X,Y,1,0.25,0.1,1")
        self.assertEqual(list(df.index), ['X', 'Y'])
        # row 1 = spillover FROM X
        self.assertEqual(df.loc['X', 'Y'], 0.25)
        self.assertEqual(df.loc['Y', 'X'], 0.1)

    def testBadMatrix(self):
        with self.assertRaises(util.CytoflowError):
            parse_spillover("3,X,Y,Z,1,0")
        with self.assertRaises(util.CytoflowError):
            parse_spillover("X,Y,1,0")

    def testKeywords(self):
        for kw in ['$SPILLOVER', 'SPILL', '$SPILL', 'spill']:
            df = fcs_spillover({kw : "2,X,Y,1,0.25,0.1,1"})
            self.assertEqual(df.loc['X', 'Y'], 0.25)
        self.assertIsNone(fcs_spillover({'$PAR' : '2'}))

    def testParameterNumbers(self):
        # the Accuri C6 lists parameter numbers instead of names
        meta = {'$PAR' : '3', '$P1N' : 'FSC', '$P2N' : 'FL1', '$P3N' : 'FL2',
                '$SPILLOVER' : '2,2,3,1,0.25,0.1,1'}
        df = fcs_spillover(meta)
        self.assertEqual(list(df.index), ['FL1', 'FL2'])
        self.assertEqual(df.loc['FL1', 'FL2'], 0.25)


class TestCompensationFCS(unittest.TestCase):
    """Read the spillover matrix from a real FCS file (BD FACSAria II, SPILL keyword)"""

    def setUp(self):
        self.cwd = os.path.dirname(os.path.abspath(__file__)) + '/data/'
        self.file = self.cwd + ARIA

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            self.ex = flow.ImportOp(tubes = [flow.Tube(file = self.file)]).apply()

        meta, _ = parse_tube(self.file, metadata_only = True)
        self.expected = parse_spillover(meta['SPILL'])

    def testMetadataLocation(self):
        fcs_meta = self.ex.metadata['fcs_metadata']
        self.assertEqual(list(fcs_meta.keys()), [self.file])
        self.assertIn('SPILL', fcs_meta[self.file])

    def testEstimate(self):
        op = flow.CompensationOp()
        op.estimate(self.ex)

        self.assertEqual(set(op.spillover.keys()), set(ARIA_CHANNELS))

        # a value from the FCS file: row B515-A (2nd), column V545-A (4th)
        self.assertAlmostEqual(op.spillover['B515-A']['V545-A'], 0.030264132476262202)
        self.assertAlmostEqual(op.spillover['V450-A']['V545-A'],
                               self.expected.loc['V450-A', 'V545-A'])
        self.assertGreater(op.spillover['V450-A']['V545-A'], 0.5)

    def testApply(self):
        op = flow.CompensationOp()
        op.estimate(self.ex)
        ex2 = op.apply(self.ex)

        s = self.expected.values
        expected = self.ex.data[ARIA_CHANNELS].values @ np.linalg.inv(s)
        np.testing.assert_allclose(ex2.data[ARIA_CHANNELS].values, expected)

        # and the inverse: compensated . S == raw
        np.testing.assert_allclose(ex2.data[ARIA_CHANNELS].values @ s,
                                   self.ex.data[ARIA_CHANNELS].values,
                                   rtol = 1e-6, atol = 1e-6)

        # non-fluorescence channels are unchanged
        pd.testing.assert_series_equal(ex2.data['FSC-A'], self.ex.data['FSC-A'])

    def testRenamedChannels(self):
        names = {n : n for n, _ in fcs_marker_names(self.file)}
        names['B515-A'] = 'FITC'
        names['V545-A'] = 'V500'

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ex = flow.ImportOp(tubes = [flow.Tube(file = self.file)],
                               channels = names).apply()

        op = flow.CompensationOp()
        op.estimate(ex)
        self.assertIn('FITC', op.spillover)
        self.assertNotIn('B515-A', op.spillover)
        self.assertAlmostEqual(op.spillover['FITC']['V500'], 0.030264132476262202)

        ex2 = op.apply(ex)
        ex_ref = flow.CompensationOp(spillover = self.expected).apply(self.ex)
        np.testing.assert_allclose(ex2.data['FITC'], ex_ref.data['B515-A'])
        np.testing.assert_allclose(ex2.data['V500'], ex_ref.data['V545-A'])

    def testMarkerNames(self):
        # channels named by $PnS
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ex = flow.ImportOp(tubes = [flow.Tube(file = self.file)],
                               name_metadata = '$PnS').apply()

        meta, _ = parse_tube(self.file, metadata_only = True)
        pns = {meta['$P{}N'.format(i)] : meta.get('$P{}S'.format(i))
               for i in range(1, int(meta['$PAR']) + 1)}
        fitc = [c for c in ex.channels if ex.metadata[c]['fcs_name'] == pns['B515-A']]
        v500 = [c for c in ex.channels if ex.metadata[c]['fcs_name'] == pns['V545-A']]
        self.assertEqual(len(fitc), 1)
        self.assertEqual(len(v500), 1)
        self.assertNotEqual(fitc[0], 'B515-A')
        
        op = flow.CompensationOp()
        op.estimate(ex)
        self.assertAlmostEqual(op.spillover[fitc[0]][v500[0]], 0.030264132476262202)

    def testSubsetOfChannels(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ex = flow.ImportOp(tubes = [flow.Tube(file = self.file)],
                               channels = {'FSC-A' : 'FSC-A',
                                           'B515-A' : 'B515-A',
                                           'V545-A' : 'V545-A'}).apply()
        op = flow.CompensationOp()
        with self.assertWarns(util.CytoflowOpWarning):
            op.estimate(ex)
        self.assertEqual(set(op.spillover.keys()), {'B515-A', 'V545-A'})

    def testNoSpillover(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ex = flow.ImportOp(tubes = [flow.Tube(file = self.cwd + "Plate01/RFP_Well_A3.fcs")]).apply()
        with self.assertRaisesRegex(util.CytoflowOpError, "spillover matrix"):
            flow.CompensationOp().estimate(ex)

    def testPlot(self):
        op = flow.CompensationOp()
        op.estimate(self.ex)
        op.default_view().plot(self.ex)
        plt.close('all')

if __name__ == "__main__":
    unittest.main()
