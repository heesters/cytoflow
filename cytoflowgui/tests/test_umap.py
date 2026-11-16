#!/usr/bin/env python3.8
# coding: latin-1

# (c) Massachusetts Institute of Technology 2015-2018
# (c) Brian Teague 2018-2025
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

'''
Tests for the UMAP GUI operation (modelled on test_pca.py)
'''

import os, unittest, tempfile
import pandas as pd

# needed for testing lambdas
from cytoflow import geom_mean, geom_sd  # @UnusedImport

from cytoflowgui.tests.test_base import ImportedDataTest
from cytoflowgui.workflow.workflow_item import WorkflowItem
from cytoflowgui.workflow.operations import UMAPWorkflowOp, UMAPChannel
from cytoflowgui.workflow.subset import CategorySubset, RangeSubset
from cytoflowgui.workflow.serialization import load_yaml, save_yaml

class TestUMAP(ImportedDataTest):
    
    # the first UMAP estimate in a fresh process is slow, because numba has
    # to compile umap-learn's (and pynndescent's) code. 
    ESTIMATE_TIMEOUT = 300
    
    def setUp(self):
        super().setUp()

        self.addTypeEqualityFunc(UMAPWorkflowOp, 'assertHasTraitsEqual')
        self.addTypeEqualityFunc(UMAPChannel, 'assertHasTraitsEqual')

        self.op = op = UMAPWorkflowOp()
        
        op.name = "UMAP"
        op.channels_list = [UMAPChannel(channel = "V2-A", scale = "log"),
                            UMAPChannel(channel = "V2-H", scale = "log"),
                            UMAPChannel(channel = "Y2-A", scale = "log"),
                            UMAPChannel(channel = "Y2-H", scale = "log")]
        op.n_components = 2
        op.sample = 0.01
        op.random_state = 0
                
        op.subset_list.append(CategorySubset(name = "Well",
                                             values = ['A', 'B']))
        op.subset_list.append(RangeSubset(name = "Dox",
                                          values = [0.0, 10.0, 100.0]))
        
        self.wi = wi = WorkflowItem(operation = op,
                                    status = 'waiting',
                                    view_error = "Not yet plotted")
        
        self.workflow.workflow.append(wi)
        self.workflow.selected = wi
        
    def _estimate(self):
        if self.wi.status != 'waiting':
            self.workflow.wi_sync(self.wi, 'status', 'waiting')
        self.op.do_estimate = True
        self.workflow.wi_waitfor(self.wi, 'status', 'valid', 
                                 timeout = self.ESTIMATE_TIMEOUT)
        self.assertTrue(self.workflow.remote_eval("self.workflow[-1].result is not None"))

    def testEstimate(self):
        self._estimate()
        self.assertTrue(self.workflow.remote_eval("'UMAP_1' in self.workflow[-1].result.channels"))
        self.assertTrue(self.workflow.remote_eval("'UMAP_2' in self.workflow[-1].result.channels"))
        self.assertFalse(self.workflow.remote_eval("'UMAP_3' in self.workflow[-1].result.channels"))
        
    def testChangeParameters(self):
        self._estimate()
        
        # changing any estimate parameter should invalidate the estimate
        changes = [lambda op: setattr(op, 'n_neighbors', 10),
                   lambda op: setattr(op, 'min_dist', 0.5),
                   lambda op: setattr(op, 'metric', 'cosine'),
                   lambda op: setattr(op, 'random_state', 1),
                   lambda op: setattr(op, 'by', ['Dox']),
                   lambda op: setattr(op.channels_list[0], 'scale', 'logicle'),
                   lambda op: op.channels_list.pop(),
                   lambda op: op.channels_list.append(UMAPChannel(channel = "B1-A", scale = "log")),
                   lambda op: setattr(op.subset_list[0], 'selected', ['A'])]
        
        for change in changes:
            self.workflow.wi_sync(self.wi, 'status', 'waiting')
            change(self.op)
            self.workflow.wi_waitfor(self.wi, 'status', 'invalid')
            self.assertTrue(self.workflow.remote_eval("self.workflow[-1].result is None"))
            self._estimate()
        
    def testComponents(self):
        self.op.n_components = 3
        self._estimate()
        self.assertTrue(self.workflow.remote_eval("'UMAP_3' in self.workflow[-1].result.channels"))

    def testSerializeOp(self):
        fh, filename = tempfile.mkstemp()
        try:
            os.close(fh)

            save_yaml(self.op, filename)
            new_op = load_yaml(filename)
        finally:
            os.unlink(filename)

        self.maxDiff = None

        self.assertEqual(self.op, new_op)
                      
    def testSerializeWorkflowItem(self):
        fh, filename = tempfile.mkstemp()
        try:
            os.close(fh)
             
            save_yaml(self.wi, filename)
            new_wi = load_yaml(filename)
             
        finally:
            os.unlink(filename)
             
        self.maxDiff = None
        
        self.assertEqual(self.wi, new_wi)
                                     
    def testNotebook(self):
        self._estimate()
        
        code = "import cytoflow as flow\n"
        for i, wi in enumerate(self.workflow.workflow):
            code = code + wi.operation.get_notebook_code(i)
            
            for view in wi.views:
                code = code + view.get_notebook_code(i)
         
        # use an explicit namespace: exec() can't reliably add variables
        # to a function's locals()
        ns = dict(globals())
        exec(code, ns)
            
        nb_data = ns['ex_{}'.format(len(self.workflow.workflow) - 1)].data
        remote_data = self.workflow.remote_eval("self.workflow[-1].result.data")
        
        pd.testing.assert_frame_equal(nb_data, remote_data)

if __name__ == "__main__":
#     import sys;sys.argv = ['', 'TestUMAP.testSerializeOp']
    unittest.main()
