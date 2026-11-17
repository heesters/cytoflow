#!/usr/bin/env python3.8

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
cytoflowgui.report_export
-------------------------

The logic behind File > Export report: which plots a workflow has, and the
dialog's model.  (The report itself is rendered in the remote process; see
`LocalWorkflow.export_report` and `cytoflow.utility.report`.)
"""

from traits.api import HasTraits, Str, List, Int, Enum, Range
from traitsui.api import View, Item, VGroup, HGroup, CheckListEditor, TextEditor, Label

from cytoflow.utility.report import PAGE_SIZES


def list_plots(workflow):
    """
    All the plots in a workflow.
    
    Parameters
    ----------
    workflow : list of WorkflowItem
    
    Returns
    -------
    list of (WorkflowItem, IWorkflowView, str)
        Each plot, with a label for it: the operation's position and name,
        and the view's title (if it has one) or type.
    """
    plots = []
    for i, wi in enumerate(workflow):
        op = wi.operation
        op_name = getattr(op, 'name', '') or getattr(op, 'friendly_id', '') or "Operation"
        for view in wi.views:
            params = getattr(view, 'plot_params', None)
            title = getattr(params, 'title', '') if params is not None else ''
            what = title or getattr(view, 'friendly_id', '') or view.id
            current = getattr(view, 'current_plot', None)
            if current not in (None, ''):
                what += " ({})".format(current)
            plots.append((wi, view, "{}. {}: {}".format(i + 1, op_name, what)))
    return plots


class ReportDialogModel(HasTraits):
    """What the user chooses in the "Export report" dialog"""
    
    # the labels of all the plots, and of the chosen ones
    labels = List(Str)
    selected = List(Str)
    
    title = Str
    rows = Range(1, 6, 2)
    columns = Range(1, 6, 2)
    page_size = Enum(*PAGE_SIZES.keys())
    dpi = Range(50, 600, 300)
    
    traits_view = View(
        VGroup(Label("Choose the plots to put in the report.  Each plot is drawn\n"
                     "as it is in Cytoflow (the same gates, scales and settings.)"),
               Item('selected',
                    editor = CheckListEditor(name = 'labels', cols = 1),
                    style = 'custom',
                    show_label = False),
               label = "Plots",
               show_border = True),
        VGroup(Item('title',
                    editor = TextEditor(auto_set = False, placeholder = "None")),
               HGroup(Item('rows'), Item('columns', label = "Columns")),
               Item('page_size', label = "Page size\n(PDF only)"),
               Item('dpi', label = "Resolution\n(dpi)"),
               Label("PDF with 1 row and 1 column: one plot per page, as a vector\n"
                     "graphic with editable text.  Otherwise, and in PowerPoint,\n"
                     "the plots are images."),
               label = "Layout",
               show_border = True),
        title = "Export report",
        buttons = ['OK', 'Cancel'],
        resizable = True,
        width = 480)
    
    def options(self, filename):
        """The options for `LocalWorkflow.export_report`"""
        opts = dict(title = self.title,
                    layout = (self.rows, self.columns),
                    dpi = self.dpi)
        if filename.lower().endswith(".pdf"):
            opts['page_size'] = self.page_size
        return opts
