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

'''
Compensation
------------

Compensate fluorescence channels using the spillover (compensation) matrix
that the cytometer stored in the FCS file.

Most cytometers' acquisition software saves the compensation matrix that was
set up during acquisition in the FCS file (in the ``$SPILLOVER``, ``SPILL``
or ``$SPILL`` keyword).  This module reads that matrix from the first tube in
the experiment and uses it to compensate the data, the same way FlowJo or
FCS Express do when you apply the "acquisition-defined" matrix.

Click **Read spillover matrix** to read the matrix, check it in the
diagnostic plot (a heat map of the spillover values), and edit any values
you'd like to change in the table.

If you'd rather use a matrix from somewhere else -- for example, one you
exported from FlowJo -- choose a **CSV file** before clicking **Read
spillover matrix**.  The first row of the file must be the channel names;
each following row is the spillover *from* one channel, in the same order
as the columns, and may start with that channel's name.

The compensated data are :math:`x = r S^{-1}`, where :math:`r` is the row of
raw values for an event and :math:`S_{ij}` is the fraction of channel
:math:`i`'s signal that spills into channel :math:`j`.

.. note::

    Use either this module or **Linear Compensation** (which estimates the
    matrix from single-color controls), not both.  Compensation should be one
    of the first operations, before (for example) autofluorescence
    correction.

.. object:: CSV file (optional)

    Read the spillover matrix from this CSV file instead of the FCS file.

.. object:: Read spillover matrix

    Read the matrix from the FCS file (or the CSV file).  This replaces any
    values you edited.

.. object:: Spillover matrix

    The spillover values: the fraction of the **From** channel's signal that
    appears in the **To** channel (eg, 0.05 is 5%.)  Click a value to edit it.

.. plot::
   :include-source: False

    import cytoflow as flow
    import_op = flow.ImportOp()
    import_op.tubes = [flow.Tube(file = "instruments/BD - FACS Aria II.fcs")]
    ex = import_op.apply()

    comp_op = flow.CompensationOp()
    comp_op.estimate(ex)
    comp_op.default_view().plot(ex)

    ex2 = comp_op.apply(ex)
'''

from traits.api import provides, List
from traitsui.api import (View, Item, VGroup, ButtonEditor, FileEditor)

from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..view_plugins import ViewHandler
from ..editors import ColorTextEditor, SpilloverMatrixEditor
from ..workflow.operations import CompensationWorkflowOp, CompensationWorkflowView

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, shared_op_traits_view, PluginHelpMixin


class CompensationHandler(OpHandler):
    operation_traits_view = \
        View(
            VGroup(
                Item('spillover_file',
                     label = "CSV file\n(optional)",
                     editor = FileEditor(dialog_style = 'open',
                                         filter = ["*.csv", "*.txt", "*"])),
                Item('do_estimate',
                     editor = ButtonEditor(value = True,
                                           label = "Read spillover matrix"),
                     show_label = False),
                VGroup(
                    Item('spillover_list',
                         editor = SpilloverMatrixEditor(),
                         show_label = False,
                         height = 220),
                    label = "Spillover Matrix",
                    show_border = False),
                shared_op_traits_view))


class CompensationViewHandler(ViewHandler):
    view_traits_view = \
        View(Item('context.view_warning',
                  resizable = True,
                  visible_when = 'context.view_warning',
                  editor = ColorTextEditor(foreground_color = "#000000",
                                          background_color = "#ffff99")),
             Item('context.view_error',
                  resizable = True,
                  visible_when = 'context.view_error',
                  editor = ColorTextEditor(foreground_color = "#000000",
                                           background_color = "#ff9191")))

    view_params_view = View()


@provides(IOperationPlugin)
class CompensationPlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.compensation'
    operation_id = 'cytoflow.operations.compensation'
    view_id = 'cytoflow.view.compensationdiagnostic'

    name = "Compensation (from FCS)"
    short_name = "Compensation"
    menu_group = "Preprocessing"

    def get_operation(self):
        return CompensationWorkflowOp()

    def get_handler(self, model, context):
        if isinstance(model, CompensationWorkflowOp):
            return CompensationHandler(model = model, context = context)
        elif isinstance(model, CompensationWorkflowView):
            return CompensationViewHandler(model = model, context = context)

    def get_icon(self):
        return ImageResource('compensation')

    plugin = List(contributes_to = OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
