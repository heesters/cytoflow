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
Calcium Flux
------------

Normalise a calcium-flux fluorescence channel to its pre-stimulation
baseline (F/F₀) and derive kinetics statistics from binned time windows.

.. object:: Name

    The output prefix for the new channel, conditions and statistics.

.. object:: Time Channel

    The channel recording time (e.g. ``Time``).

.. object:: Signal Channel

    The fluorescence channel to analyse (e.g. ``Fluo4_A``).

.. object:: Ratio Denominator (optional)

    For ratiometric dyes (e.g. Indo-1) set this to the second channel.
    When non-empty the signal is computed as
    *signal_channel / ratio_denominator* before normalisation.
    Leave empty for single-channel dyes such as Fluo-4.

.. object:: Baseline End (s)

    Events with time ≤ this value (seconds) form the pre-stimulation
    baseline used to compute F₀.

.. object:: Response Threshold

    F/F₀ value above which a cell is classified as "Responding".

.. object:: Bin Width (s)

    Width of each time bin in seconds.

'''

from traits.api import provides, List
from traitsui.api import View, Item, EnumEditor, VGroup, TextEditor
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..view_plugins import ViewHandler
from ..editors import SubsetListEditor, ColorTextEditor, InstanceHandlerEditor
from ..workflow.operations.calcium_flux import (CalciumFluxWorkflowOp,
                                                CalciumFluxWorkflowView)
from ..subset_controllers import subset_handler_factory

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, shared_op_traits_view, PluginHelpMixin


class CalciumFluxHandler(OpHandler):
    operation_traits_view = \
        View(Item('name',
                  editor=TextEditor(auto_set=False,
                                    placeholder="None"),
                  label="Name"),
             Item('time_channel',
                  editor=EnumEditor(name='context_handler.previous_channels'),
                  label="Time Channel"),
             Item('signal_channel',
                  editor=EnumEditor(name='context_handler.previous_channels'),
                  label="Signal Channel"),
             Item('ratio_denominator',
                  editor=TextEditor(auto_set=False,
                                    placeholder="(optional)"),
                  label="Ratio Denom."),
             Item('baseline_end',
                  editor=TextEditor(auto_set=False,
                                    evaluate=float),
                  label="Baseline End (s)"),
             Item('response_threshold',
                  editor=TextEditor(auto_set=False,
                                    evaluate=float),
                  label="Resp. Threshold"),
             Item('bin_width',
                  editor=TextEditor(auto_set=False,
                                    evaluate=float),
                  label="Bin Width (s)"),
             shared_op_traits_view)


class CalciumFluxViewHandler(ViewHandler):
    view_traits_view = \
        View(Item('context.view_warning',
                  resizable=True,
                  visible_when='context.view_warning',
                  editor=ColorTextEditor(foreground_color="#000000",
                                         background_color="#ffff99")),
             Item('context.view_error',
                  resizable=True,
                  visible_when='context.view_error',
                  editor=ColorTextEditor(foreground_color="#000000",
                                          background_color="#ff9191")))

    view_params_view = View()


@provides(IOperationPlugin)
class CalciumFluxPlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.calcium_flux'
    operation_id = 'cytoflow.operations.calcium_flux'
    view_id = 'cytoflow.views.calcium_flux'

    name = "Calcium Flux"
    short_name = "Calcium\nFlux"
    menu_group = "TOOLS"

    def get_operation(self):
        return CalciumFluxWorkflowOp()

    def get_handler(self, model, context):
        if isinstance(model, CalciumFluxWorkflowOp):
            return CalciumFluxHandler(model=model, context=context)
        elif isinstance(model, CalciumFluxWorkflowView):
            return CalciumFluxViewHandler(model=model, context=context)

    def get_icon(self):
        return ImageResource('calcium_flux')

    plugin = List(contributes_to=OP_PLUGIN_EXT)

    def _plugin_default(self):
        return [self]
