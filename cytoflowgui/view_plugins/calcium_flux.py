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
Calcium Flux View
-----------------

Plots mean F/F₀ vs time for a calcium-flux experiment.

This view is associated with **Calcium Flux** operation and is available via
``op.default_view()``.

.. object:: Hue Facet

    Colour separate lines by this condition.

.. object:: Row Facet

    Arrange plots in multiple rows, one per unique value of this condition.

.. object:: Column Facet

    Arrange plots in multiple columns, one per unique value of this condition.

.. object:: Show Error Band

    Shade ±SEM around the mean line.

.. object:: Show Stimulation Line

    Draw a vertical dashed line at the end of the baseline window.

"""

from traits.api import provides, List
from traitsui.api import View, Item, EnumEditor, VGroup, CheckBoxEditor
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..editors import ColorTextEditor
from ..workflow.views.calcium_flux import CalciumFluxWorkflowView

from .i_view_plugin import IViewPlugin, VIEW_PLUGIN_EXT
from .view_plugin_base import ViewHandler, PluginHelpMixin


class CalciumFluxViewHandler(ViewHandler):
    view_traits_view = \
        View(VGroup(
             Item('huefacet',
                  editor=EnumEditor(name='context_handler.previous_conditions_names'),
                  label="Hue Facet"),
             Item('rowfacet',
                  editor=EnumEditor(name='context_handler.previous_conditions_names'),
                  label="Row Facet"),
             Item('colfacet',
                  editor=EnumEditor(name='context_handler.previous_conditions_names'),
                  label="Col Facet"),
             Item('show_error_band',
                  label="Error Band"),
             Item('show_stimulation_line',
                  label="Stim. Line"),
             label="Calcium Flux Plot",
             show_border=False),
             Item('context.view_warning',
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


@provides(IViewPlugin)
class CalciumFluxViewPlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.view_plugins.calcium_flux'
    view_id = 'cytoflow.views.calcium_flux'

    name = short_name = "Calcium Flux"

    def get_view(self):
        return CalciumFluxWorkflowView()

    def get_handler(self, model, context):
        if isinstance(model, CalciumFluxWorkflowView):
            return CalciumFluxViewHandler(model=model, context=context)

    def get_icon(self):
        return ImageResource('calcium_flux')

    plugin = List(contributes_to=VIEW_PLUGIN_EXT)

    def _plugin_default(self):
        return [self]
