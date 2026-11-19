#!/usr/bin/env python3.8
# coding: latin-1

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
Dot Plot
--------

A scatter plot where every event is colored by which boolean gate it belongs
to.  Events belonging to no listed gate are shown in light grey.  Later gates
in the list draw on top of earlier ones.

This is useful for visualising the composition of a scatter plot: you can see
at a glance where your CD4+ cells sit relative to your CD8+ cells on
FSC/SSC, or where live/dead cells cluster.

.. object:: X Channel, Y Channel

    The channels to plot on the X and Y axes.

.. object:: X Scale, Y Scale

    How to scale the axes.

.. object:: Populations

    Tick the boolean gate conditions to color.  Order matters: later entries
    are drawn on top.

.. object:: Horizontal Facet / Vertical Facet / Tab Facet

    Standard facetting controls.

.. object:: Subset

    Plot only a subset of the data.
"""

from traits.api import provides, List, observe
from traitsui.api import (View, Item, EnumEditor, VGroup, CheckListEditor,
                           TextEditor, Controller)
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..workflow.views import DotPlotWorkflowView, DotPlotParams
from ..editors import SubsetListEditor, ColorTextEditor, ExtendableEnumEditor, InstanceHandlerEditor
from ..subset_controllers import subset_handler_factory

from .i_view_plugin import IViewPlugin, VIEW_PLUGIN_EXT
from .view_plugin_base import ViewHandler, PluginHelpMixin, Data2DPlotParamsView


class DotPlotParamsHandler(Controller):

    view_params_view = \
        View(Item('alpha',
                  editor=TextEditor(auto_set=False)),
             Item('size',
                  editor=TextEditor(auto_set=False)),
             Data2DPlotParamsView.content)


class DotPlotHandler(ViewHandler):

    @observe('model.xchannel', post_init=True)
    def _xchannel_changed(self, _):
        self._apply_default_scale('xchannel', 'xscale')

    @observe('model.ychannel', post_init=True)
    def _ychannel_changed(self, _):
        self._apply_default_scale('ychannel', 'yscale')

    view_traits_view = \
        View(VGroup(
             VGroup(Item('xchannel',
                         editor=EnumEditor(name='context_handler.channels'),
                         label="X Channel"),
                    Item('xscale',
                         label="X Scale"),
                    Item('ychannel',
                         editor=EnumEditor(name='context_handler.channels'),
                         label="Y Channel"),
                    Item('yscale',
                         label="Y Scale"),
                    Item('populations',
                         editor=CheckListEditor(cols=2,
                                               name='context_handler.bool_conditions_names'),
                         label="Populations",
                         style='custom',
                         tooltip="Color events by which gate they belong to"),
                    Item('xfacet',
                         editor=ExtendableEnumEditor(name='context_handler.conditions_names',
                                                     extra_items={"None": ""}),
                         label="Horizontal\nFacet"),
                    Item('yfacet',
                         editor=ExtendableEnumEditor(name='context_handler.conditions_names',
                                                     extra_items={"None": ""}),
                         label="Vertical\nFacet"),
                    Item('plotfacet',
                         editor=ExtendableEnumEditor(name='context_handler.conditions_names',
                                                     extra_items={"None": ""}),
                         label="Tab\nFacet"),
                    label="Dot Plot",
                    show_border=False),
             VGroup(Item('subset_list',
                         show_label=False,
                         editor=SubsetListEditor(conditions="context_handler.conditions",
                                                 editor=InstanceHandlerEditor(view='subset_view',
                                                                              handler_factory=subset_handler_factory),
                                                 mutable=False)),
                    label="Subset",
                    show_border=False,
                    show_labels=False),
             Item('context.view_warning',
                  resizable=True,
                  visible_when='context.view_warning',
                  editor=ColorTextEditor(foreground_color="#000000",
                                         background_color="#ffff99")),
             Item('context.view_error',
                  resizable=True,
                  visible_when='context.view_error',
                  editor=ColorTextEditor(foreground_color="#000000",
                                          background_color="#ff9191"))))

    view_params_view = \
        View(Item('plot_params',
                  editor=InstanceHandlerEditor(view='view_params_view',
                                               handler_factory=DotPlotParamsHandler),
                  style='custom',
                  show_label=False))


@provides(IViewPlugin)
class DotPlotPlugin(Plugin, PluginHelpMixin):
    id = 'cytoflowgui.view.dot_plot'
    view_id = 'cytoflow.view.dot_plot'
    name = "Dot Plot"
    short_name = "Dot Plot"

    def get_view(self):
        return DotPlotWorkflowView()

    def get_handler(self, model, context):
        if isinstance(model, DotPlotWorkflowView):
            return DotPlotHandler(model=model, context=context)
        elif isinstance(model, DotPlotParams):
            return DotPlotParamsHandler(model=model, context=context)

    def get_icon(self):
        return ImageResource('scatterplot')

    plugin = List(contributes_to=VIEW_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
