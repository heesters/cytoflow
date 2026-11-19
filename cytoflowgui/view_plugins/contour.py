#!/usr/bin/env python3.8
# coding: latin-1

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
Contour Plot
------------

A 2-D density contour map using fast ``numpy.histogram2d`` binning.  Unlike
the 2D Kernel Density Estimate, this view handles very large datasets (> 100 k
events) in a fraction of a second.

.. object:: X Channel, Y Channel

    The channels to plot on the X and Y axes.

.. object:: X Scale, Y Scale

    How to scale the axes.

.. object:: Filled contours

    Draw filled contour regions instead of iso-lines.

.. object:: Show scatter

    Draw individual events as a faint grey scatter underneath the contours.

.. object:: Num levels

    How many contour levels to draw.

.. object:: Grid size

    Resolution of the underlying 2-D histogram (pixels per axis).

.. object:: Horizontal Facet / Vertical Facet / Color Facet / Tab Facet

    Standard facetting controls.

.. object:: Subset

    Plot only a subset of the data.
"""

from traits.api import provides, List, observe
from traitsui.api import View, Item, EnumEditor, VGroup, TextEditor, Controller
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..workflow.views import ContourWorkflowView, ContourPlotParams
from ..editors import SubsetListEditor, ColorTextEditor, ExtendableEnumEditor, InstanceHandlerEditor
from ..subset_controllers import subset_handler_factory

from .i_view_plugin import IViewPlugin, VIEW_PLUGIN_EXT
from .view_plugin_base import ViewHandler, PluginHelpMixin, Data2DPlotParamsView


class ContourParamsHandler(Controller):

    view_params_view = \
        View(Item('filled',
                  label="Filled\ncontours"),
             Item('show_scatter',
                  label="Show\nscatter"),
             Item('n_levels',
                  editor=TextEditor(auto_set=False),
                  label="Num\nlevels"),
             Item('gridsize',
                  editor=TextEditor(auto_set=False),
                  label="Grid\nsize"),
             Item('alpha',
                  editor=TextEditor(auto_set=False)),
             Data2DPlotParamsView.content)


class ContourHandler(ViewHandler):

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
                    Item('xfacet',
                         editor=ExtendableEnumEditor(name='context_handler.conditions_names',
                                                     extra_items={"None": ""}),
                         label="Horizontal\nFacet"),
                    Item('yfacet',
                         editor=ExtendableEnumEditor(name='context_handler.conditions_names',
                                                     extra_items={"None": ""}),
                         label="Vertical\nFacet"),
                    Item('huefacet',
                         editor=ExtendableEnumEditor(name='context_handler.conditions_names',
                                                     extra_items={"None": ""}),
                         label="Color\nFacet"),
                    Item('plotfacet',
                         editor=ExtendableEnumEditor(name='context_handler.conditions_names',
                                                     extra_items={"None": ""}),
                         label="Tab\nFacet"),
                    label="Contour Plot",
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
                                               handler_factory=ContourParamsHandler),
                  style='custom',
                  show_label=False))


@provides(IViewPlugin)
class ContourPlugin(Plugin, PluginHelpMixin):
    id = 'cytoflowgui.view.contour'
    view_id = 'cytoflow.view.contour'
    name = "Contour Plot"
    short_name = "Contour"

    def get_view(self):
        return ContourWorkflowView()

    def get_handler(self, model, context):
        if isinstance(model, ContourWorkflowView):
            return ContourHandler(model=model, context=context)
        elif isinstance(model, ContourPlotParams):
            return ContourParamsHandler(model=model, context=context)

    def get_icon(self):
        return ImageResource('kde_2d')

    plugin = List(contributes_to=VIEW_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
