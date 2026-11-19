#!/usr/bin/env python3.8
# coding: latin-1

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
Overlay Histogram
-----------------

Overlay histograms of multiple gated populations on one panel using the
Okabe-Ito colorblind-safe palette.  An "All events" background trace is always
shown in grey.

.. object:: Channel

    The channel to histogram.

.. object:: Scale

    Axis scale.

.. object:: Populations

    Boolean conditions (gate names) to overlay.  Tick the gates whose events
    you want to show as separate colored traces.

.. object:: Density

    Normalise each trace to a probability density so populations of different
    sizes are visually comparable.

.. object:: Subset

    Plot only a subset of the data.
"""

from traits.api import provides, List, observe
from traitsui.api import (View, Item, EnumEditor, VGroup, CheckListEditor,
                           ColorTextEditor, Controller)
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..workflow.views import OverlayHistogramWorkflowView, OverlayHistogramPlotParams
from ..editors import SubsetListEditor, ColorTextEditor, InstanceHandlerEditor
from ..subset_controllers import subset_handler_factory

from .i_view_plugin import IViewPlugin, VIEW_PLUGIN_EXT
from .view_plugin_base import ViewHandler, PluginHelpMixin


class OverlayHistogramParamsHandler(Controller):

    view_params_view = \
        View(Item('bins'),
             Item('alpha'),
             Item('density'))


class OverlayHistogramHandler(ViewHandler):

    @observe('model.channel', post_init=True)
    def _channel_changed(self, _):
        self._apply_default_scale('channel', 'scale')

    view_traits_view = \
        View(VGroup(
             VGroup(Item('channel',
                         editor=EnumEditor(name='context_handler.channels'),
                         label="Channel"),
                    Item('scale',
                         label="Scale"),
                    Item('populations',
                         editor=CheckListEditor(cols=2,
                                               name='context_handler.bool_conditions_names'),
                         label="Populations",
                         style='custom',
                         tooltip="Tick the gate populations to overlay"),
                    label="Overlay Histogram",
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
                                               handler_factory=OverlayHistogramParamsHandler),
                  style='custom',
                  show_label=False))


@provides(IViewPlugin)
class OverlayHistogramPlugin(Plugin, PluginHelpMixin):
    id = 'cytoflowgui.view.overlay_histogram'
    view_id = 'cytoflow.view.overlay_histogram'
    name = "Overlay Histogram"
    short_name = "Overlay"

    def get_view(self):
        return OverlayHistogramWorkflowView()

    def get_handler(self, model, context):
        if isinstance(model, OverlayHistogramWorkflowView):
            return OverlayHistogramHandler(model=model, context=context)
        elif isinstance(model, OverlayHistogramPlotParams):
            return OverlayHistogramParamsHandler(model=model, context=context)

    def get_icon(self):
        return ImageResource('histogram')

    plugin = List(contributes_to=VIEW_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
