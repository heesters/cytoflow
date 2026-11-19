#!/usr/bin/env python3.8
# coding: latin-1

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
Backgating
----------

Show all events in gray in the background, with a selected gate population
highlighted in color in the foreground.  This is the standard way to verify
gate accuracy and visualise what cells you are selecting.

.. object:: X Channel, Y Channel

    The channels to plot on the X and Y axes.

.. object:: X Scale, Y Scale

    How to scale the axes.

.. object:: Gate (highlight)

    The boolean gate condition whose events are drawn on top in color.
    All other events are drawn in gray behind them.

.. object:: Subset

    Plot only a subset of the data.

"""

from traits.api import provides, List
from traitsui.api import (View, Item, EnumEditor, VGroup)
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..workflow.views import ScatterplotWorkflowView
from ..editors import SubsetListEditor, ExtendableEnumEditor, InstanceHandlerEditor
from ..subset_controllers import subset_handler_factory

from .i_view_plugin import IViewPlugin, VIEW_PLUGIN_EXT
from .view_plugin_base import ViewHandler, PluginHelpMixin


class BackgateHandler(ViewHandler):

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
                    Item('highlight',
                         editor=ExtendableEnumEditor(
                             name='context_handler.bool_conditions_names',
                             extra_items={"None": ""}),
                         label="Gate\n(highlight)",
                         tooltip="Events where this gate is True are drawn on top "
                                 "in color; all other events are drawn in gray."),
                    label="Backgating",
                    show_border=False),
             VGroup(Item('subset_list',
                         show_label=False,
                         editor=SubsetListEditor(
                             conditions="context_handler.conditions",
                             editor=InstanceHandlerEditor(
                                 view='subset_view',
                                 handler_factory=subset_handler_factory),
                             mutable=False)),
                    label="Subset",
                    show_border=False,
                    show_labels=False)))


@provides(IViewPlugin)
class BackgatingPlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.view.backgate'
    view_id = 'cytoflow.view.scatterplot'
    name = "Backgating"
    short_name = "Backgate"

    def get_view(self):
        return ScatterplotWorkflowView()

    def get_handler(self, model, context):
        return BackgateHandler(model=model, context=context)

    def get_icon(self):
        return ImageResource('scatterplot')

    plugin = List(contributes_to=VIEW_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
