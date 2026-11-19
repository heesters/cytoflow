#!/usr/bin/env python3.8

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
Figure Layout
-------------

Arrange multiple plots in a rows × cols grid.

Add panels using the **+ Panel** button, choose the view type for each cell,
and pick the channels.  All panels use the data from the current workflow step.
When the grid has more cells than panels, the extra cells are left blank.

.. object:: Rows / Columns
    Grid dimensions.

.. object:: Figure title
    Optional title printed above the grid.

.. object:: Panel width / height (inches)
    Size of each panel cell in the rendered figure.

.. object:: Panels
    One per cell (row-major).  For each:

    * **Type** — which kind of plot.
    * **X channel / Y channel** — the channels.
    * **X scale / Y scale** — axis scales.
    * **Color facet** — color by condition.
    * **Subset** — restrict to a subset.
    * **Panel title** — override the auto-generated title.
"""

from traits.api import provides, List, Instance, Button, observe
from traitsui.api import (View, Item, VGroup, HGroup, TextEditor,
                          EnumEditor, ButtonEditor, Controller, Label, Spring)
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from cytoflow.views.figure_layout import FigurePanel

from ..workflow.views import FigureLayoutWorkflowView
from ..editors import VerticalListEditor, InstanceHandlerEditor, ColorTextEditor

from .i_view_plugin import IViewPlugin, VIEW_PLUGIN_EXT
from .view_plugin_base import ViewHandler, PluginHelpMixin

PANEL_VIEW_TYPES = ['Histogram', 'Scatterplot', 'Density', 'KDE 1D',
                    'KDE 2D', 'Histogram 2D', 'Bar Chart']


class FigurePanelHandler(Controller):
    """Handler for one FigurePanel row in the list editor."""
    default_view = View(
        HGroup(
            Item('view_type',
                 editor=EnumEditor(values=PANEL_VIEW_TYPES),
                 show_label=False, width=100,
                 tooltip="Plot type for this panel"),
            Item('channel_x',
                 editor=TextEditor(auto_set=False, placeholder="X channel"),
                 show_label=False, width=90,
                 tooltip="X channel (or single channel for 1D plots)"),
            Item('channel_y',
                 editor=TextEditor(auto_set=False, placeholder="Y channel"),
                 show_label=False, width=90,
                 tooltip="Y channel (2D plots only)"),
            Item('scale_x',
                 editor=EnumEditor(values=['logicle', 'log', 'linear']),
                 show_label=False, width=70,
                 tooltip="X axis scale"),
            Item('title',
                 editor=TextEditor(auto_set=False, placeholder="Panel title"),
                 show_label=False, width=100,
                 tooltip="Optional panel title"),
        )
    )


class FigureLayoutHandler(ViewHandler):

    add_panel = Button(label="+ Panel")

    view_traits_view = \
        View(VGroup(
             VGroup(
                 HGroup(
                     Item('rows',
                          editor=TextEditor(auto_set=False, evaluate=int,
                                            format_func=str),
                          label="Rows", width=50),
                     Item('cols',
                          editor=TextEditor(auto_set=False, evaluate=int,
                                            format_func=str),
                          label="Cols", width=50),
                 ),
                 Item('title',
                      editor=TextEditor(auto_set=False, placeholder="Figure title"),
                      label="Title"),
                 HGroup(
                     Item('figsize_w',
                          editor=TextEditor(auto_set=False, evaluate=int,
                                            format_func=str),
                          label="Panel width (in)", width=50),
                     Item('figsize_h',
                          editor=TextEditor(auto_set=False, evaluate=int,
                                            format_func=str),
                          label="height (in)", width=50),
                 ),
                 label="Layout",
                 show_border=False),
             VGroup(
                 HGroup(
                     Label("Type"),
                     Spring(),
                     Label("X channel"),
                     Label("Y channel"),
                     Label("X scale"),
                     Label("Panel title"),
                     visible_when='object.panels'),
                 Item('panels',
                      editor=VerticalListEditor(
                          editor=InstanceHandlerEditor(
                              handler_factory=FigurePanelHandler,
                              view='default_view'),
                          style='custom',
                          mutable=False,
                          deletable=True),
                      show_label=False),
                 Item('handler.add_panel',
                      editor=ButtonEditor(value=True),
                      show_label=False),
                 label="Panels",
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
                                         background_color="#ff9191"))))

    view_params_view = View()

    @observe('add_panel')
    def _on_add_panel(self, _):
        self.model.panels.append(FigurePanel())


@provides(IViewPlugin)
class FigureLayoutPlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.view.figure_layout'
    view_id = 'cytoflow.view.figure_layout'

    name = "Figure Layout"
    short_name = "Layout"

    def get_view(self):
        return FigureLayoutWorkflowView()

    def get_handler(self, model, context):
        return FigureLayoutHandler(model=model, context=context)

    def get_icon(self):
        return ImageResource('scatterplot')

    plugin = List(contributes_to=VIEW_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
