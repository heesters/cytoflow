#!/usr/bin/env python3.8

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

'''
Time Gate
---------

Gate on the instrument **Time** channel to exclude unstable acquisition
periods — debris at the start of a run, clogs, or pressure fluctuations.

Draw the time window by clicking and dragging on the histogram.  Drag
either edge to adjust, or drag the middle to move the whole window.

The time channel is auto-detected from the experiment (looks for channels
named ``Time``, ``TIME``, ``time``, etc.).  You can override it with the
**Channel** dropdown.

.. object:: Name

    The operation name.  Used to name the new metadata field created by
    this module.

.. object:: Channel

    The time channel.  Auto-detected if left as the first entry.

.. object:: Low / High

    The time window to keep (in the channel's native units, usually ms).

.. object:: Scale

    The axis scale.  Normally **linear** for time.

.. object:: Hue facet

    Show different experimental conditions in different colours.

.. object:: Subset

    Show only a subset of the data.
'''

from traits.api import provides, List
from traitsui.api import View, Item, EnumEditor, VGroup, TextEditor
from envisage.api import Plugin
from pyface.api import ImageResource  # @UnresolvedImport

from ..view_plugins import ViewHandler
from ..view_plugins.histogram import HistogramParamsHandler
from ..editors import SubsetListEditor, ColorTextEditor, ExtendableEnumEditor, InstanceHandlerEditor
from ..workflow.operations import TimeGateWorkflowOp
from ..workflow.operations.time_gate import TimeGateSelectionView
from ..subset_controllers import subset_handler_factory

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, shared_op_traits_view, PluginHelpMixin


class TimeGateHandler(OpHandler):

    operation_traits_view = \
        View(Item('name',
                  editor=TextEditor(auto_set=False)),
             Item('channel',
                  editor=ExtendableEnumEditor(
                      name='context_handler.previous_channels',
                      extra_items={'Auto-detect': ''}),
                  label='Channel',
                  tooltip='Time channel (auto-detected if blank)'),
             Item('low',
                  editor=TextEditor(auto_set=False,
                                    evaluate=float,
                                    format_func=lambda x: '' if x is None else str(x),
                                    placeholder='None'),
                  label='Low'),
             Item('high',
                  editor=TextEditor(auto_set=False,
                                    evaluate=float,
                                    format_func=lambda x: '' if x is None else str(x),
                                    placeholder='None'),
                  label='High'),
             shared_op_traits_view)


class TimeGateViewHandler(ViewHandler):
    view_traits_view = \
        View(VGroup(
            VGroup(Item('channel',
                        label='Channel',
                        style='readonly'),
                   Item('scale'),
                   Item('huefacet',
                        editor=ExtendableEnumEditor(
                            name='context_handler.previous_conditions_names',
                            extra_items={'None': ''}),
                        label='Color\nFacet'),
                   label='Plot settings',
                   show_border=False),
            VGroup(Item('subset_list',
                        show_label=False,
                        editor=SubsetListEditor(
                            conditions='context_handler.previous_conditions',
                            editor=InstanceHandlerEditor(
                                view='subset_view',
                                handler_factory=subset_handler_factory),
                            mutable=False)),
                   label='Subset',
                   show_border=False,
                   show_labels=False),
            Item('context.view_warning',
                 resizable=True,
                 visible_when='context.view_warning',
                 editor=ColorTextEditor(foreground_color='#000000',
                                        background_color='#ffff99')),
            Item('context.view_error',
                 resizable=True,
                 visible_when='context.view_error',
                 editor=ColorTextEditor(foreground_color='#000000',
                                        background_color='#ff9191'))))

    view_params_view = \
        View(Item('plot_params',
                  editor=InstanceHandlerEditor(view='view_params_view',
                                               handler_factory=HistogramParamsHandler),
                  style='custom',
                  show_label=False))


@provides(IOperationPlugin)
class TimeGatePlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.time_gate'
    operation_id = 'cytoflow.operations.time_gate'
    view_id = 'cytoflow.views.range'   # reuses range view machinery

    name = 'Time Gate'
    short_name = 'Time'
    menu_group = 'Gates'

    def get_operation(self):
        return TimeGateWorkflowOp()

    def get_handler(self, model, context):
        if isinstance(model, TimeGateWorkflowOp):
            return TimeGateHandler(model=model, context=context)
        elif isinstance(model, TimeGateSelectionView):
            return TimeGateViewHandler(model=model, context=context)

    def get_icon(self):
        return ImageResource('range')

    plugin = List(contributes_to=OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
