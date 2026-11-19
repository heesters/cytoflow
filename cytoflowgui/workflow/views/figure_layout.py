#!/usr/bin/env python3.8

# (c) Heesters Lab, Utrecht University 2024
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.

"""
cytoflowgui.workflow.views.figure_layout
-----------------------------------------
"""

from traits.api import provides, Str, Int, List, Instance, Bool

from cytoflow.views.figure_layout import FigureLayoutView, FigurePanel

from cytoflowgui.workflow.serialization import camel_registry, cytoflow_class_repr
from .view_base import IWorkflowView, WorkflowView

FigureLayoutView.__repr__ = cytoflow_class_repr
FigurePanel.__repr__ = cytoflow_class_repr


@provides(IWorkflowView)
class FigureLayoutWorkflowView(WorkflowView, FigureLayoutView):
    rows      = Int(2, apply=True)
    cols      = Int(2, apply=True)
    title     = Str('', apply=True)
    panels    = List(Instance(FigurePanel), apply=True)
    figsize_w = Int(5, apply=True)
    figsize_h = Int(4, apply=True)

    def get_notebook_code(self, idx):
        panels_repr = "[" + ", ".join(repr(p) for p in self.panels) + "]"
        return (f"\nfrom cytoflow.views.figure_layout import FigureLayoutView, FigurePanel\n"
                f"view_{idx} = FigureLayoutView(rows={self.rows}, cols={self.cols},\n"
                f"    title={self.title!r}, panels={panels_repr})\n"
                f"view_{idx}.plot(ex_{idx - 1})")


# --- serialization ---

@camel_registry.dumper(FigurePanel, 'FigurePanel', version=1)
def _dump_panel(p):
    return dict(view_type=p.view_type, channel_x=p.channel_x,
                channel_y=p.channel_y, scale_x=p.scale_x, scale_y=p.scale_y,
                huefacet=p.huefacet, subset=p.subset, title=p.title)


@camel_registry.loader('FigurePanel', version=1)
def _load_panel(data, version):
    return FigurePanel(**{k: v for k, v in data.items() if v})


@camel_registry.dumper(FigureLayoutWorkflowView, 'FigureLayoutWorkflowView', version=1)
def _dump(v):
    return dict(rows=v.rows, cols=v.cols, title=v.title,
                panels=v.panels, figsize_w=v.figsize_w, figsize_h=v.figsize_h)


@camel_registry.loader('FigureLayoutWorkflowView', version=1)
def _load(data, version):
    return FigureLayoutWorkflowView(
        rows=data.get('rows', 2),
        cols=data.get('cols', 2),
        title=data.get('title', ''),
        panels=data.get('panels', []),
        figsize_w=data.get('figsize_w', 5),
        figsize_h=data.get('figsize_h', 4))
