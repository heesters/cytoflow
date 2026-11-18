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
cytoflowgui.workflow.views.calcium_flux
---------------------------------------

"""

from traits.api import provides, Instance

from cytoflow.operations.calcium_flux import CalciumFluxView

from .view_base import IWorkflowView, WorkflowView
from ..serialization import camel_registry, traits_str, dedent


@provides(IWorkflowView)
class CalciumFluxWorkflowView(WorkflowView, CalciumFluxView):
    plot_params = Instance(object, ())

    def should_plot(self, changed, payload):
        from .. import Changed
        if changed in (Changed.APPLY, Changed.PREV_RESULT):
            return True
        return False

    def get_notebook_code(self, idx):
        view = CalciumFluxView()
        view.copy_traits(self, view.copyable_trait_names())

        return dedent("""
        op_{idx}.default_view({traits}).plot(ex_{idx})
        """
        .format(traits=traits_str(view),
                idx=idx))


### Serialization
@camel_registry.dumper(CalciumFluxWorkflowView, 'calcium-flux-view', version=1)
def _dump_view(view):
    return dict(op=view.op,
                huefacet=view.huefacet,
                rowfacet=view.rowfacet,
                colfacet=view.colfacet,
                show_error_band=view.show_error_band,
                show_stimulation_line=view.show_stimulation_line)


@camel_registry.loader('calcium-flux-view', version=1)
def _load_view(data, version):
    return CalciumFluxWorkflowView(**data)
