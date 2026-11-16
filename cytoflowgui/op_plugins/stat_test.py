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
Statistical Test
----------------

Test whether a quantity differs between the groups defined by an
experimental condition -- one-way ANOVA, Kruskal-Wallis, Welch's t-test or
Mann-Whitney U -- optionally followed by a Tukey HSD comparison of every pair
of groups.  The plot shows the mean +/- SEM of each group, the individual
replicates, and brackets with stars over the pairs of groups that differ
significantly (``*`` p < 0.05, ``**`` p < 0.01, ``***`` p < 0.001).

**Which data are compared?**

- **Statistic (replicates)** -- the usual choice.  Compare the values of a
  statistic you computed earlier, eg with a **Channel Statistic**.  The
  groups are the values of **Groups (By)**; the observations in each group
  are the statistic's values for all its *other* conditions -- the
  replicates.  For example, to compare the geometric mean of FITC-A between
  Dox concentrations, where each concentration was measured in several
  wells, compute a **Channel Statistic** (Geo.Mean of FITC-A) grouped by
  *both* Dox and Well, then test it with **Groups (By)** = Dox.  Each
  group needs at least two replicates: a statistic grouped by Dox alone has
  only one value per concentration, and can't be tested.  Use **Subset** to
  test within one value of another condition.

- **Individual events** -- compare the raw values of a **Channel** in the
  individual events (cells) of each group.  **Be careful**: with tens of
  thousands of events per group, almost every difference, however small, is
  "statistically significant" -- and cells from the same tube are not
  independent replicates.  Pay attention to the effect size, not the
  p-value.

.. object:: Name

    The operation name.  The results are added as the statistics
    ``<Name>`` (the overall test: ``statistic``, ``p_value``, degrees of
    freedom, ``effect_size``), ``<Name>_Groups`` (the mean, SD, SEM etc.
    of each group) and ``<Name>_Tukey`` (the pairwise comparisons, if you
    chose the Tukey HSD post-hoc test).  View them with the **Table View**:
    for example, ``<Name>_Tukey`` with **Rows** = group1, **Columns** =
    group2 and **Feature** = p_adj; or ``<Name>`` with **Rows** = Test and
    **Feature** = p_value.

.. object:: Compare

    Compare the replicate values of a statistic, or individual events?

.. object:: Statistic

    The statistic to compare (when comparing a statistic).

.. object:: Feature

    The feature (column) of the statistic to compare.

.. object:: Channel

    The channel to compare (when comparing individual events).

.. object:: Groups (By)

    The condition whose values are the groups to compare.

.. object:: Test

    - **One-way ANOVA** -- do the group means differ?  Assumes roughly
      normal data with similar variances.  Effect size: eta-squared.
    - **Kruskal-Wallis** -- the rank-based (non-parametric) version of the
      ANOVA.  Effect size: epsilon-squared.
    - **Welch t-test** -- compare the means of exactly two groups, without
      assuming equal variances.  Effect size: Cohen's d.
    - **Mann-Whitney U** -- the rank-based version of the t-test, exactly
      two groups.  Effect size: rank-biserial correlation.

.. object:: Post-hoc test

    **Tukey HSD** compares every pair of groups.  Its p-values (``p_adj``)
    and confidence intervals are **adjusted for multiple comparisons**:
    they keep the chance of *any* false positive among all the pairs at
    **alpha**, so don't correct them again.  Tukey HSD is the follow-up to
    a one-way ANOVA.

.. object:: alpha

    The significance level, used for the Tukey HSD ``reject`` column and
    confidence intervals, and for the significance stars.

.. object:: Subset

    Only test a subset of the data.

.. object:: Show non-significant brackets

    Also draw brackets (marked "ns") for the pairs that don't differ
    significantly.

.. plot::
   :include-source: False

    import cytoflow as flow
    import_op = flow.ImportOp()
    import_op.tubes = [flow.Tube(file = "Plate01/CFP_Well_A4.fcs",
                                 conditions = {'Dox' : 1.0, 'Well' : 'A'}),
                       flow.Tube(file = "Plate01/RFP_Well_A3.fcs",
                                 conditions = {'Dox' : 1.0, 'Well' : 'B'}),
                       flow.Tube(file = "Plate01/YFP_Well_A7.fcs",
                                 conditions = {'Dox' : 10.0, 'Well' : 'A'}),
                       flow.Tube(file = "Plate01/CFP_Well_B4.fcs",
                                 conditions = {'Dox' : 10.0, 'Well' : 'B'}),
                       flow.Tube(file = "Plate01/RFP_Well_A6.fcs",
                                 conditions = {'Dox' : 100.0, 'Well' : 'A'}),
                       flow.Tube(file = "Plate01/YFP_Well_C7.fcs",
                                 conditions = {'Dox' : 100.0, 'Well' : 'B'})]
    import_op.conditions = {'Dox' : 'float', 'Well' : 'category'}
    ex = import_op.apply()

    ex2 = flow.ChannelStatisticOp(name = 'GeoMean',
                                  channel = 'Y2-A',
                                  function = flow.geom_mean,
                                  by = ['Dox', 'Well']).apply(ex)

    test_op = flow.StatisticalTestOp(name = 'DoxTest',
                                     statistic = 'GeoMean',
                                     feature = 'Y2-A',
                                     by = 'Dox',
                                     test = 'anova',
                                     posthoc = 'tukey')
    ex3 = test_op.apply(ex2)
    test_op.default_view(show_ns = True).plot(ex3)

'''

import numpy as np
import pandas as pd

from traits.api import provides, Property, List
from traitsui.api import (View, Item, TextEditor, VGroup, EnumEditor)
from pyface.api import ImageResource  # @UnresolvedImport
from envisage.api import Plugin

from ..workflow.operations import StatisticalTestWorkflowOp, StatisticalTestWorkflowView
from ..view_plugins import ViewHandler
from ..editors import SubsetListEditor, InstanceHandlerEditor, ColorTextEditor
from ..subset_controllers import subset_handler_factory

from .i_op_plugin import IOperationPlugin, OP_PLUGIN_EXT
from .op_plugin_base import OpHandler, PluginHelpMixin, shared_op_traits_view


class StatisticalTestHandler(OpHandler):

    # the choices for "by"
    indices = Property(depends_on = "context.previous_wi.statistics, "
                                    "context.previous_wi.conditions, "
                                    "model.statistic, model.observations")

    # the conditions (and their values) for the subset editor
    levels = Property(depends_on = "context.previous_wi.statistics, "
                                   "context.previous_wi.conditions, "
                                   "model.statistic, model.observations")

    features = Property(depends_on = "context.previous_wi.statistics, model.statistic")

    def _get_stat(self):
        if not (self.context
                and self.context.previous_wi
                and self.context.previous_wi.statistics
                and self.model
                and self.model.statistic
                and self.model.statistic in self.context.previous_wi.statistics):
            return None

        return self.context.previous_wi.statistics[self.model.statistic]

    def _get_conditions(self):
        if (self.context 
            and self.context.previous_wi 
            and self.context.previous_wi.conditions):
            return self.context.previous_wi.conditions
        
        return {}
    
    # MAGIC: gets the value for the property 'indices'
    def _get_indices(self):
        if self.model and self.model.observations == "events":
            return sorted(self._get_conditions().keys())

        stat = self._get_stat()
        if stat is None:
            return []

        return list(stat.index.names)

    # MAGIC: gets the value for the property 'levels'
    # returns a Dict(Str, pd.Series)
    def _get_levels(self):
        if self.model and self.model.observations == "events":
            return self._get_conditions()

        stat = self._get_stat()
        if stat is None:
            return {}

        return {name : pd.Series(np.unique(stat.index.get_level_values(name)))
                for name in stat.index.names}

    # MAGIC: gets the value for the property 'features'
    def _get_features(self):
        stat = self._get_stat()
        if stat is None:
            return []

        if isinstance(stat, pd.Series):
            return [stat.name] if stat.name else []

        return stat.columns.to_list()

    operation_traits_view = \
        View(Item('name',
                  editor = TextEditor(auto_set = False,
                                      placeholder = "None")),
             Item('observations',
                  editor = EnumEditor(values = {"statistic" : "1:Statistic (replicates)",
                                                "events" : "2:Individual events"}),
                  label = "Compare"),
             Item('statistic',
                  editor = EnumEditor(name = 'context_handler.previous_statistics_names'),
                  visible_when = 'observations == "statistic"',
                  label = "Statistic"),
             Item('feature',
                  editor = EnumEditor(name = 'handler.features'),
                  visible_when = 'observations == "statistic"',
                  label = "Feature"),
             Item('channel',
                  editor = EnumEditor(name = 'context_handler.previous_channels'),
                  visible_when = 'observations == "events"',
                  label = "Channel"),
             Item('by',
                  editor = EnumEditor(name = 'handler.indices'),
                  label = "Groups\n(By)"),
             Item('test',
                  editor = EnumEditor(values = {"anova" : "1:One-way ANOVA",
                                                "kruskal" : "2:Kruskal-Wallis",
                                                "welch" : "3:Welch t-test (2 groups)",
                                                "mannwhitney" : "4:Mann-Whitney U (2 groups)"}),
                  label = "Test"),
             Item('posthoc',
                  editor = EnumEditor(values = {"none" : "1:None",
                                                "tukey" : "2:Tukey HSD"}),
                  label = "Post-hoc\ntest"),
             Item('alpha',
                  editor = TextEditor(auto_set = False,
                                      evaluate = float,
                                      format_func = lambda x: "" if x is None else str(x)),
                  label = "alpha"),
             VGroup(Item('subset_list',
                         show_label = False,
                         editor = SubsetListEditor(conditions = "handler.levels",
                                                   editor = InstanceHandlerEditor(view = 'subset_view',
                                                                                  handler_factory = subset_handler_factory))),
                    label = "Subset",
                    show_border = False,
                    show_labels = False),
             shared_op_traits_view)


class StatisticalTestViewHandler(ViewHandler):
    view_traits_view = \
        View(Item('show_ns',
                  label = "Show non-\nsignificant\nbrackets"),
             Item('context.view_warning',
                  resizable = True,
                  visible_when = 'context.view_warning',
                  editor = ColorTextEditor(foreground_color = "#000000",
                                           background_color = "#ffff99")),
             Item('context.view_error',
                  resizable = True,
                  visible_when = 'context.view_error',
                  editor = ColorTextEditor(foreground_color = "#000000",
                                           background_color = "#ff9191")))

    view_params_view = View()


@provides(IOperationPlugin)
class StatisticalTestPlugin(Plugin, PluginHelpMixin):

    id = 'cytoflowgui.op_plugins.statistical_test'
    operation_id = 'cytoflow.operations.statistical_test'
    view_id = 'cytoflow.view.statisticaltestview'

    name = "Statistical Test"
    short_name = "Stat\nTest"
    menu_group = "Statistics"

    def get_operation(self):
        return StatisticalTestWorkflowOp()

    def get_handler(self, model, context):
        if isinstance(model, StatisticalTestWorkflowOp):
            return StatisticalTestHandler(model = model, context = context)
        elif isinstance(model, StatisticalTestWorkflowView):
            return StatisticalTestViewHandler(model = model, context = context)

    def get_icon(self):
        return ImageResource('stat_test')

    plugin = List(contributes_to = OP_PLUGIN_EXT)
    def _plugin_default(self):
        return [self]
