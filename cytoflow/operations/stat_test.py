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
cytoflow.operations.stat_test
-----------------------------

Statistical hypothesis tests (ANOVA, Kruskal-Wallis, Welch's t-test,
Mann-Whitney U, Tukey HSD) comparing groups defined by an experimental
condition.  `stat_test` has two classes:

`StatisticalTestOp` -- runs the test and adds the results to
`Experiment.statistics`.

`StatisticalTestView` -- plots the group means +/- SEM with significance
brackets.
"""

import warnings
from collections import OrderedDict

import numpy as np
import pandas as pd
import scipy.stats

from traits.api import (HasStrictTraits, Str, Constant, Enum, Float, Instance,
                        Bool, provides)

import cytoflow.utility as util

from .i_operation import IOperation

#: The tests `StatisticalTestOp` can run, and their human-readable names
TESTS = OrderedDict([("anova", "One-way ANOVA"),
                     ("kruskal", "Kruskal-Wallis"),
                     ("welch", "Welch's t-test"),
                     ("mannwhitney", "Mann-Whitney U")])

#: The post-hoc tests `StatisticalTestOp` can run
POSTHOC_TESTS = OrderedDict([("none", "None"),
                             ("tukey", "Tukey HSD")])

# the tests that compare exactly two groups
_TWO_GROUP_TESTS = ("welch", "mannwhitney")


def significance_stars(p, alpha = 0.05):
    """
    Convert a p-value to the usual "star" notation: ``***`` for p < 0.001,
    ``**`` for p < 0.01, ``*`` for p < `alpha` (0.05 by default) and ``ns``
    otherwise.
    """
    if p is None or np.isnan(p):
        return ""
    if p < 0.001 and p < alpha:
        return "***"
    if p < 0.01 and p < alpha:
        return "**"
    if p < alpha:
        return "*"
    return "ns"


@provides(IOperation)
class StatisticalTestOp(HasStrictTraits):
    """
    Test whether a quantity differs between the groups defined by an
    experimental condition: one-way ANOVA, Kruskal-Wallis, Welch's t-test or
    Mann-Whitney U, optionally followed by Tukey's HSD post-hoc test of every
    pair of groups.

    **What is compared.**  There are two modes, set by `observations`:

    - ``"statistic"`` (the default, and usually what you want): compare the
      values of a feature of an existing statistic.  The groups are the
      levels of the condition `by`; the observations in each group are the
      statistic's values for every *other* index level -- ie, the replicates.
      For example, if the statistic ``GeoMean`` is the geometric mean of
      ``FITC-A`` computed ``by = ["Dox", "Well"]``, and each ``Dox``
      concentration was measured in several wells, then a test with
      ``by = "Dox"`` compares the per-well geometric means between ``Dox``
      concentrations, with the wells as replicates.  Each group needs at
      least two replicates: a statistic computed only ``by = ["Dox"]`` has
      one value per group, and can't be tested.  All other index levels of
      the statistic are pooled as replicates; use `subset` to test within
      one level of another condition.

    - ``"events"``: compare the raw values of `channel` for the individual
      events (cells) in each group.

      .. warning::
          With 10^4 - 10^6 events per group, almost every difference,
          however tiny, is "statistically significant".  Look at the
          effect size (and the size of the difference itself) rather than
          the p-value, and remember that cells from the same tube are not
          independent biological replicates.  Testing a per-tube statistic
          in ``"statistic"`` mode is almost always the better choice.

    **Tests.**

    - ``"anova"``: one-way ANOVA (`scipy.stats.f_oneway`).  Effect size:
      eta-squared (between-group sum of squares / total sum of squares).
    - ``"kruskal"``: Kruskal-Wallis H test (`scipy.stats.kruskal`), the
      rank-based (non-parametric) alternative to ANOVA.  Effect size:
      epsilon-squared, H / (N - 1).
    - ``"welch"``: Welch's unequal-variances t-test
      (`scipy.stats.ttest_ind` with ``equal_var = False``).  Exactly two
      groups.  Effect size: Cohen's d, (mean1 - mean2) / sqrt((sd1^2 + sd2^2) / 2).
    - ``"mannwhitney"``: Mann-Whitney U test (`scipy.stats.mannwhitneyu`,
      two-sided).  Exactly two groups.  Effect size: rank-biserial
      correlation, 2 U / (n1 n2) - 1.

    If `posthoc` is ``"tukey"``, every pair of groups is also compared with
    Tukey's honestly-significant-difference test
    (`statsmodels.stats.multicomp.pairwise_tukeyhsd`).  Its p-values and
    confidence intervals are adjusted for multiple comparisons (they control
    the family-wise error rate at `alpha` across all the pairs), so you don't
    need to correct them further.  Tukey's HSD is the post-hoc test that
    goes with the one-way ANOVA; it assumes roughly normal data with similar
    variances in each group.

    **Results.**  `apply` adds two or three statistics:

    - `name` -- the overall test.  One row, indexed by ``Test`` (the test's
      name), with the columns ``statistic`` (F, H, t or U), ``p_value``,
      ``df_between`` and ``df_within`` (the degrees of freedom: ``df_between``
      is the numerator / only degrees of freedom, ``NaN`` if the test has
      none), ``effect_size``, ``n_groups`` and ``n_obs``.
    - `name` + ``_Groups`` -- descriptive statistics for each group, indexed
      by `by`: ``n``, ``mean``, ``sd``, ``sem`` and ``median``.
    - `name` + ``_Tukey`` (only if `posthoc` is ``"tukey"``) -- one row per
      pair of groups, indexed by ``group1`` and ``group2``, with the columns
      ``meandiff`` (mean of ``group2`` minus mean of ``group1``), ``p_adj``
      (the adjusted p-value), ``lower`` and ``upper`` (the adjusted confidence
      interval of ``meandiff``), ``reject`` (``True`` if ``p_adj`` < `alpha`)
      and ``significance`` (``***``, ``**``, ``*`` or ``ns``).

    Attributes
    ----------
    name : Str
        The operation name.  Becomes the name of the new statistic(s).

    observations : Enum("statistic", "events") (default = "statistic")
        What to compare: the values of a statistic (replicates), or the
        individual events.  See above.

    statistic : Str
        In ``"statistic"`` mode, the statistic to test.

    feature : Str
        In ``"statistic"`` mode, the feature (column) of `statistic` to test.

    channel : Str
        In ``"events"`` mode, the channel whose values are compared.

    by : Str
        The condition whose levels define the groups to compare.  In
        ``"statistic"`` mode it must be one of the statistic's index levels.

    test : Enum("anova", "kruskal", "welch", "mannwhitney") (default = "anova")
        Which test to run.

    posthoc : Enum("none", "tukey") (default = "none")
        Which post-hoc pairwise test to run.

    alpha : Float (default = 0.05)
        The significance level used for the post-hoc test's ``reject``
        column, its confidence intervals and the significance stars.

    subset : Str
        A Python expression used to select a subset of the data to test.
        In ``"statistic"`` mode, it is applied to the statistic (with
        `pandas.DataFrame.query`) and may refer to its index levels; in
        ``"events"`` mode, it is passed to `Experiment.query`.

    Examples
    --------

    .. plot::
        :context: close-figs

        Make a little data set, with two replicate wells for each
        ``Dox`` concentration.

        >>> import cytoflow as flow
        >>> import_op = flow.ImportOp()
        >>> import_op.tubes = [flow.Tube(file = "Plate01/CFP_Well_A4.fcs",
        ...                              conditions = {'Dox' : 1.0, 'Well' : 'A'}),
        ...                    flow.Tube(file = "Plate01/RFP_Well_A3.fcs",
        ...                              conditions = {'Dox' : 1.0, 'Well' : 'B'}),
        ...                    flow.Tube(file = "Plate01/YFP_Well_A7.fcs",
        ...                              conditions = {'Dox' : 10.0, 'Well' : 'A'}),
        ...                    flow.Tube(file = "Plate01/CFP_Well_B4.fcs",
        ...                              conditions = {'Dox' : 10.0, 'Well' : 'B'}),
        ...                    flow.Tube(file = "Plate01/RFP_Well_A6.fcs",
        ...                              conditions = {'Dox' : 100.0, 'Well' : 'A'}),
        ...                    flow.Tube(file = "Plate01/YFP_Well_C7.fcs",
        ...                              conditions = {'Dox' : 100.0, 'Well' : 'B'})]
        >>> import_op.conditions = {'Dox' : 'float', 'Well' : 'category'}
        >>> ex = import_op.apply()

    Compute the geometric mean of each well (tube).

    .. plot::
        :context: close-figs

        >>> ex2 = flow.ChannelStatisticOp(name = 'GeoMean',
        ...                               channel = 'Y2-A',
        ...                               function = flow.geom_mean,
        ...                               by = ['Dox', 'Well']).apply(ex)

    Test whether the per-well geometric means differ between ``Dox``
    concentrations, then compare each pair of concentrations.

    .. plot::
        :context: close-figs

        >>> test_op = flow.StatisticalTestOp(name = 'DoxTest',
        ...                                  statistic = 'GeoMean',
        ...                                  feature = 'Y2-A',
        ...                                  by = 'Dox',
        ...                                  test = 'anova',
        ...                                  posthoc = 'tukey')
        >>> ex3 = test_op.apply(ex2)
        >>> test_op.default_view().plot(ex3)

    The results are new statistics.

    .. plot::
        :context: close-figs

        >>> ex3.statistics['DoxTest']
        >>> ex3.statistics['DoxTest_Tukey']
    """

    id = Constant('cytoflow.operations.statistical_test')
    friendly_id = Constant("Statistical Test")

    name = Str
    observations = Enum("statistic", "events")
    statistic = Str
    feature = Str
    channel = Str
    by = Str
    test = Enum(*TESTS.keys())
    posthoc = Enum(*POSTHOC_TESTS.keys())
    alpha = Float(0.05)
    subset = Str

    def apply(self, experiment):
        """
        Run the test.

        Parameters
        ----------
        experiment : `Experiment`
            The `Experiment` containing the statistic (or events) to test.

        Returns
        -------
        Experiment
            A new `Experiment` (a shallow copy of the old one) with the
            results added to `Experiment.statistics` (see the class
            documentation).
        """

        if experiment is None:
            raise util.CytoflowOpError('experiment',
                                       "Must specify an experiment")

        if not self.name:
            raise util.CytoflowOpError('name',
                                       "Must specify a name")

        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError('name',
                                       "Name can only contain letters, numbers and underscores.")

        new_names = [self.name, self.name + "_Groups"]
        if self.posthoc == "tukey":
            new_names.append(self.name + "_Tukey")

        for n in new_names:
            if n in experiment.statistics:
                raise util.CytoflowOpError('name',
                                           "{} is already in the experiment's statistics"
                                           .format(n))

        if not (0.0 < self.alpha < 1.0):
            raise util.CytoflowOpError('alpha',
                                       "alpha must be between 0 and 1")

        groups = self._get_groups(experiment)
        k = len(groups)
        values = list(groups.values())
        n_obs = int(sum(len(v) for v in values))

        # the overall test
        if self.test == "anova":
            res = scipy.stats.f_oneway(*values)
            stat, p = res.statistic, res.pvalue
            df_between, df_within = k - 1, n_obs - k
            grand_mean = np.concatenate(values).mean()
            ss_between = sum(len(v) * (v.mean() - grand_mean) ** 2 for v in values)
            ss_total = ((np.concatenate(values) - grand_mean) ** 2).sum()
            effect_size = ss_between / ss_total if ss_total > 0 else np.nan
        elif self.test == "kruskal":
            try:
                res = scipy.stats.kruskal(*values)
            except ValueError as e:
                raise util.CytoflowOpError('test',
                                           "Kruskal-Wallis test failed: {}"
                                           .format(e)) from e
            stat, p = res.statistic, res.pvalue
            df_between, df_within = k - 1, np.nan
            effect_size = stat / (n_obs - 1)
        elif self.test == "welch":
            a, b = values
            res = scipy.stats.ttest_ind(a, b, equal_var = False)
            stat, p = res.statistic, res.pvalue
            df_between, df_within = res.df, np.nan
            sd_pooled = np.sqrt((a.var(ddof = 1) + b.var(ddof = 1)) / 2.0)
            effect_size = (a.mean() - b.mean()) / sd_pooled if sd_pooled > 0 else np.nan
        elif self.test == "mannwhitney":
            a, b = values
            res = scipy.stats.mannwhitneyu(a, b, alternative = 'two-sided')
            stat, p = res.statistic, res.pvalue
            df_between, df_within = np.nan, np.nan
            effect_size = 2.0 * stat / (len(a) * len(b)) - 1.0

        if np.isnan(p):
            warnings.warn("The test's p-value is NaN -- are all the "
                          "observations identical?",
                          util.CytoflowOpWarning)

        overall = pd.DataFrame({"statistic" : [float(stat)],
                                "p_value" : [float(p)],
                                "df_between" : [float(df_between)],
                                "df_within" : [float(df_within)],
                                "effect_size" : [float(effect_size)],
                                "n_groups" : [k],
                                "n_obs" : [n_obs]},
                               index = pd.Index([TESTS[self.test]], name = "Test"))

        # descriptive statistics for each group
        desc = pd.DataFrame({"n" : [len(v) for v in values],
                             "mean" : [v.mean() for v in values],
                             "sd" : [v.std(ddof = 1) for v in values],
                             "sem" : [scipy.stats.sem(v) for v in values],
                             "median" : [np.median(v) for v in values]},
                            index = pd.Index(list(groups.keys()), name = self.by))

        new_experiment = experiment.clone(deep = False)
        new_experiment.statistics[self.name] = overall
        new_experiment.statistics[self.name + "_Groups"] = desc

        # post-hoc test
        if self.posthoc == "tukey":
            from statsmodels.stats.multicomp import pairwise_tukeyhsd

            endog = np.concatenate(values)
            labels = np.concatenate([np.full(len(v), i) for i, v in enumerate(values)])

            res = pairwise_tukeyhsd(endog, labels, alpha = self.alpha)

            keys = list(groups.keys())
            i1, i2 = res._multicomp.pairindices
            g1 = [keys[int(res.groupsunique[i])] for i in i1]
            g2 = [keys[int(res.groupsunique[i])] for i in i2]

            pvalues = np.asarray(res.pvalues, dtype = float)
            confint = np.asarray(res.confint, dtype = float)

            tukey = pd.DataFrame({"meandiff" : np.asarray(res.meandiffs, dtype = float),
                                  "p_adj" : pvalues,
                                  "lower" : confint[:, 0],
                                  "upper" : confint[:, 1],
                                  "reject" : np.asarray(res.reject, dtype = bool),
                                  "significance" : [significance_stars(x, self.alpha) for x in pvalues]},
                                 index = pd.MultiIndex.from_arrays([g1, g2], names = ["group1", "group2"]))
            new_experiment.statistics[self.name + "_Tukey"] = tukey

        new_experiment.history.append(self.clone_traits(transient = lambda _: True))
        return new_experiment

    def _get_groups(self, experiment):
        """
        Returns an `OrderedDict` mapping each group (a level of `by`, in sorted
        order) to a 1D `numpy.ndarray` of its observations.  Raises
        `CytoflowOpError` if the parameters are invalid or there are too few
        groups or observations.
        """

        if not self.by:
            raise util.CytoflowOpError('by',
                                       "Must specify the condition whose levels "
                                       "define the groups ('by')")

        if self.observations == "statistic":
            if not self.statistic:
                raise util.CytoflowOpError('statistic',
                                           "Statistic not set")

            if self.statistic not in experiment.statistics:
                raise util.CytoflowOpError('statistic',
                                           "Can't find the statistic {} in the experiment"
                                           .format(self.statistic))

            stat = experiment.statistics[self.statistic]
            if isinstance(stat, pd.Series):
                stat = stat.to_frame()

            if not self.feature:
                raise util.CytoflowOpError('feature',
                                           "Must set a feature")

            if self.feature not in stat:
                raise util.CytoflowOpError('feature',
                                           "Can't find feature {} in statistic {}. "
                                           "Must be one of {}"
                                           .format(self.feature, self.statistic,
                                                   stat.columns.to_list()))

            if self.by not in stat.index.names:
                raise util.CytoflowOpError('by',
                                           "{} is not an index level of statistic {}; "
                                           "must be one of {}"
                                           .format(self.by, self.statistic,
                                                   list(stat.index.names)))

            if self.subset:
                try:
                    stat = stat.query(self.subset)
                except Exception as e:
                    raise util.CytoflowOpError('subset',
                                               "Subset string '{}' isn't valid"
                                               .format(self.subset)) from e
                if len(stat) == 0:
                    raise util.CytoflowOpError('subset',
                                               "Subset string '{}' returned no values"
                                               .format(self.subset))

            data = pd.to_numeric(stat[self.feature], errors = 'coerce')
            group_values = stat.index.get_level_values(self.by)
            unit = "replicate"

        else:
            if not self.channel:
                raise util.CytoflowOpError('channel',
                                           "Must set a channel")

            if self.channel not in experiment.channels:
                raise util.CytoflowOpError('channel',
                                           "Channel {} not in the experiment"
                                           .format(self.channel))

            if self.by not in experiment.conditions:
                raise util.CytoflowOpError('by',
                                           "{} is not a condition in the experiment"
                                           .format(self.by))

            if self.subset:
                try:
                    experiment = experiment.query(self.subset)
                except Exception as e:
                    raise util.CytoflowOpError('subset',
                                               "Subset string '{}' isn't valid"
                                               .format(self.subset)) from e

                if len(experiment) == 0:
                    raise util.CytoflowOpError('subset',
                                               "Subset string '{}' returned no events"
                                               .format(self.subset))

            data = experiment.data[self.channel].astype(float)
            group_values = experiment.data[self.by]
            unit = "event"

        data = pd.Series(np.asarray(data, dtype = float))
        group_values = pd.Series(np.asarray(group_values))

        nan = data.isna() | np.isinf(data)
        if nan.any():
            warnings.warn("Dropping {} NaN or infinite values".format(int(nan.sum())),
                          util.CytoflowOpWarning)
            data = data[~nan]
            group_values = group_values[~nan]

        groups = OrderedDict()
        for g in sorted(group_values.unique()):
            groups[g] = data[group_values == g].to_numpy()

        if len(groups) < 2:
            raise util.CytoflowOpError('by',
                                       "Need at least two groups to compare, but "
                                       "{} has only {} value(s)"
                                       .format(self.by, len(groups)))

        if self.test in _TWO_GROUP_TESTS and len(groups) != 2:
            raise util.CytoflowOpError('test',
                                       "{} compares exactly two groups, but {} has "
                                       "{} values.  Use {} or {} instead, or use "
                                       "'subset' to choose two groups."
                                       .format(TESTS[self.test], self.by, len(groups),
                                               TESTS['anova'], TESTS['kruskal']))

        for g, v in groups.items():
            if len(v) < 2:
                if self.observations == "statistic":
                    raise util.CytoflowOpError('statistic',
                                               "Group {} = {} has only {} {}; each group "
                                               "needs at least two replicates.  Was "
                                               "statistic {} computed with only '{}' in its 'by'?  "
                                               "Compute it by the replicate condition too "
                                               "(eg, the well or tube)."
                                               .format(self.by, g, len(v), unit,
                                                       self.statistic, self.by))
                else:
                    raise util.CytoflowOpError('by',
                                               "Group {} = {} has only {} {}; each group "
                                               "needs at least two."
                                               .format(self.by, g, len(v), unit))

        return groups

    def default_view(self, **kwargs):
        """
        Returns a `StatisticalTestView` that plots the group means +/- SEM
        with significance brackets.
        """
        v = StatisticalTestView(op = self)
        v.trait_set(**kwargs)
        return v


class StatisticalTestView(HasStrictTraits):
    """
    Plots the mean +/- SEM of each group compared by a `StatisticalTestOp`,
    the individual replicates (in ``"statistic"`` mode), and significance
    brackets: for every pair of groups that the Tukey HSD test found
    significantly different (``*`` p < 0.05, ``**`` p < 0.01,
    ``***`` p < 0.001), or, for a two-group test, for the overall test.
    The title reports the overall test.

    Call `plot` with the `Experiment` returned by `StatisticalTestOp.apply`.

    Attributes
    ----------
    op : Instance(`StatisticalTestOp`)
        The operation whose results to plot.

    show_ns : Bool (default = False)
        Also draw brackets for comparisons that aren't significant ("ns")?
    """

    id = Constant('cytoflow.view.statisticaltestview')
    friendly_id = Constant("Statistical Test View")

    op = Instance(StatisticalTestOp)
    show_ns = Bool(False)

    def plot(self, experiment, **kwargs):
        """
        Plot the group means and the significance brackets.
        """

        if experiment is None:
            raise util.CytoflowViewError('experiment', "No experiment specified")

        if self.op is None:
            raise util.CytoflowViewError('op', "No operation specified")

        op = self.op

        if not op.name or op.name not in experiment.statistics \
           or op.name + "_Groups" not in experiment.statistics:
            raise util.CytoflowViewError('op',
                                         "Can't find the test results in the "
                                         "experiment.  Did you call apply()?")

        overall = experiment.statistics[op.name]
        desc = experiment.statistics[op.name + "_Groups"]
        tukey = experiment.statistics.get(op.name + "_Tukey", None) \
                if op.posthoc == "tukey" else None

        # the replicates, if they're few enough to show
        obs = None
        if op.observations == "statistic":
            try:
                obs = op._get_groups(experiment)
            except util.CytoflowOpError:
                obs = None

        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)

        keys = list(desc.index)
        x = np.arange(len(keys))
        means = desc["mean"].to_numpy(dtype = float)
        sems = np.nan_to_num(desc["sem"].to_numpy(dtype = float))

        ax.bar(x, means, width = 0.6, color = "#a6c8e8", edgecolor = "#3a6e9e",
               zorder = 1)
        ax.errorbar(x, means, yerr = sems, fmt = 'none', ecolor = "black",
                    capsize = 5, zorder = 3)

        tops = means + sems
        if obs is not None:
            for xi, k in zip(x, keys):
                if k not in obs:
                    continue
                v = obs[k]
                jitter = np.linspace(-0.12, 0.12, len(v)) if len(v) > 1 else [0.0]
                ax.scatter(xi + np.asarray(jitter), v, color = "black", s = 16,
                           zorder = 4)
                tops[xi] = max(tops[xi], np.max(v))

        ax.set_xticks(x)
        ax.set_xticklabels([_fmt_value(k) for k in keys])
        ax.set_xlabel(op.by)
        ax.set_ylabel(op.feature if op.observations == "statistic" else op.channel)

        # the brackets
        brackets = []
        if tukey is not None:
            for (g1, g2), row in tukey.iterrows():
                if not row["reject"] and not self.show_ns:
                    continue
                brackets.append((keys.index(g1), keys.index(g2), row["significance"]))
        elif op.test in _TWO_GROUP_TESTS and len(keys) == 2:
            p = overall["p_value"].iloc[0]
            stars = significance_stars(p, op.alpha)
            if stars != "ns" or self.show_ns:
                brackets.append((0, 1, stars))

        lo = min(0.0, np.nanmin(means - sems))
        hi = np.nanmax(tops)
        span = (hi - lo) if hi > lo else 1.0
        step = 0.08 * span
        level_top = hi + step * 0.5

        # narrow brackets first, so they nest under the wide ones
        for i, j, text in sorted(brackets, key = lambda b: (abs(b[1] - b[0]), b[0])):
            i, j = min(i, j), max(i, j)
            y = max(level_top, np.nanmax(tops[i:j + 1]) + step * 0.5)
            ax.plot([x[i], x[i], x[j], x[j]], [y, y + step * 0.3, y + step * 0.3, y],
                    color = "black", linewidth = 1, zorder = 5)
            ax.text((x[i] + x[j]) / 2.0, y + step * 0.35, text,
                    ha = "center", va = "bottom", zorder = 5)
            level_top = y + step * 1.2

        ax.set_ylim(top = max(ax.get_ylim()[1], level_top + step * 0.5))

        # the title
        row = overall.iloc[0]
        stat_name = {"anova" : "F", "kruskal" : "H",
                     "welch" : "t", "mannwhitney" : "U"}[op.test]
        if op.test == "anova":
            stat_str = "F({:g}, {:g}) = {:.3g}".format(row["df_between"], row["df_within"], row["statistic"])
        elif op.test == "welch":
            stat_str = "t({:.3g}) = {:.3g}".format(row["df_between"], row["statistic"])
        else:
            stat_str = "{} = {:.3g}".format(stat_name, row["statistic"])
        title = "{}: {}, p = {:.3g}".format(overall.index[0], stat_str, row["p_value"])
        if tukey is not None:
            title += "\nbrackets: Tukey HSD (adjusted p)"
        ax.set_title(title)


def _fmt_value(v):
    try:
        return "{:g}".format(v)
    except (TypeError, ValueError):
        return str(v)
