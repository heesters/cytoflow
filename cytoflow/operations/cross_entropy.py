#!/usr/bin/env python3
# coding: latin-1

# (c) Brian Teague 2018-2026
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
cytoflow.operations.cross_entropy
---------------------------------

The *cross entropy test*: a statistical test of whether groups of events
(eg, samples or treatments) differ in a t-SNE embedding.  `cross_entropy`
has:

`CrossEntropyTestOp` -- the `IOperation` that does the test.

`pointwise_cross_entropy` -- computes the cross entropy of each event.

`CrossEntropyTestView` -- plots the distributions of cross entropy (or a
dendrogram of the groups.)

This is an independent implementation of the method published in [1]_,
written from the equations in the paper.  It is not derived from (and may
differ numerically from) the authors' R scripts.

References
----------
.. [1] Roca CP, Burton OT, Neumann J, et al.  A cross entropy test allows
   quantitative statistical comparison of t-SNE and UMAP representations.
   Cell Reports Methods 3, 100390 (2023).
   https://doi.org/10.1016/j.crmeth.2022.100390
"""

import itertools
import warnings

from traits.api import (HasStrictTraits, Str, List, Dict, Float, Int, Enum,
                        Any, Bool, Constant, Instance, provides)

import numpy as np
import pandas as pd
import scipy.sparse
import scipy.stats
import scipy.spatial.distance
import scipy.cluster.hierarchy
import matplotlib.pyplot as plt

import cytoflow.views
import cytoflow.utility as util

from .i_operation import IOperation


def _conditional_probabilities(sq_dists, perplexity, tol = 1e-5, max_iter = 100):
    """
    For each row of squared distances to a point's neighbors, find the
    Gaussian kernel width that gives the conditional distribution
    p_{j|i} the requested perplexity (by bisection, as in t-SNE), and
    return those distributions (rows sum to 1.)
    """
    n, k = sq_dists.shape
    target = np.log(perplexity)

    # work with distances relative to the nearest neighbor, for stability
    d = sq_dists - sq_dists.min(axis = 1, keepdims = True)

    beta = np.ones(n)
    lo = np.zeros(n)
    hi = np.full(n, np.inf)

    for _ in range(max_iter):
        p = np.exp(-d * beta[:, None])
        s = p.sum(axis = 1)
        p /= s[:, None]
        # the Shannon entropy of each row
        h = -(p * np.log(np.maximum(p, 1e-300))).sum(axis = 1)
        diff = h - target
        if np.all(np.abs(diff) < tol):
            break

        # entropy too high -> narrower kernel (larger beta)
        too_high = diff > 0
        lo = np.where(too_high, beta, lo)
        hi = np.where(~too_high, beta, hi)
        beta = np.where(np.isinf(hi), beta * 2, (lo + hi) / 2)

    return p


def pointwise_cross_entropy(high, low, perplexity = 30.0, num_neighbors = None):
    """
    The cross entropy of each event between the neighbor probabilities in
    the original (high-dimensional) space and in the embedding, as defined
    for t-SNE in Roca et al. (2023):

    - :math:`p_{j|i}` -- a Gaussian kernel on the high-dimensional
      distances from event *i*, with the width chosen so that the
      distribution has the given ``perplexity``;
    - :math:`p_{ij} = (p_{j|i} + p_{i|j}) / 2n`;
    - :math:`q_{ij} = (1 + |y_i - y_j|^2)^{-1} / \\sum_{k \\ne l} (1 + |y_k - y_l|^2)^{-1}`
      in the embedding;
    - :math:`H_i = -\\sum_{j \\ne i} n\\,p_{ij} \\log(n\\,q_{ij})`.

    As in most t-SNE implementations, :math:`p_{j|i}` is computed for the
    ``num_neighbors`` nearest neighbors of each event only (it is
    negligible for the others.)

    Parameters
    ----------
    high : numpy.ndarray, (events x channels)
        The (scaled) data that was embedded.

    low : numpy.ndarray, (events x dimensions)
        The embedding (eg, the two t-SNE coordinates.)

    perplexity : float
        Use the perplexity the embedding was made with.

    num_neighbors : int
        Default: ``3 * perplexity``.

    Returns
    -------
    numpy.ndarray
        The cross entropy of each event.
    """
    from sklearn.neighbors import NearestNeighbors

    high = np.asarray(high, dtype = float)
    low = np.asarray(low, dtype = float)
    n = len(high)

    if len(low) != n:
        raise util.CytoflowOpError(None,
                                   "The data and the embedding must have the "
                                   "same number of events")

    if num_neighbors is None:
        num_neighbors = int(3 * perplexity)
    k = min(num_neighbors, n - 1)

    if k < 2 or perplexity >= k:
        raise util.CytoflowOpError('perplexity',
                                   "Too few events ({}) for a perplexity of {}"
                                   .format(n, perplexity))

    nn = NearestNeighbors(n_neighbors = k + 1).fit(high)
    dist, idx = nn.kneighbors(high)

    # drop each point itself.  (with duplicate points, "itself" isn't
    # necessarily the first neighbor.)
    rows = np.arange(n)[:, None]
    is_self = idx == rows
    no_self = ~is_self.any(axis = 1)
    is_self[no_self, -1] = True
    keep = ~is_self
    idx = idx[keep].reshape(n, k)
    dist = dist[keep].reshape(n, k)

    cond = _conditional_probabilities(dist ** 2, perplexity)

    pc = scipy.sparse.csr_matrix((cond.ravel(),
                                  (np.repeat(np.arange(n), k), idx.ravel())),
                                 shape = (n, n))
    p = ((pc + pc.T) / (2.0 * n)).tocoo()

    # the normalization of q: the sum over all pairs, in chunks
    z = 0.0
    chunk = max(1, int(2e7 // n))
    for start in range(0, n, chunk):
        d2 = scipy.spatial.distance.cdist(low[start:start + chunk], low, 'sqeuclidean')
        z += (1.0 / (1.0 + d2)).sum()
    z -= n      # the k == l terms

    d2 = ((low[p.row] - low[p.col]) ** 2).sum(axis = 1)
    q = (1.0 / (1.0 + d2)) / z

    terms = -(n * p.data) * np.log(n * q)
    return np.bincount(p.row, weights = terms, minlength = n)


def _adjust(pvalues, method):
    p = np.asarray(pvalues, dtype = float)
    m = len(p)
    if method == 'none' or m == 0:
        return p
    if method == 'bonferroni':
        return np.minimum(p * m, 1.0)
    # holm
    order = np.argsort(p)
    adj = np.empty(m)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (m - rank) * p[i])
        adj[i] = min(running, 1.0)
    return adj


@provides(IOperation)
class CrossEntropyTestOp(HasStrictTraits):
    """
    The *cross entropy test* (Roca et al., Cell Reports Methods 2023): test
    whether groups of events -- eg, samples, or treated vs control -- are
    distributed differently in a t-SNE embedding.

    Looking at t-SNE plots of two samples and deciding "they're different"
    is subjective.  This test makes it quantitative.  For every event, it
    computes the *cross entropy* between the event's neighborhood in the
    original data and in the embedding (the quantity t-SNE minimizes.)
    Then it compares the distributions of cross entropy between each pair
    of groups with a two-sample Kolmogorov-Smirnov test.

    Embed all the groups **together**, in one t-SNE (`tSNEOp` with no
    ``by``), then run this operation with the same `channels`, `scale` and
    `perplexity`.

    `estimate` does the computation, and `apply` adds these statistics:

    - `name` -- one row for each pair of groups (indexed by ``group1`` and
      ``group2``): ``D``, the Kolmogorov-Smirnov statistic (the largest
      difference between the two cumulative distributions; between 0 and 1,
      and useful as a *distance* between the groups); ``p``; ``p_adj``,
      adjusted for the number of pairs with `correction`; and ``reject``,
      whether ``p_adj < alpha``.

    - `name` + ``_Groups`` -- for each group: the number of events used,
      and the median and mean cross entropy.

    .. note::
        With thousands of events per group, even very small differences
        are "significant".  Look at ``D`` too, and consider a strict
        threshold: the authors suggest p < 0.001 for data sets with many
        events.

    .. note::
        The test compares how each group's events sit within the structure
        of the data.  Groups that differ only in the proportions of 
        populations with a similar shape can look alike to it.
        
    .. note::
        This is the t-SNE form of the test.  It can be run on any
        embedding's coordinates (eg, UMAP's), where it still measures how
        each group's neighborhoods are preserved, but the paper's Methods
        don't define a UMAP-specific cross entropy, so none is implemented.

    Attributes
    ----------
    name : Str
        The operation name; used to name the new statistics.

    channels : List(Str)
        The channels that were embedded.

    scale : Dict(Str : {"linear", "logicle", "log", "arcsinh"})
        The scale of each channel, as used for the embedding.  Channels not
        in the dict are linear.

    embedding : List(Str)
        The channels holding the embedding's coordinates, eg
        ``["tSNE_1", "tSNE_2"]``.

    by : Str
        The condition whose values are the groups to compare.

    perplexity : Float (default = 30)
        The perplexity the t-SNE was made with.

    max_events : Int (default = 20000)
        The computation needs all pairs of events, so at most this many
        events (a random sample, in proportion to each group) are used.

    correction : {"holm", "bonferroni", "none"} (default = "holm")
        The correction for testing multiple pairs of groups.  (The paper
        doesn't specify one.)

    alpha : Float (default = 0.05)
        The significance threshold for ``reject``.

    random_state : Int (default = 0)
        The seed for the random sample of events.

    Examples
    --------

    >>> tsne_op = flow.tSNEOp(name = "tSNE", channels = ["B1-A", "V2-A", "Y2-A"],
    ...                       scale = {c : "logicle" for c in ["B1-A", "V2-A", "Y2-A"]},
    ...                       perplexity = 30)
    >>> tsne_op.estimate(ex)
    >>> ex2 = tsne_op.apply(ex)
    >>> ce_op = flow.CrossEntropyTestOp(name = "CE", by = "Dox",
    ...                                 channels = tsne_op.channels,
    ...                                 scale = tsne_op.scale,
    ...                                 embedding = ["tSNE_1", "tSNE_2"],
    ...                                 perplexity = 30)
    >>> ce_op.estimate(ex2)
    >>> ex3 = ce_op.apply(ex2)
    >>> ex3.statistics["CE"]
    >>> ce_op.default_view().plot(ex3)
    """

    id = Constant('cytoflow.operations.cross_entropy')
    friendly_id = Constant("Cross Entropy Test")

    name = Str
    channels = List(Str)
    scale = Dict(Str, util.ScaleEnum)
    embedding = List(Str)
    by = Str
    perplexity = Float(30.0)
    max_events = Int(20000)
    correction = Enum("holm", "bonferroni", "none")
    alpha = Float(0.05)
    random_state = Int(0)

    # the results of estimate()
    _cross_entropy = Dict(Any, Any, transient = True)   # group -> array
    _pairs = Instance(pd.DataFrame, transient = True)
    _groups = Instance(pd.DataFrame, transient = True)

    def estimate(self, experiment, subset = None):
        """
        Compute each event's cross entropy and test the groups against each
        other.

        Parameters
        ----------
        experiment : Experiment
            The experiment, with the embedding already in it.

        subset : Str
            Only use these events.
        """
        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")

        if len(self.channels) < 2:
            raise util.CytoflowOpError('channels',
                                       "Must set at least two channels")

        if len(self.embedding) < 2:
            raise util.CytoflowOpError('embedding',
                                       "Must set the embedding's (at least "
                                       "two) channels")

        for c in self.channels:
            if c not in experiment.channels:
                raise util.CytoflowOpError('channels',
                                           "Channel {} isn't in the experiment".format(c))

        for c in self.embedding:
            if c not in experiment.channels:
                raise util.CytoflowOpError('embedding',
                                           "Channel {} isn't in the experiment. "
                                           "Did you apply the t-SNE operation?".format(c))
            if c in self.channels:
                raise util.CytoflowOpError('embedding',
                                           "{} can't be both a data channel and "
                                           "an embedding channel".format(c))

        for c in self.scale:
            if c not in self.channels:
                raise util.CytoflowOpError('scale',
                                           "Scale set for channel {}, but it isn't "
                                           "in 'channels'".format(c))

        if not self.by:
            raise util.CytoflowOpError('by', "Must set 'by', the condition to compare")

        if self.by not in experiment.conditions:
            raise util.CytoflowOpError('by',
                                       "Condition {} isn't in the experiment".format(self.by))

        if not self.perplexity > 1:
            raise util.CytoflowOpError('perplexity', "perplexity must be > 1")

        if self.correction not in ("holm", "bonferroni", "none"):
            raise util.CytoflowOpError('correction', "Unknown correction")

        self._cross_entropy = {}
        self._pairs = None
        self._groups = None

        if subset:
            try:
                experiment = experiment.query(subset)
            except Exception as e:
                raise util.CytoflowOpError('subset',
                                           "Subset string '{0}' isn't valid"
                                           .format(subset)) from e
            if len(experiment) == 0:
                raise util.CytoflowOpError('subset',
                                           "Subset string '{0}' returned no events"
                                           .format(subset))

        data = experiment.data[self.channels + self.embedding + [self.by]].copy()
        for c in self.channels:
            s = util.scale_factory(self.scale.get(c, util.get_default_scale()),
                                   experiment, channel = c)
            data[c] = s(data[c])
        data = data.replace([np.inf, -np.inf], np.nan).dropna()

        if self.max_events and len(data) > self.max_events:
            data = data.sample(n = self.max_events, random_state = self.random_state)

        groups = [g for g in sorted(data[self.by].unique())]
        if len(groups) < 2:
            raise util.CytoflowOpError('by',
                                       "Need at least two groups to compare; "
                                       "condition {} has {}"
                                       .format(self.by, len(groups)))

        counts = data[self.by].value_counts()
        if counts.min() < 10:
            raise util.CytoflowOpError('by',
                                       "Group {} has fewer than 10 events"
                                       .format(counts.idxmin()))
        if counts.min() < 3 * self.perplexity:
            warnings.warn("Group {} has only {} events, which is few for a "
                          "perplexity of {}"
                          .format(counts.idxmin(), counts.min(), self.perplexity),
                          util.CytoflowOpWarning)

        try:
            ce = pointwise_cross_entropy(data[self.channels].values,
                                         data[self.embedding].values,
                                         perplexity = self.perplexity)
        except util.CytoflowOpError as e:
            raise util.CytoflowOpError('perplexity', e.args[-1]) from e

        labels = data[self.by].values
        cross_entropy = {g : ce[labels == g] for g in groups}

        rows = []
        for g1, g2 in itertools.combinations(groups, 2):
            res = scipy.stats.ks_2samp(cross_entropy[g1], cross_entropy[g2])
            rows.append((g1, g2, float(res.statistic), float(res.pvalue)))

        pairs = pd.DataFrame(rows, columns = ["group1", "group2", "D", "p"])
        pairs["p_adj"] = _adjust(pairs["p"].values, self.correction)
        pairs["reject"] = pairs["p_adj"] < self.alpha
        pairs = pairs.set_index(["group1", "group2"])

        desc = pd.DataFrame({"n" : [len(cross_entropy[g]) for g in groups],
                             "median" : [float(np.median(cross_entropy[g])) for g in groups],
                             "mean" : [float(np.mean(cross_entropy[g])) for g in groups]},
                            index = pd.Index(groups, name = self.by))

        self._groups = desc
        self._pairs = pairs
        # set last: the GUI watches this one
        self._cross_entropy = cross_entropy

    def apply(self, experiment):
        """
        Add the test's results to the experiment's statistics.

        Returns
        -------
        Experiment
            A new experiment (sharing ``experiment``'s data) with the
            statistics `name` and `name` + ``_Groups``.
        """
        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")

        if not self.name:
            raise util.CytoflowOpError('name', "Must specify a name")

        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError('name',
                                       "Name can only contain letters, numbers and underscores.")

        if self._pairs is None:
            raise util.CytoflowOpError(None, "No results. Did you forget to call estimate()?")

        for n in (self.name, self.name + "_Groups"):
            if n in experiment.statistics:
                raise util.CytoflowOpError('name',
                                           "Experiment already has a statistic named {}".format(n))

        new_experiment = experiment.clone(deep = False)
        new_experiment.statistics[self.name] = self._pairs.copy()
        new_experiment.statistics[self.name + "_Groups"] = self._groups.copy()
        new_experiment.history.append(self.clone_traits(transient = lambda _: True))
        return new_experiment

    def distance_matrix(self):
        """
        The matrix of Kolmogorov-Smirnov distances (``D``) between the
        groups, as a `pandas.DataFrame`.
        """
        if self._pairs is None:
            raise util.CytoflowOpError(None, "No results. Did you forget to call estimate()?")
        groups = list(self._groups.index)
        m = pd.DataFrame(0.0, index = groups, columns = groups)
        for (g1, g2), row in self._pairs.iterrows():
            m.loc[g1, g2] = m.loc[g2, g1] = row["D"]
        return m

    def default_view(self, **kwargs):
        """
        Returns a `CrossEntropyTestView`: the cumulative distribution of
        cross entropy in each group, or (``dendrogram = True``) a
        dendrogram of the groups.
        """
        v = CrossEntropyTestView(op = self)
        v.trait_set(**kwargs)
        return v


@provides(cytoflow.views.IView)
class CrossEntropyTestView(HasStrictTraits):
    """
    Plots the results of a `CrossEntropyTestOp`.

    Attributes
    ----------
    op : Instance(CrossEntropyTestOp)
        The operation whose results to plot.

    dendrogram : Bool (default = False)
        If ``False``, plot the cumulative distribution of the events' cross
        entropy in each group (the curves the Kolmogorov-Smirnov test
        compares: ``D`` is the largest vertical gap between two of them.)
        If ``True``, plot a dendrogram of the groups, clustered (average
        linkage) by ``D``.
    """

    id = Constant('cytoflow.view.crossentropytestview')
    friendly_id = Constant("Cross Entropy Test")

    op = Instance(CrossEntropyTestOp)
    dendrogram = Bool(False)

    def plot(self, experiment = None, **kwargs):
        """
        Plot the view.  ``experiment`` isn't used.

        Parameters
        ----------
        title : str
            The plot title
        """
        if self.op is None or self.op._pairs is None:
            raise util.CytoflowViewError('op',
                                         "No results to plot. Did you forget "
                                         "to call estimate()?")

        plt.figure()
        ax = plt.gca()

        if self.dendrogram:
            m = self.op.distance_matrix()
            if len(m) < 3:
                raise util.CytoflowViewError('dendrogram',
                                             "Need at least three groups for a dendrogram")
            z = scipy.cluster.hierarchy.linkage(
                    scipy.spatial.distance.squareform(m.values, checks = False),
                    method = 'average')
            scipy.cluster.hierarchy.dendrogram(z, labels = [str(g) for g in m.index],
                                               ax = ax, color_threshold = 0,
                                               above_threshold_color = 'k')
            ax.set_ylabel("Kolmogorov-Smirnov distance (D)")
            ax.set_xlabel(self.op.by)
        else:
            for g, ce in self.op._cross_entropy.items():
                x = np.sort(ce)
                ax.step(x, np.arange(1, len(x) + 1) / len(x), where = 'post',
                        label = "{} (n = {})".format(g, len(x)))
            ax.set_xlabel("Cross entropy")
            ax.set_ylabel("Cumulative fraction of events")
            ax.set_ylim(0, 1)
            ax.legend(title = self.op.by)

        if kwargs.get('title'):
            ax.set_title(kwargs['title'])


util.expand_class_attributes(CrossEntropyTestOp)
util.expand_method_parameters(CrossEntropyTestOp, CrossEntropyTestOp.estimate)
