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
cytoflow.operations.phenograph
------------------------------

Graph-based clustering in the style of PhenoGraph: build a k-nearest-neighbor
graph of the events, weight its edges by how many neighbors the two events
share, then find communities in the graph. `phenograph` has three classes:

`PhenographOp` -- the `IOperation` to perform the clustering.

`Phenograph1DView` -- a diagnostic view of the clustering (1D, using a histogram)

`Phenograph2DView` -- a diagnostic view of the clustering (2D, using a scatterplot)
"""

import random
from warnings import warn
from traits.api import (HasStrictTraits, Str, Dict, Any, Instance,
                        Constant, List, Enum, provides)

import numpy as np
import scipy.sparse
import sklearn.neighbors
import igraph

import pandas as pd

from cytoflow.views import IView, HistogramView, ScatterplotView
import cytoflow.utility as util

from .i_operation import IOperation
from .base_op_views import By1DView, By2DView, AnnotatingView

# how many rows to process at once when computing shared neighbors and when
# assigning events to clusters.  keeps memory use bounded.
_CHUNK = 2000

# above this many events in one graph, things get slow and memory-hungry.
_WARN_EVENTS = 200000

def _jaccard_graph(x, num_neighbors):
    """
    Build the PhenoGraph graph: a k-nearest-neighbors graph whose edges are
    weighted by the Jaccard index of the two events' neighborhoods.

    Parameters
    ----------
    x : numpy.ndarray, shape (n_events, n_channels)
        The (scaled) events

    num_neighbors : int
        The number of neighbors of each event.  Must be less than the number
        of events.

    Returns
    -------
    (edges, weights, nn) : (numpy.ndarray, numpy.ndarray, sklearn.neighbors.NearestNeighbors)
        The edges (shape ``(n_edges, 2)``, each undirected edge listed once)
        and their weights, and the fitted nearest-neighbors model.
    """

    n = x.shape[0]
    k = num_neighbors

    nn = sklearn.neighbors.NearestNeighbors(n_neighbors = k).fit(x)

    # with no arguments, an event is not its own neighbor.
    neighbors = nn.kneighbors(return_distance = False)

    # adj[i, j] == 1 if j is one of i's nearest neighbors
    adj = scipy.sparse.csr_matrix((np.ones(n * k, dtype = np.int32),
                                   neighbors.ravel(),
                                   np.arange(0, n * k + 1, k)),
                                  shape = (n, n))
    adj_t = adj.T.tocsr()

    # the number of shared neighbors of i and j is (adj @ adj.T)[i, j]. we
    # only want it where j is a neighbor of i -- and adj @ adj.T is MUCH
    # denser than adj, so compute it a few rows at a time.
    chunks = []
    for start in range(0, n, _CHUNK):
        rows = adj[start:start + _CHUNK]
        shared = (rows @ adj_t).multiply(rows).tocoo()

        # Jaccard index: |intersection| / |union|, and each neighborhood
        # has exactly k events in it.
        weights = shared.data / (2.0 * k - shared.data)
        chunks.append((shared.row + start, shared.col, weights))

    row = np.concatenate([c[0] for c in chunks])
    col = np.concatenate([c[1] for c in chunks])
    weights = np.concatenate([c[2] for c in chunks])

    # symmetrize the same way PhenoGraph does: an edge that is in the kNN graph
    # in only one direction gets half the weight.  edges between events that
    # share no neighbors are dropped.
    jaccard = scipy.sparse.csr_matrix((weights, (row, col)), shape = (n, n))
    jaccard = scipy.sparse.triu((jaccard + jaccard.T) / 2.0, k = 1).tocoo()

    keep = jaccard.data > 0
    edges = np.column_stack((jaccard.row[keep], jaccard.col[keep]))

    return edges, jaccard.data[keep], nn


def _find_communities(num_vertices, edges, weights, method, resolution, random_state):
    """
    Partition the graph by optimizing modularity.  Returns an array with one
    cluster label per vertex; labels are ordered by cluster size, largest first.
    """

    graph = igraph.Graph(n = num_vertices, edges = edges.tolist(), directed = False)
    weights = weights.tolist()

    # igraph uses python's global random number generator unless told
    # otherwise.  don't disturb (or depend on) the global generator's state.
    if random_state is not None:
        igraph.set_random_number_generator(random.Random(random_state))

    try:
        if method == "leiden":
            partition = graph.community_leiden(objective_function = "modularity",
                                               weights = weights,
                                               resolution = resolution,
                                               n_iterations = -1)
        else:
            partition = graph.community_multilevel(weights = weights,
                                                   resolution = resolution)
    finally:
        if random_state is not None:
            igraph.set_random_number_generator(random)

    membership = np.asarray(partition.membership, dtype = int)

    # number the clusters by size, largest first.  break ties with the
    # position of each cluster's first event so that the numbering is stable.
    _, first, counts = np.unique(membership, return_index = True, return_counts = True)
    order = np.lexsort((first, -counts))
    relabel = np.empty(len(order), dtype = int)
    relabel[order] = np.arange(len(order))

    return relabel[membership]


@provides(IOperation)
class PhenographOp(HasStrictTraits):
    """
    Use graph-based clustering to cluster events, in the style of PhenoGraph
    (Levine et al, 2015).  Unlike `KMeansOp`, you don't have to say how many
    clusters there are -- the algorithm finds out -- and the clusters don't
    need to be round, or all the same size.  It is particularly good at
    finding rare populations in high-dimensional data.

    The algorithm has three steps:

    1. Find each event's `num_neighbors` nearest neighbors (in the scaled
       channels), and make a graph where each event is connected to its
       neighbors.

    2. Weight each connection by how similar the two events' neighborhoods are:
       the number of neighbors they share divided by the number of neighbors
       that either has (the Jaccard index.)

    3. Split the graph into communities -- groups of events that are much more
       strongly connected to each other than to the rest of the graph -- using
       the Leiden or Louvain algorithm to optimize the graph's modularity.

    Call `estimate` to build the graph and find the clusters. Because this
    can be slow and memory-intensive for large data sets, you can ask
    `estimate` to use only a random subset of the events with `sample`.

    Calling `apply` creates a new categorical metadata variable
    named ``{name}``, with possible values ``0`` .... ``n - 1``
    where ``n`` is the number of clusters that were found. Clusters are
    numbered by size: cluster ``0`` has the most events in it (of the events
    used by `estimate`), cluster ``1`` the next most, and so on.

    `apply` assigns **every** event to a cluster, not just the ones that
    were used to build the graph. To do so, it finds the event's
    `num_neighbors` nearest neighbors among the events that `estimate`
    clustered, and assigns the event to the cluster that most of those
    neighbors belong to (ties go to the larger cluster.) An event that is
    identical to one of the clustered events is assigned to that event's
    cluster, so if `sample` is ``1.0`` and you apply the operation to the same
    data that you estimated it with, you get precisely the graph's clusters.

    The same model may not be appropriate for different subsets of the data set.
    If this is the case, you can use the `by` attribute to specify
    metadata by which to aggregate the data before estimating (and applying) a
    model.  Each subset is clustered separately, so **the number of clusters
    may differ from subset to subset, and cluster 0 in one subset has nothing
    to do with cluster 0 in another.**

    Attributes
    ----------
    name : Str
        The operation name; determines the name of the new metadata column

    channels : List(Str)
        The channels to apply the clustering algorithm to.

    scale : Dict(Str : {"linear", "logicle", "log"})
        Re-scale the data in the specified channels before fitting.  If a
        channel is in `channels` but not in `scale`, the current
        package-wide default (set with `set_default_scale`) is used.

        .. note::
           Sometimes you may see events labeled ``-1`` -- this results
           from events for which the selected scale is invalid. For example, if
           an event has a negative measurement in a channel and that channel's
           scale is set to "log", this event will be set to ``-1``.

    num_neighbors : Int (default = 30)
        How many nearest neighbors to connect each event to. Smaller values
        find smaller, more numerous clusters (and follow noise more closely);
        larger values find fewer, bigger clusters and take more time and
        memory.  If a subset has fewer events than this, all of them are
        neighbors.

    resolution : Float (default = 1.0)
        The resolution of the community detection. Values above 1 give
        more, smaller clusters; values below 1 give fewer, larger ones.

    method : {"leiden", "louvain"} (default = "leiden")
        Which algorithm to use to find communities in the graph.  Leiden
        is newer and generally finds better-connected communities; Louvain
        is what the original PhenoGraph used.

    random_state : Int or None (default = 0)
        The seed for the random number generators (both for choosing the
        `sample` and for the community detection.)  With an integer,
        the results are reproducible. If ``None``, they can change from run
        to run.

    sample : Float (default = 1.0)
        What proportion of the events (in each subset in `by`) to build the
        graph with?  Must be greater than 0 and at most 1.  Reduce this if
        `estimate` is too slow or uses too much memory -- a few tens of
        thousands of events is usually plenty.

    by : List(Str)
        A list of metadata attributes to aggregate the data before estimating
        the model.  For example, if the experiment has two pieces of metadata,
        ``Time`` and ``Dox``, setting `by` to ``["Time", "Dox"]`` will
        fit the model separately to each subset of the data with a unique
        combination of ``Time`` and ``Dox``.

    Statistics
    ----------

    Adds a statistic whose name is the name of the operation, whose columns are
    the channels used for clustering, and whose values are the centroids for
    each cluster in that channel (the mean, on the channel's scale, of the
    events that `estimate` put in the cluster). Useful for hierarchical
    clustering, minimum spanning tree visualizations, etc. The index has
    levels from `by`, plus a new level with the same name as the operation.

    The new statistic also has a feature named ``Proportion``, which has the
    proportion of events in each cluster.

    Notes
    -----

    The graph is stored sparsely and the shared-neighbor computation is done
    in pieces, so memory use is proportional to (number of events) x
    `num_neighbors`: 50,000 events with 30 neighbors needs a few hundred
    megabytes.

    References
    ----------

    .. [1] Levine JH, Simonds EF, Bendall SC, et al.
       Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells
       that Correlate with Prognosis.
       Cell 2015; 162(1):184-197.  doi:10.1016/j.cell.2015.05.047

    .. [2] Traag VA, Waltman L, van Eck NJ.
       From Louvain to Leiden: guaranteeing well-connected communities.
       Sci Rep 2019; 9:5233.  doi:10.1038/s41598-019-41695-z

    .. [3] Blondel VD, Guillaume J-L, Lambiotte R, Lefebvre E.
       Fast unfolding of communities in large networks.
       J Stat Mech 2008; P10008.

    Examples
    --------

    .. plot::
        :context: close-figs

        Make a little data set.

        >>> import cytoflow as flow
        >>> import_op = flow.ImportOp()
        >>> import_op.tubes = [flow.Tube(file = "Plate01/RFP_Well_A3.fcs",
        ...                              conditions = {'Dox' : 10.0}),
        ...                    flow.Tube(file = "Plate01/CFP_Well_A4.fcs",
        ...                              conditions = {'Dox' : 1.0})]
        >>> import_op.conditions = {'Dox' : 'float'}
        >>> ex = import_op.apply()

    Create and parameterize the operation.

    .. plot::
        :context: close-figs

        >>> pg_op = flow.PhenographOp(name = 'Phenograph',
        ...                           channels = ['V2-A', 'Y2-A'],
        ...                           scale = {'V2-A' : 'log',
        ...                                    'Y2-A' : 'log'},
        ...                           resolution = 0.1)

    Estimate the clusters

    .. plot::
        :context: close-figs

        >>> pg_op.estimate(ex)

    Plot a diagnostic view

    .. plot::
        :context: close-figs

        >>> pg_op.default_view().plot(ex)

    Apply the gate

    .. plot::
        :context: close-figs

        >>> ex2 = pg_op.apply(ex)

    Plot a diagnostic view with the event assignments

    .. plot::
        :context: close-figs

        >>> pg_op.default_view().plot(ex2)
    """

    id = Constant('cytoflow.operations.phenograph')
    friendly_id = Constant("Graph Clustering (PhenoGraph)")

    name = Str
    channels = List(Str)
    scale = Dict(Str, util.ScaleEnum)
    num_neighbors = util.PositiveInt(30, allow_zero = False)
    resolution = util.PositiveFloat(1.0, allow_zero = False)
    method = Enum("leiden", "louvain")
    random_state = util.IntOrNone(0)
    sample = util.UnitFloat(1.0)
    by = List(Str)

    # for each group: the nearest-neighbors model of the events that were
    # clustered, those events' clusters, and the cluster centroids (in the
    # scaled space; shape (n_clusters, n_channels))
    _nn = Dict(Any, Instance(sklearn.neighbors.NearestNeighbors), transient = True)
    _labels = Dict(Any, Instance(np.ndarray), transient = True)
    _centers = Dict(Any, Instance(np.ndarray), transient = True)
    _scale = Dict(Str, Instance(util.IScale), transient = True)

    def estimate(self, experiment, subset = None):
        """
        Build the nearest-neighbors graph and find the clusters.

        Parameters
        ----------
        experiment : Experiment
            The `Experiment` to use to estimate the clusters

        subset : str (default = None)
            A Python expression that specifies a subset of the data in
            ``experiment`` to use to parameterize the operation.
        """

        if experiment is None:
            raise util.CytoflowOpError('experiment',
                                       "No experiment specified")

        if self.num_neighbors < 2:
            raise util.CytoflowOpError('num_neighbors',
                                       "num_neighbors must be >= 2")

        if self.sample <= 0:
            raise util.CytoflowOpError('sample',
                                       "sample must be > 0")

        if len(self.channels) == 0:
            raise util.CytoflowOpError('channels',
                                       "Must set at least one channel")

        if len(self.channels) != len(set(self.channels)):
            raise util.CytoflowOpError('channels',
                                       "Must not duplicate channels")

        for c in self.channels:
            if c not in experiment.data:
                raise util.CytoflowOpError('channels',
                                           "Channel {0} not found in the experiment"
                                      .format(c))

        for c in self.scale:
            if c not in self.channels:
                raise util.CytoflowOpError('scale',
                                           "Scale set for channel {0}, but it isn't "
                                           "in the experiment"
                                           .format(c))

        for b in self.by:
            if b not in experiment.data:
                raise util.CytoflowOpError('by',
                                           "Aggregation metadata {} not found, "
                                           "must be one of {}"
                                           .format(b, experiment.conditions))

        if subset:
            try:
                experiment = experiment.query(subset)
            except:
                raise util.CytoflowOpError('subset',
                                            "Subset string '{0}' isn't valid"
                                            .format(subset))

            if len(experiment) == 0:
                raise util.CytoflowOpError('subset',
                                           "Subset string '{0}' returned no events"
                                           .format(subset))

        if self.by:
            groupby = experiment.data.groupby(self.by, observed = False)
        else:
            # use a lambda expression to return a group that contains
            # all the events
            groupby = experiment.data.groupby(lambda _: True, observed = False)

        # get the scale. estimate the scale params for the ENTIRE data set,
        # not subsets we get from groupby().  And we need to save it so that
        # the data is transformed the same way when we apply()
        scale = {}
        for c in self.channels:
            if c in self.scale:
                scale[c] = util.scale_factory(self.scale[c], experiment, channel = c)
            else:
                scale[c] = util.scale_factory(util.get_default_scale(), experiment, channel = c)

        nn = {}
        labels = {}
        centers = {}
        for group, data_subset in groupby:
            if len(data_subset) == 0:
                warn("Group {} had no data".format(group),
                     util.CytoflowOpWarning)
                continue

            if self.sample < 1.0:
                data_subset = data_subset.sample(frac = self.sample,
                                                 random_state = self.random_state)

            x = data_subset.loc[:, self.channels[:]]
            for c in self.channels:
                x[c] = scale[c](x[c])

            # drop data that isn't in the scale range
            for c in self.channels:
                x = x[~(np.isnan(x[c]))]
            x = x.values

            if x.shape[0] < 3:
                raise util.CytoflowOpError('sample',
                                           "Group {} has fewer than three events "
                                           "to cluster{}"
                                           .format(group,
                                                   " -- try increasing 'sample'"
                                                   if self.sample < 1.0 else ""))

            if x.shape[0] > _WARN_EVENTS:
                warn("Building a graph with {} events, which may be slow and "
                     "use a lot of memory. Consider decreasing 'sample'."
                     .format(x.shape[0]),
                     util.CytoflowOpWarning)

            edges, weights, nn[group] = \
                _jaccard_graph(x, min(self.num_neighbors, x.shape[0] - 1))

            labels[group] = l = _find_communities(x.shape[0],
                                                  edges,
                                                  weights,
                                                  self.method,
                                                  self.resolution,
                                                  self.random_state)

            centers[group] = np.vstack([x[l == i].mean(axis = 0)
                                        for i in range(l.max() + 1)])

        # do this last so the UI can pick up that the estimate changed
        self._scale = scale
        self._nn = nn
        self._labels = labels
        self._centers = centers

    def _predict(self, group, x):
        """
        Assign the (scaled) events in ``x`` to the clusters that were
        estimated for ``group``, by majority vote of their nearest neighbors.
        """
        nn = self._nn[group]
        labels = self._labels[group]
        num_clusters = self._centers[group].shape[0]

        predicted = np.empty(x.shape[0], dtype = int)
        for start in range(0, x.shape[0], _CHUNK):
            dist, idx = nn.kneighbors(x[start:start + _CHUNK])
            neighbor_labels = labels[idx]

            # count the votes for each cluster.  argmax returns the first
            # maximum, so ties go to the lowest-numbered (ie, largest) cluster
            n = neighbor_labels.shape[0]
            votes = np.bincount((neighbor_labels +
                                 np.arange(n)[:, np.newaxis] * num_clusters).ravel(),
                                minlength = n * num_clusters).reshape((n, num_clusters))
            p = votes.argmax(axis = 1)

            # an event that was in the graph keeps the cluster it got there.
            same = dist[:, 0] == 0
            p[same] = neighbor_labels[same, 0]

            predicted[start:start + _CHUNK] = p

        return predicted

    def apply(self, experiment):
        """
        Assign events to the clusters found by `estimate`.

        Returns
        -------
        Experiment
            a new Experiment with one additional entry in `Experiment.conditions`
            named `name`, of type ``category``.  The new category has
            values  ``0``, ``1``, etc to indicate which cluster
            an event is a member of.

            The new `Experiment` also has one new statistic named `name`,
            which has the centroids of each cluster and the proportion of
            events in it.
        """

        if experiment is None:
            raise util.CytoflowOpError('experiment',
                                       "No experiment specified")

        # make sure name got set!
        if not self.name:
            raise util.CytoflowOpError('name',
                                       "You have to set the gate's name "
                                       "before applying it!")

        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError('name',
                                       "Name can only contain letters, numbers and underscores."
                                       .format(self.name))

        if self.name in experiment.data.columns:
            raise util.CytoflowOpError('name',
                                       "Experiment already has a column named {0}"
                                       .format(self.name))

        if not self._centers:
            raise util.CytoflowOpError(None,
                                       "No clusters found.  Did you forget to "
                                       "call estimate()?")

        if len(self.channels) == 0:
            raise util.CytoflowOpError('channels',
                                       "Must set at least one channel")

        for c in self.channels:
            if c not in experiment.data:
                raise util.CytoflowOpError('channels',
                                           "Channel {0} not found in the experiment"
                                      .format(c))

        for c in self.scale:
            if c not in self.channels:
                raise util.CytoflowOpError('scale',
                                           "Scale set for channel {0}, but it isn't "
                                           "in the experiment"
                                           .format(c))

        for b in self.by:
            if b not in experiment.data:
                raise util.CytoflowOpError('by',
                                           "Aggregation metadata {} not found, "
                                           "must be one of {}"
                                           .format(b, experiment.conditions))

        if self.name in self.by:
            raise util.CytoflowOpError('by',
                                       "'{}' is going to be added as an "
                                       "index level to the new statistic, so you "
                                       "can't use it to aggregate events."
                                       .format(self.name))


        if self.by:
            groupby = experiment.data.groupby(self.by, observed = False)
        else:
            # use a lambda expression to return a group that contains
            # all the events
            groupby = experiment.data.groupby(lambda _: True, observed = False)

        event_assignments = pd.Series([-1] * len(experiment), dtype = "int")

        # the rows of the new statistic.  unlike KMeans, the number of clusters
        # can be different in each group.
        stat_idx = []
        stat_rows = []

        for group, data_subset in groupby:
            if len(data_subset) == 0:
                warn("Group {} had no data".format(group),
                     util.CytoflowOpWarning)
                continue

            if group not in self._centers:
                raise util.CytoflowOpError('by',
                                           "Group {} not found in the estimated model. "
                                           "Do you need to re-run estimate()?"
                                           .format(group))

            x = data_subset.loc[:, self.channels[:]]
            for c in self.channels:
                x[c] = self._scale[c](x[c])

            # which values are missing?

            x_na = pd.Series([False] * len(x))
            for c in self.channels:
                x_na[np.isnan(x[c]).values] = True

            x = x.values
            x_na = x_na.values
            group_idx = data_subset.index

            centers = self._centers[group]

            predicted = pd.Series([-1] * len(x))
            predicted[~x_na] = self._predict(group, x[~x_na])
            predicted.index = group_idx

            event_assignments.iloc[group_idx] = predicted

            for c in range(centers.shape[0]):
                if len(self.by) == 0:
                    g = tuple([c])
                elif isinstance(group, tuple):
                    g = group + tuple([c])
                else:
                    g = tuple([group, c])

                stat_idx.append(g)
                stat_rows.append(list(centers[c]) + [sum(predicted == c) / len(predicted)])

        idx = pd.MultiIndex.from_tuples(stat_idx, names = list(self.by) + [self.name])
        new_stat = pd.DataFrame(stat_rows,
                                index = idx,
                                columns = list(self.channels) + ["Proportion"],
                                dtype = 'float').sort_index()

        new_experiment = experiment.clone(deep = False)
        new_experiment.add_condition(self.name, "category", event_assignments)

        new_experiment.statistics[self.name] = new_stat

        new_experiment.history.append(self.clone_traits(transient = lambda _: True))
        return new_experiment

    def default_view(self, **kwargs):
        """
        Returns a diagnostic plot of the clustering.

        Returns
        -------
            IView : an IView, call `Phenograph1DView.plot` to see the diagnostic plot.
        """
        channels = kwargs.pop('channels', self.channels)
        scale = kwargs.pop('scale', self.scale)

        for c in channels:
            if c not in self.channels:
                raise util.CytoflowViewError('channels',
                                             "Channel {} isn't in the operation's channels"
                                             .format(c))

        for s in scale:
            if s not in self.channels:
                raise util.CytoflowViewError('scale',
                                             "Channel {} isn't in the operation's channels"
                                             .format(s))

        for c in channels:
            if c not in scale:
                scale[c] = util.get_default_scale()

        if len(channels) == 0:
            raise util.CytoflowViewError('channels',
                                         "Must specify at least one channel for a default view")
        elif len(channels) == 1:
            v = Phenograph1DView(op = self)
            v.trait_set(channel = channels[0],
                        scale = scale[channels[0]],
                        **kwargs)
            return v

        elif len(channels) == 2:
            v = Phenograph2DView(op = self)
            v.trait_set(xchannel = channels[0],
                        ychannel = channels[1],
                        xscale = scale[channels[0]],
                        yscale = scale[channels[1]],
                        **kwargs)
            return v

        else:
            raise util.CytoflowViewError('channels',
                                         "Can't specify more than two channels for a default view")


@provides(IView)
class Phenograph1DView(By1DView, AnnotatingView, HistogramView):
    """
    A diagnostic view for `PhenographOp` (1D, using a histogram)
    
    Attributes
    ----------    
    op : Instance(PhenographOp)
        The op whose parameters we're viewing.
    """
    
    id = Constant('cytoflow.views.phenograph1dview')
    friendly_id = Constant("1D Phenograph Diagnostic Plot")
    
    channel = Str
    scale = util.ScaleEnum
    
    def plot(self, experiment, **kwargs):
        """
        Plot the plots.
        
        Parameters
        ----------
        """
        
        if experiment is None:
            raise util.CytoflowViewError('experiment', "No experiment specified")
                
        view, trait_name = self._strip_trait(self.op.name)
        
        
        if self.channel in self.op._scale:
            scale = self.op._scale[self.channel]
        else:
            scale = util.scale_factory(self.scale, experiment, channel = self.channel)
    
        super(Phenograph1DView, view).plot(experiment,
                                       annotation_facet = self.op.name,
                                       annotation_trait = trait_name,
                                       annotations = self.op._centers,
                                       scale = scale,
                                       **kwargs)
 
    def _annotation_plot(self, axes, annotation, annotation_facet, 
                         annotation_value, annotation_color, **kwargs):
                                                        
        # plot the cluster centers
            
        centers = annotation
        
        kwargs.setdefault('orientation', 'vertical')
        
        if kwargs['orientation'] == 'horizontal':
            scale = kwargs['yscale']
            cidx = self.op.channels.index(self.channel)
            for k in range(centers.shape[0]):
                c = scale.inverse(centers[k][cidx])
                axes.axhline(c, linewidth=3, color='blue')         
        else:
            scale = kwargs['xscale']
            cidx = self.op.channels.index(self.channel)
            for k in range(centers.shape[0]):
                c = scale.inverse(centers[k][cidx])
                axes.axvline(c, linewidth=3, color='blue')                      

     
@provides(IView)
class Phenograph2DView(By2DView, AnnotatingView, ScatterplotView):
    """
    A diagnostic view for `PhenographOp` (2D, using a scatterplot).
    
    Attributes
    ----------
    op : Instance(PhenographOp)
        The op whose parameters we're viewing.        
    """
     
    id = Constant('cytoflow.view.phenograph2dview')
    friendly_id = Constant("2D Phenograph Diagnostic Plot")
    
    xchannel = Str
    ychannel = Str
    
    xscale = util.ScaleEnum
    yscale = util.ScaleEnum
    
    def plot(self, experiment, **kwargs):
        """
        Plot the plots.
        
        Parameters
        ----------
        """

        if experiment is None:
            raise util.CytoflowViewError('experiment', "No experiment specified")
                
        view, trait_name = self._strip_trait(self.op.name)
        
        if self.xchannel in self.op._scale:
            xscale = self.op._scale[self.xchannel]
        else:
            xscale = util.scale_factory(self.xscale, experiment, channel = self.xchannel)

        if self.ychannel in self.op._scale:
            yscale = self.op._scale[self.ychannel]
        else:
            yscale = util.scale_factory(self.yscale, experiment, channel = self.ychannel)
    
        super(Phenograph2DView, view).plot(experiment,
                                       annotation_facet = self.op.name,
                                       annotation_trait = trait_name,
                                       annotations = self.op._centers,
                                       xscale = xscale,
                                       yscale = yscale,
                                       **kwargs)
 
    def _annotation_plot(self, axes, annotation, annotation_facet, 
                         annotation_value, annotation_color, **kwargs):
                                                        
        # plot the cluster centers
            
        centers = annotation
        xscale = kwargs['xscale']
        yscale = kwargs['yscale']
        
        ix = self.op.channels.index(self.xchannel)
        iy = self.op.channels.index(self.ychannel)
        
        for k in range(centers.shape[0]):
            x = xscale.inverse(centers[k][ix])
            y = yscale.inverse(centers[k][iy])
            
            axes.plot(x, y, '*', color = 'blue')

util.expand_class_attributes(Phenograph1DView)
util.expand_method_parameters(Phenograph1DView, Phenograph1DView.plot)

util.expand_class_attributes(Phenograph2DView)
util.expand_method_parameters(Phenograph2DView, Phenograph2DView.plot)