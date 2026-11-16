#!/usr/bin/env python3.11
# coding: latin-1

# (c) Massachusetts Institute of Technology 2015-2018
# (c) Brian Teague 2018-2025
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
cytoflow.operations.umap
------------------------

Apply Uniform Manifold Approximation and Projection (UMAP) to flow data. Like
tSNE, this is a non-linear dimensionality reduction algorithm that tries to
preserve the local structure of the data; unlike tSNE, it is generally
faster and (arguably) better at preserving global structure as well.
`umap` has one class:

`UMAPOp` -- the `IOperation` that applies UMAP to an `Experiment`.
"""

import inspect
import logging
from warnings import warn
from traits.api import (HasStrictTraits, Str, Dict, Any, Instance,
                        Constant, List, Enum, provides)

import numpy as np
import pandas as pd

import cytoflow.utility as util
from .i_operation import IOperation

UMAP_METRICS = ("euclidean", "manhattan", "chebyshev", "minkowski",
                "canberra", "braycurtis", "cosine", "correlation")

@provides(IOperation)
class UMAPOp(HasStrictTraits):
    """
    Use Uniform Manifold Approximation and Projection (UMAP) to reduce the
    dimensionality of the data set.

    Call `estimate` to compute the embedding.

    Calling `apply` creates new "channels" named ``{name}_1``, ``{name}_2``,
    ..., ``{name}_n``, where ``name`` is the `name` attribute and ``n`` is
    `n_components` (by default, 2.)

    The same embedding may not be appropriate for different subsets of the data set.
    If this is the case, you can use the `by` attribute to specify
    metadata by which to aggregate the data before estimating (and applying) a
    model.  The UMAP parameters such as the number of neighbors are the same
    across each subset, though.

    `estimate` fits the embedding to a random `sample` of the events; `apply`
    then uses ``UMAP.transform`` to project *all* the events (including those
    in the training sample) into that embedding. This is the same strategy
    that `tSNEOp` uses: it keeps the fitting step fast on large data sets,
    and every event is placed by the same procedure, so the result doesn't
    depend on which events happened to be in the training sample.

    UMAP is provided by the optional ``umap-learn`` package. If it isn't
    installed, `estimate` raises a `CytoflowOpError`.

    Attributes
    ----------
    name : Str
        The operation name; determines the name of the new columns.

    channels : List(Str)
        The channels to apply the embedding to.

    scale : Dict(Str : {"linear", "logicle", "log"})
        Re-scale the data in the specified channels before fitting.  If a
        channel is in `channels` but not in `scale`, the current
        package-wide default (set with `set_default_scale`) is used.

    by : List(Str)
        A list of metadata attributes to aggregate the data before estimating
        the model.  For example, if the experiment has two pieces of metadata,
        ``Time`` and ``Dox``, setting `by` to ``["Time", "Dox"]`` will
        fit the model separately to each subset of the data with a unique
        combination of ``Time`` and ``Dox``.

    n_components : Int (default = 2)
        The number of dimensions to embed the data into.

    n_neighbors : Int (default = 15)
        The size of the local neighborhood UMAP looks at when learning the
        structure of the data. Smaller values emphasize local structure;
        larger values emphasize global structure. Must be at least 2.

    min_dist : Float (default = 0.1)
        How tightly UMAP is allowed to pack points together in the
        embedding. Smaller values give tighter, "clumpier" clusters; larger
        values spread points out more evenly.

    metric : Enum (default = "euclidean")
        How to compute "distance"? One of ``euclidean``, ``manhattan``,
        ``chebyshev``, ``minkowski``, ``canberra``, ``braycurtis``,
        ``cosine`` or ``correlation``. If using many channels, try
        ``cosine``.

    random_state : Int (default = None)
        If set, seeds the random number generators (both for choosing the
        training sample and for UMAP itself) so that the embedding is
        reproducible. Setting this makes UMAP run single-threaded, so it
        will be slower.

    sample : Float (default = 0.01)
        What proportion of the data set to use for training? Defaults to 1%
        of the dataset to help with runtime.


    Notes
    -----
    Uses ``umap-learn`` by Leland McInnes, John Healy and James Melville [1]_

    References
    ----------

    .. [1] McInnes, L, Healy, J and Melville, J (2018).
           UMAP: Uniform Manifold Approximation and Projection for Dimension
           Reduction. arXiv:1802.03426
           https://arxiv.org/abs/1802.03426


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

        >>> umap = flow.UMAPOp(name = 'UMAP',
        ...                    channels = ['V2-A', 'V2-H', 'Y2-A', 'Y2-H'],
        ...                    scale = {'V2-A' : 'log',
        ...                             'V2-H' : 'log',
        ...                             'Y2-A' : 'log',
        ...                             'Y2-H' : 'log'},
        ...                    by = ["Dox"],
        ...                    random_state = 0)

    Estimate the embedding

    .. plot::
        :context: close-figs

        >>> umap.estimate(ex)

    Apply the operation

    .. plot::
        :context: close-figs

        >>> ex2 = umap.apply(ex)

    Plot a scatterplot of the embedding.  Compare to a scatterplot of
    the underlying channels.

    .. plot::
        :context: close-figs

        >>> flow.ScatterplotView(xchannel = "V2-A",
        ...                      xscale = "log",
        ...                      ychannel = "Y2-A",
        ...                      yscale = "log",
        ...                      subset = "Dox == 1.0").plot(ex2)

        >>> flow.ScatterplotView(xchannel = "UMAP_1",
        ...                      ychannel = "UMAP_2",
        ...                      subset = "Dox == 1.0").plot(ex2)

    .. plot::
        :context: close-figs

        >>> flow.ScatterplotView(xchannel = "V2-A",
        ...                      xscale = "log",
        ...                      ychannel = "Y2-A",
        ...                      yscale = "log",
        ...                      subset = "Dox == 10.0").plot(ex2)

        >>> flow.ScatterplotView(xchannel = "UMAP_1",
        ...                      ychannel = "UMAP_2",
        ...                      subset = "Dox == 10.0").plot(ex2)
    """

    id = Constant('cytoflow.operations.umap')
    friendly_id = Constant("Uniform Manifold Approximation and Projection")

    name = Str
    channels = List(Str)
    scale = Dict(Str, util.ScaleEnum)
    n_components = util.PositiveInt(2)
    n_neighbors = util.PositiveInt(15)
    min_dist = util.PositiveFloat(0.1, allow_zero = True)
    metric = Enum(*UMAP_METRICS)
    random_state = util.IntOrNone(None)
    sample = util.UnitFloat(0.01)
    by = List(Str)

    _umap = Dict(Any, Any, transient = True)
    _scale = Dict(Str, Instance(util.IScale), transient = True)

    def estimate(self, experiment, subset = None):
        """
        Estimate the embedding

        Parameters
        ----------
        experiment : Experiment
            The `Experiment` to use to estimate the embedding

        subset : str (default = None)
            A Python expression that specifies a subset of the data in
            ``experiment`` to use to parameterize the operation.

        """

        if experiment is None:
            raise util.CytoflowOpError('experiment',
                                       "No experiment specified")

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
                                           "in `channels`"
                                           .format(c))

        for b in self.by:
            if b not in experiment.data:
                raise util.CytoflowOpError('by',
                                           "Aggregation metadata {} not found, "
                                           "must be one of {}"
                                           .format(b, experiment.conditions))

        if self.n_neighbors < 2:
            raise util.CytoflowOpError('n_neighbors',
                                       "n_neighbors must be at least 2")

        # umap-learn is an optional dependency, and it's slow to import
        # (it compiles a bunch of numba code), so only import it when we
        # actually need it.
        try:
            from umap import UMAP
        except ImportError as e:
            raise util.CytoflowOpError(None,
                                       "UMAP requires the 'umap-learn' package, "
                                       "which isn't installed. Install it with "
                                       "'pip install umap-learn' or "
                                       "'conda install -c conda-forge umap-learn'.") from e

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
        for c in self.channels:
            if c in self.scale:
                self._scale[c] = util.scale_factory(self.scale[c], experiment, channel = c)
            else:
                self._scale[c] = util.scale_factory(util.get_default_scale(), experiment, channel = c)

        # umap-learn 0.6 added a new (different) layout pipeline, and warns
        # that it will become the default. Keep the classic (0.5.x) behavior
        # so that embeddings don't change depending on which umap-learn
        # version is installed.
        umap_kwargs = {}
        if 'compatibility_layout' in inspect.signature(UMAP).parameters:
            umap_kwargs['compatibility_layout'] = True

        umaps = {}
        for group, data_subset in groupby:
            if len(data_subset) == 0:
                raise util.CytoflowOpError('by',
                                           "Group {} had no data"
                                           .format(group))
            x = data_subset.sample(frac = self.sample,
                                   random_state = self.random_state).loc[:, self.channels[:]]
            for c in self.channels:
                x[c] = self._scale[c](x[c])

            # drop data that isn't in the scale range
            for c in self.channels:
                x = x[~(np.isnan(x[c]))]

            # UMAP needs more training points than neighbors (and its
            # spectral initialization needs more points than components)
            if len(x) <= max(self.n_neighbors, self.n_components + 1):
                raise util.CytoflowOpError('sample',
                                           "Group {} only had {} events to train "
                                           "on, which isn't enough for "
                                           "n_neighbors = {}. Increase the "
                                           "sample fraction or decrease "
                                           "n_neighbors."
                                           .format(group, len(x), self.n_neighbors))

            logging.info("Fitting UMAP to {} events in group {}"
                         .format(len(x), group))

            umaps[group] = \
                UMAP(n_components = self.n_components,
                     n_neighbors = self.n_neighbors,
                     min_dist = self.min_dist,
                     metric = self.metric,
                     random_state = self.random_state,
                     # setting random_state forces n_jobs = 1; say so
                     # explicitly so umap doesn't warn about it.
                     n_jobs = 1 if self.random_state is not None else -1,
                     # by default, umap uses exact (brute-force) neighbor
                     # search when fitting fewer than 4096 events. that
                     # makes transform()ing the rest of the data set
                     # (which is most of it!) extremely slow, so always
                     # use the approximate nearest-neighbor index instead.
                     force_approximation_algorithm = True,
                     **umap_kwargs).fit(x.values)

        # set this atomically to support GUI
        self._umap = umaps

    def apply(self, experiment):
        """
        Apply the UMAP embedding to the data.

        Returns
        -------
        Experiment
            a new Experiment with additional `Experiment.channels`
            named ``name_1`` ... ``name_n``, where ``n`` is `n_components`.

        """

        if experiment is None:
            raise util.CytoflowOpError('experiment',
                                       "No experiment specified")

        if not self._umap:
            raise util.CytoflowOpError(None,
                                       "No UMAP found.  Did you forget to call estimate()?")

        # make sure name got set!
        if not self.name:
            raise util.CytoflowOpError('name',
                                       "You have to set the operation's name "
                                       "before applying it!")

        if self.name != util.sanitize_identifier(self.name):
            raise util.CytoflowOpError('name',
                                       "Name can only contain letters, numbers and underscores."
                                       .format(self.name))

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

        if self.by:
            groupby = experiment.data.groupby(self.by, observed = False)
        else:
            # use a lambda expression to return a group that contains
            # all the events
            groupby = experiment.data.groupby(lambda _: True, observed = False)

        # use the number of components the models were actually fit with,
        # in case n_components changed after estimate() was called.
        n_components = next(iter(self._umap.values())).n_components

        # need deep = True because of the data.dropna below
        new_experiment = experiment.clone(deep = True)
        new_channels = []
        for i in range(n_components):
            cname = "{}_{}".format(self.name, i + 1)
            if cname in experiment.data:
                raise util.CytoflowOpError('name',
                                           "Channel {} is already in the experiment"
                                           .format(cname))

            new_experiment.add_channel(cname, pd.Series(index = experiment.data.index))
            new_channels.append(cname)

        for group, data_subset in groupby:
            if len(data_subset) == 0:
                warn("Group {} had no data".format(group),
                     util.CytoflowOpWarning)
                continue

            if group not in self._umap:
                raise util.CytoflowOpError('by',
                                           "No UMAP model was estimated for "
                                           "group {}".format(group))

            x = data_subset.loc[:, self.channels[:]]
            for c in self.channels:
                x[c] = self._scale[c](x[c])

            # which values are missing (say, negative with a log transform?)
            # we will set them to 0 because we aren't estimating the embedding;
            # and set them back to NaN after.

            x_na = pd.Series([False] * len(x))
            for c in self.channels:
                x_na[np.isnan(x[c]).values] = True
            x_na = x_na.values
            x[x_na] = 0

            group_idx = data_subset.index

            logging.info("Applying UMAP to {} events in group {}"
                         .format(len(x), group))

            # project ALL the events (not just the training sample) into
            # the embedding.
            x_tf = self._umap[group].transform(x.values)
            x_tf[x_na] = np.nan

            for ci, c in enumerate(new_channels):
                new_experiment.data.loc[group_idx, c] = x_tf[:, ci]

        new_experiment.data.dropna(inplace = True)
        new_experiment.data.reset_index(drop = True, inplace = True)

        new_experiment.history.append(self.clone_traits(transient = lambda _: True))
        return new_experiment
