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
cytoflow.operations.compensation
--------------------------------

The `compensation` module applies the compensation (spillover) matrix that
the cytometer stored in the FCS file -- or one the user supplies, for example
from a CSV file exported from FlowJo.  It contains:

`CompensationOp` -- applies a spillover matrix to an `Experiment`.

`CompensationDiagnostic` -- a heat-map of the spillover matrix.

`parse_spillover` -- parse the value of an FCS ``$SPILLOVER`` / ``SPILL`` /
``$SPILL`` keyword into a `pandas.DataFrame`.

`fcs_spillover` -- find and parse the spillover keyword in a dictionary of
FCS metadata.

`read_spillover_csv` -- read a spillover matrix from a CSV file.

The math
~~~~~~~~

The spillover matrix :math:`S` has one row and one column per fluorescence
channel; :math:`S_{ij}` is the fraction of the signal from the fluorophore
measured in channel :math:`i` that is detected in channel :math:`j`, so
the diagonal is (usually) 1.  This is the convention of the FCS 3.1
standard's ``$SPILLOVER`` keyword (the values are listed row by row,
row :math:`i` being the spillover *from* parameter :math:`i`), and the one
used by FlowJo, FlowKit / FlowUtils and flowCore.

If :math:`x` is an event's vector of (true) fluorophore signals, the
cytometer measures the row vector :math:`r = x S`; so the compensated
data are

.. math::

    x = r S^{-1}

which is what FlowUtils (``np.linalg.solve(S.T, r.T).T``) and flowCore
(``r %*% solve(S)``) compute.
"""

import os, warnings

from natsort import natsorted

from traits.api import (HasStrictTraits, Str, File, Dict, Float, Constant,
                        Instance, provides)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import cytoflow.views
import cytoflow.utility as util

from .i_operation import IOperation

# the FCS keywords that may hold the spillover matrix, in order of preference.
SPILLOVER_KEYWORDS = ['$SPILLOVER', 'SPILL', '$SPILL', 'SPILLOVER']

# prefixes FlowJo et al. add to compensated channel names
_COMP_PREFIXES = ['FJComp-', 'Comp-', 'Comp ', '<', '[']


class SpilloverMatrix(Dict):
    """
    A trait holding a spillover matrix as a dict-of-dicts:
    ``spillover[from_channel][to_channel] = value``.  Also accepts a
    `pandas.DataFrame` (rows are the "from" channels, columns are
    the "to" channels), which it converts to a dict-of-dicts.
    """

    def __init__(self, **metadata):
        super().__init__(Str, Dict(Str, Float), **metadata)

    def validate(self, obj, name, value):
        if isinstance(value, pd.DataFrame):
            value = spillover_to_dict(value)
        return super().validate(obj, name, value)


def spillover_to_dict(df):
    """Convert a spillover `pandas.DataFrame` to a dict-of-dicts"""
    return {str(f) : {str(t) : float(df.loc[f, t]) for t in df.columns}
            for f in df.index}


def spillover_to_dataframe(spillover):
    """
    Convert a dict-of-dicts spillover matrix to a square `pandas.DataFrame`.
    Missing entries are 1.0 on the diagonal and 0.0 elsewhere.
    """
    channels = list(spillover.keys())
    for f in spillover:
        for t in spillover[f]:
            if t not in channels:
                channels.append(t)

    return pd.DataFrame([[spillover.get(f, {}).get(t, 1.0 if f == t else 0.0)
                          for t in channels]
                         for f in channels],
                        index = channels,
                        columns = channels,
                        dtype = float)


def parse_spillover(text):
    """
    Parse the value of an FCS spillover keyword
    (``n,ch_1,...,ch_n,s_11,s_12,...,s_nn``).

    Returns
    -------
    pandas.DataFrame
        The spillover matrix; the index is the channel the signal spills
        *from* and the columns are the channels it spills *to*.
    """
    items = [x.strip() for x in str(text).split(',')]
    try:
        n = int(float(items[0]))
    except ValueError as e:
        raise util.CytoflowError("Can't parse spillover matrix '{}'"
                                 .format(text)) from e

    if n < 1 or len(items) < 1 + n + n * n:
        raise util.CytoflowError("Spillover matrix '{}' is too short "
                                 "(expected {} channels and {} values)"
                                 .format(text, n, n * n))

    channels = items[1 : n + 1]
    try:
        values = np.array([float(x) for x in items[n + 1 : n + 1 + n * n]])
    except ValueError as e:
        raise util.CytoflowError("Can't parse the values of spillover matrix '{}'"
                                 .format(text)) from e

    # row-major: row i is the spillover from channel i
    return pd.DataFrame(values.reshape((n, n)),
                        index = channels,
                        columns = channels)


def fcs_spillover(fcs_metadata):
    """
    Find and parse the spillover matrix in a dictionary of FCS metadata
    (keywords ``$SPILLOVER``, ``SPILL``, ``$SPILL``).  Channel names in the
    matrix are the FCS ``$PnN`` names; some instruments (eg the Accuri C6)
    write parameter numbers instead, which are converted to ``$PnN`` names.

    Returns
    -------
    pandas.DataFrame or None
        The spillover matrix, or ``None`` if there isn't one.
    """
    keys = {k.upper() : k for k in fcs_metadata.keys()}
    for kw in SPILLOVER_KEYWORDS:
        if kw in keys and str(fcs_metadata[keys[kw]]).strip():
            df = parse_spillover(fcs_metadata[keys[kw]])
            break
    else:
        return None

    # are the channels parameter numbers instead of names?
    pnn = _pnn_names(fcs_metadata)
    if pnn and not set(df.index) & set(pnn.values()) \
           and all(c.isdigit() and int(c) in pnn for c in df.index):
        names = [pnn[int(c)] for c in df.index]
        df.index = names
        df.columns = names

    return df


def read_spillover_csv(filename):
    """
    Read a spillover matrix from a CSV file, such as one exported from
    FlowJo's compensation editor.

    The first row must be the channel names.  Each following row is the
    spillover *from* one channel, in the same order as the columns; it may
    optionally start with that channel's name.  If the diagonal is ~100,
    the values are taken to be percentages and are divided by 100.

    Returns
    -------
    pandas.DataFrame
        The spillover matrix (rows are "from", columns are "to").
    """
    if not os.path.isfile(filename):
        raise util.CytoflowError("Can't find file {}".format(filename))

    df = pd.read_csv(filename, header = 0, skipinitialspace = True,
                     sep = None, engine = 'python')

    # drop empty columns / rows (eg trailing commas)
    df = df.dropna(axis = 'columns', how = 'all').dropna(axis = 'index', how = 'all')

    # is the first column the row names?
    first = df.columns[0]
    if not pd.api.types.is_numeric_dtype(df[first]):
        df = df.set_index(first)
    elif len(df.columns) == len(df) + 1:
        df = df.set_index(first)
    else:
        df.index = df.columns[0:len(df)]

    df.index = [str(x).strip() for x in df.index]
    df.columns = [str(x).strip() for x in df.columns]

    if len(df.index) != len(df.columns):
        raise util.CytoflowError("The spillover matrix in {} isn't square "
                                 "({} rows, {} columns)"
                                 .format(filename, len(df.index), len(df.columns)))

    if set(df.index) != set(df.columns):
        if set(df.index) and all(i.startswith("Unnamed") for i in df.index):
            df.index = df.columns
        else:
            raise util.CytoflowError("The row names and column names in {} "
                                     "don't match".format(filename))

    try:
        df = df.astype(float)
    except ValueError as e:
        raise util.CytoflowError("Non-numeric values in the spillover matrix "
                                 "in {}".format(filename)) from e

    diag = np.diag(df.loc[df.columns, df.columns].values)
    if np.all(np.abs(diag - 100.0) < 1e-6):
        df = df / 100.0

    return df


def _pnn_names(fcs_metadata):
    """Return a dict of parameter number --> $PnN"""
    try:
        par = int(fcs_metadata['$PAR'])
    except (KeyError, ValueError):
        return {}

    return {i : str(fcs_metadata["$P{}N".format(i)]).strip()
            for i in range(1, par + 1) if "$P{}N".format(i) in fcs_metadata}


def _channel_map(experiment, fcs_metadata):
    """
    Return a dict mapping FCS names ($PnN and $PnS) to the experiment's channel
    names.  The experiment's channels may have been renamed when they were
    imported; the original name is in ``experiment.metadata[channel]['fcs_name']``,
    and is either the $PnN or the $PnS name depending on
    ``experiment.metadata['name_metadata']``.
    """
    fcs_names = {experiment.metadata[c].get('fcs_name', c) : c
                 for c in experiment.channels}

    name_metadata = experiment.metadata.get('name_metadata', '$PnN')
    ret = {}

    if fcs_metadata:
        try:
            par = int(fcs_metadata['$PAR'])
        except (KeyError, ValueError):
            par = 0

        for i in range(1, par + 1):
            pnn = str(fcs_metadata.get("$P{}N".format(i), "")).strip()
            pns = str(fcs_metadata.get("$P{}S".format(i), "")).strip()

            preferred = pns if name_metadata == '$PnS' else pnn
            other = pnn if name_metadata == '$PnS' else pns

            if preferred and preferred in fcs_names:
                c = fcs_names[preferred]
            elif other and other in fcs_names:
                c = fcs_names[other]
            else:
                continue

            if pnn: ret.setdefault(pnn, c)
            if pns: ret.setdefault(pns, c)

    # fall back to the fcs names and channel names themselves
    for fcs_name, c in fcs_names.items():
        ret.setdefault(fcs_name, c)

    # channel names always map to themselves
    for c in experiment.channels:
        ret[c] = c

    return ret


def _lookup(name, channel_map):
    """Find a matrix channel name in the channel map, stripping prefixes"""
    if name in channel_map:
        return channel_map[name]

    stripped = name
    for prefix in _COMP_PREFIXES:
        if stripped.startswith(prefix):
            stripped = stripped[len(prefix):]

    stripped = stripped.rstrip(']>').strip()
    return channel_map.get(stripped, None)


@provides(IOperation)
class CompensationOp(HasStrictTraits):
    """
    Compensate fluorescence channels using a spillover matrix -- usually the
    one the cytometer stored in the FCS file.

    Most cytometers' acquisition software write the compensation matrix that
    was set up during acquisition into the FCS file (in the ``$SPILLOVER``,
    ``SPILL`` or ``$SPILL`` keyword).  `estimate` reads that matrix from the
    experiment's first tube (or from a CSV file, if `spillover_file` is set),
    and `apply` compensates the data.  You can also set (or edit) `spillover`
    yourself.

    The compensated data are :math:`x = r S^{-1}`, where :math:`r` is the
    row vector of raw channel values and :math:`S` is the spillover matrix
    (:math:`S_{ij}` is the fraction of channel :math:`i`'s fluorophore that is
    detected in channel :math:`j`).  This is the FCS 3.1 convention, and the
    one used by FlowJo, FlowKit and flowCore.

    .. note::
        Don't use this *and* `BleedthroughLinearOp` -- that would compensate
        the data twice!  The data should not have been compensated before,
        and compensation should come before any other operation that changes
        the fluorescence channels (eg autofluorescence correction).

    Attributes
    ----------
    spillover : Dict(Str, Dict(Str, Float))
        The spillover matrix, as a dictionary of dictionaries:
        ``spillover[from_channel][to_channel]`` is the fraction of
        ``from_channel``'s signal that appears in ``to_channel``.  You can
        also assign a `pandas.DataFrame` whose index is the "from" channels
        and whose columns are the "to" channels.  Missing entries are 1 on the
        diagonal and 0 elsewhere.  Channel names are the experiment's channel
        names.  Set by `estimate`, or set it yourself.

    spillover_file : File
        (Optional) If set, `estimate` reads the spillover matrix from this
        CSV file instead of the FCS metadata.  See `read_spillover_csv` for
        the format.

    Examples
    --------

    Import some data from a BD FACSAria II, which stored its spillover matrix
    in the FCS file:

    .. plot::
        :context: close-figs

        >>> import cytoflow as flow
        >>> import_op = flow.ImportOp()
        >>> import_op.tubes = [flow.Tube(file = "instruments/BD - FACS Aria II.fcs")]
        >>> ex = import_op.apply()

    Read the spillover matrix from the FCS file

    .. plot::
        :context: close-figs

        >>> comp_op = flow.CompensationOp()
        >>> comp_op.estimate(ex)

    Plot the spillover matrix

    .. plot::
        :context: close-figs

        >>> comp_op.default_view().plot(ex)

    Apply the compensation

    .. plot::
        :context: close-figs

        >>> ex2 = comp_op.apply(ex)
    """

    # traits
    id = Constant('cytoflow.operations.compensation')
    friendly_id = Constant("Compensation")

    name = Constant("Compensation")

    spillover = SpilloverMatrix()
    spillover_file = File(exists = False)

    def estimate(self, experiment, subset = None):
        """
        Read the spillover matrix from `spillover_file` (if it's set), or
        from the FCS metadata of the experiment's first tube, and map its
        channel names to the experiment's channel names.

        Matrix channels that aren't in the experiment are dropped, with a
        warning.

        Parameters
        ----------
        experiment : Experiment
            The `Experiment` whose metadata to read.

        subset : str
            Ignored -- this operation doesn't need any data to estimate its
            parameters.
        """

        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")

        fcs_metadata = experiment.metadata.get('fcs_metadata', {})
        tube_meta = next(iter(fcs_metadata.values())) if fcs_metadata else {}

        if self.spillover_file:
            try:
                df = read_spillover_csv(self.spillover_file)
            except util.CytoflowError as e:
                raise util.CytoflowOpError('spillover_file', str(e)) from e
            source = self.spillover_file
        else:
            if not fcs_metadata:
                raise util.CytoflowOpError('experiment',
                                           "The experiment has no FCS metadata. "
                                           "Was it created with ImportOp?")
            filename = next(iter(fcs_metadata.keys()))
            df = fcs_spillover(tube_meta)
            if df is None:
                raise util.CytoflowOpError('experiment',
                                           "Can't find a spillover matrix "
                                           "({}) in the metadata of {}. You "
                                           "can load one from a CSV file "
                                           "with 'spillover_file' instead."
                                           .format(", ".join(SPILLOVER_KEYWORDS),
                                                   filename))
            source = filename

            # warn if the tubes have different matrices
            for other_file, other_meta in list(fcs_metadata.items())[1:]:
                other_df = fcs_spillover(other_meta)
                if other_df is None or not other_df.equals(df):
                    warnings.warn("Tube {} has a different spillover matrix "
                                  "than {}; using the one from {}"
                                  .format(other_file, filename, filename),
                                  util.CytoflowOpWarning)
                    break

        channel_map = _channel_map(experiment, tube_meta)

        names = {}
        missing = []
        for c in df.index:
            new_name = _lookup(c, channel_map)
            if new_name is None:
                missing.append(c)
            else:
                names[c] = new_name

        if len(set(names.values())) != len(names):
            raise util.CytoflowOpError('spillover',
                                       "More than one channel in the spillover "
                                       "matrix maps to the same experiment channel")

        if not names:
            raise util.CytoflowOpError('spillover',
                                       "None of the channels in the spillover "
                                       "matrix from {} ({}) are in the experiment"
                                       .format(source, ", ".join(df.index)))

        if missing:
            warnings.warn("Channels {} are in the spillover matrix from {} "
                          "but not in the experiment; ignoring them"
                          .format(", ".join(missing), source),
                          util.CytoflowOpWarning)

        df = df.loc[list(names.keys()), list(names.keys())]
        df = df.rename(index = names, columns = names)

        mat = df.values
        if np.allclose(mat, np.eye(len(mat))):
            warnings.warn("The spillover matrix from {} is the identity "
                          "matrix -- compensation won't change the data"
                          .format(source),
                          util.CytoflowOpWarning)

        # set atomically to support the GUI
        self.spillover = df


    def apply(self, experiment):
        """
        Compensate the data in an `Experiment`.

        Parameters
        ----------
        experiment : `Experiment`
            The experiment to which this operation is applied

        Returns
        -------
        Experiment
            A new `Experiment` with the channels in `spillover` compensated.
            The compensated channels have the following metadata added:

            - **compensated** : Bool
              ``True``

            - **spillover** : Dict(Str : Float)
              The spillover from the other channels into this one.
        """

        if experiment is None:
            raise util.CytoflowOpError('experiment', "No experiment specified")

        if not self.spillover:
            raise util.CytoflowOpError('spillover',
                                       "Spillover matrix isn't set. "
                                       "Did you forget to call estimate()?")

        df = spillover_to_dataframe(self.spillover)
        channels = list(df.index)

        for c in channels:
            if c not in experiment.channels:
                raise util.CytoflowOpError('spillover',
                                           "Channel {} is in the spillover "
                                           "matrix but not in the experiment"
                                           .format(c))

            if experiment.metadata[c].get('compensated', False) or \
               'bleedthrough_channels' in experiment.metadata[c]:
                warnings.warn("Channel {} has already been compensated!"
                              .format(c),
                              util.CytoflowOpWarning)

        s = df.values
        if not np.all(np.isfinite(s)):
            raise util.CytoflowOpError('spillover',
                                       "The spillover matrix has missing or "
                                       "non-finite values")

        try:
            s_inv = np.linalg.inv(s)
        except np.linalg.LinAlgError as e:
            raise util.CytoflowOpError('spillover',
                                       "The spillover matrix is singular and "
                                       "can't be inverted") from e

        new_experiment = experiment.clone(deep = True)

        compensated = experiment.data[channels].values @ s_inv

        for i, c in enumerate(channels):
            new_experiment.data[c] = compensated[:, i]
            new_experiment.metadata[c]['compensated'] = True
            new_experiment.metadata[c]['spillover'] = \
                {f : float(df.loc[f, c]) for f in channels if f != c}

        new_experiment.history.append(self.clone_traits(transient = lambda _: True))
        return new_experiment


    def default_view(self, **kwargs):
        """
        Returns a diagnostic view: a heat map of the spillover matrix.

        Returns
        -------
        `IView`
            Call `CompensationDiagnostic.plot` to see the matrix.
        """
        v = CompensationDiagnostic(op = self)
        v.trait_set(**kwargs)
        return v


@provides(cytoflow.views.IView)
class CompensationDiagnostic(HasStrictTraits):
    """
    Plots the spillover matrix as a heat map.  Rows are the channels the
    signal spills *from*; columns are the channels it spills *into*.  The
    off-diagonal values are printed in each cell (as percentages).

    Attributes
    ----------
    op : Instance(CompensationOp)
        The operation whose spillover matrix to plot.  If you made the
        instance with `CompensationOp.default_view`, this is set for you.
    """

    id = Constant("cytoflow.view.compensationdiagnostic")
    friendly_id = Constant("Compensation Diagnostic")

    op = Instance(IOperation)

    # unused, but part of the interface the GUI expects
    subset = Str

    def plot(self, experiment = None, **kwargs):
        """
        Plot the spillover matrix.  ``experiment`` is optional; if it's given,
        the matrix's channels are sorted the same way as the experiment's
        channels.
        """

        if self.op is None or not self.op.spillover:
            raise util.CytoflowViewError('op',
                                         "No spillover matrix specified. "
                                         "Did you forget to call estimate()?")

        df = spillover_to_dataframe(self.op.spillover)
        channels = natsorted(df.index)
        df = df.loc[channels, channels]

        n = len(channels)
        values = df.values * 100.0
        offdiag = values.copy()
        np.fill_diagonal(offdiag, np.nan)
        vmax = np.nanmax(np.abs(offdiag)) if n > 1 else 1.0
        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1.0

        ax = plt.gca()

        cmap = plt.get_cmap(kwargs.pop('cmap', 'Reds')).copy()
        cmap.set_bad('lightgrey')
        im = ax.imshow(np.abs(offdiag), cmap = cmap, vmin = 0, vmax = vmax,
                       interpolation = 'nearest')

        if n <= 30:
            for i in range(n):
                for j in range(n):
                    v = values[i, j]
                    if i == j:
                        text = "{:g}".format(v / 100.0) if not np.isclose(v, 100.0) else ""
                    elif v == 0:
                        text = "0"
                    else:
                        text = "{:.3g}".format(v)
                    ax.text(j, i, text, ha = 'center', va = 'center',
                            fontsize = 'small' if n <= 10 else 'xx-small',
                            color = 'white' if abs(offdiag[i, j]) > vmax * 0.6 else 'black')

        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(channels, rotation = 45, ha = 'right')
        ax.set_yticklabels(channels)
        ax.set_xlabel("Spillover into")
        ax.set_ylabel("Spillover from")
        ax.set_title("Spillover (%)")
        ax.grid(False)
        plt.colorbar(im, ax = ax, label = "Spillover (%)", shrink = 0.8)
