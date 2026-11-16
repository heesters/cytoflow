.. _user_gate_interchange:


HOWTO: Exchange gates with FlowJo and other software
====================================================

``Cytoflow`` can share gates with other flow cytometry software in two ways:

* It can **export** and **import** `Gating-ML 2.0 <https://isac-net.org/page/Data-Standards>`_
  files.  Gating-ML is the ISAC standard for describing gates, and it is
  supported (to varying degrees) by FlowJo, Cytobank, FlowKit, the 
  R/Bioconductor packages (``CytoML``, ``flowCore``, ``openCyto``) and others.
  The files that ``Cytoflow`` writes have been checked with FlowKit.

* It can **import** the gates from a FlowJo 10 workspace (a ``.wsp`` file)
  directly.

A few things work differently in ``Cytoflow`` than they do in those programs,
so please read the section on :ref:`what to expect <user_gate_interchange_expect>`,
and **always check a few event counts** against the other program.


Exporting gates
---------------

#. Set up your workflow as usual.  Make sure that all of the gates you want
   to export are "valid" (they have a green check-mark.)

#. Choose **File → Export gates (Gating-ML)...** and choose a file name.

#. ``Cytoflow`` tells you how many gates it wrote, and warns you about any
   that it couldn't.

Every gate in the workflow is exported:

====================================  ==================================
``Cytoflow`` gate                     Gating-ML gate
====================================  ==================================
Threshold, Range, 2D Range            ``RectangleGate``
Polygon, Freehand                     ``PolygonGate``
Ellipse                               ``EllipsoidGate``
Quad                                  ``QuadrantGate`` (the quadrants are
                                      named ``<name>_1`` to ``<name>_4``,
                                      like the gate's categories)
Boolean                               ``BooleanGate``
====================================  ==================================

The gate's name in ``Cytoflow`` is its ID in the file.  Other kinds of
operations that make "gates" -- density gates, Gaussian mixture models,
clustering -- can't be represented in Gating-ML and are skipped.

From a script or notebook, use `cytoflow.export_gatingml`::

    flow.export_gatingml([threshold_op, polygon_op, boolean_op], ex, "gates.xml")


Importing gates
---------------

#. Import your data, and do any compensation or other calibration.  (The
   gates are matched to the channels of the *last* operation in the workflow.)

#. Choose **File → Import gates...** and choose a Gating-ML file (``.xml``)
   or a FlowJo workspace (``.wsp``).

#. *FlowJo workspaces only:* FlowJo keeps a separate set of gates for every
   sample (and a template for each group), but ``Cytoflow`` applies one set of
   gates to all of your data.  If the workspace has more than one sample,
   choose the sample or group whose gates you want.  (A group's gates are 
   interpreted with the axis settings of the first sample in the group.)

#. The new gates are added to the end of the workflow.  They're ordinary
   ``Cytoflow`` operations: you can look at them, move their vertices, rename
   them or delete them.

#. **Read the message that pops up.**  It lists any gates that couldn't be
   imported, and anything that was approximated.

From a script or notebook, use `cytoflow.import_gatingml` or 
`cytoflow.import_flowjo_workspace`; both return a list of operations::

    for op in flow.import_flowjo_workspace("experiment.wsp", ex, sample = "Tube_001.fcs"):
        ex = op.apply(ex)


.. _user_gate_interchange_expect:

What to expect
--------------

Nested gates
    In FlowJo (and Gating-ML) gates are *nested*: the ``Singlets`` gate
    is drawn "inside" the ``Cells`` gate, and the ``Singlets`` population only
    has events that are in both.  In ``Cytoflow``, every gate applies to *all* 
    of the events, and you combine them afterwards.

    So, a nested gate becomes **two** operations.  ``Singlets`` is the gate
    itself, applied to all the events -- you can see and edit its outline.
    ``Singlets_gated`` is a **Boolean Gate** with the expression 
    ``Singlets and Cells``; it has the same events as FlowJo's ``Singlets``
    population, and it's the one you usually want to use for subsets and 
    statistics.  If ``Cells`` were nested too, the expression would be 
    ``Singlets and Cells_gated``.

    When you *export*, there's no nesting to write: export the Boolean Gates
    along with the gates they combine.

Names
    ``Cytoflow`` names can only have letters, numbers and underscores.
    ``CD4+`` becomes ``CD4pos``, ``CD8-`` becomes ``CD8neg``, other characters
    become ``_``, and a name that's already used in your experiment gets a
    number added to it.

Channel names
    Gates are matched to channels by the channel's name *in the FCS file*,
    even if you renamed the channel when you imported your data.  Gates on 
    channels that aren't in your experiment are skipped.

Compensation
    **Compensation is not exchanged.**  Imported gates are applied to your
    channels as they are, so if the gates were drawn on compensated data
    (FlowJo's ``Comp-FITC-A``), compensate your data in ``Cytoflow`` first -- 
    with the same matrix!  ``Cytoflow`` warns you when this is the case.  
    Similarly, exported gates describe the data as it is in ``Cytoflow`` 
    *at that point in the workflow*.

Axis scales
    Polygon and ellipse gates depend on the scale of the plot they were drawn 
    on: a straight line on a logicle plot isn't straight on a linear one.  
    ``Cytoflow`` writes its scales as the equivalent Gating-ML transforms
    (linear → ``flin``, log → ``flog``, logicle → ``logicle``, 
    arcsinh → ``fasinh``) with the parameters it uses for that channel, and 
    converts the gates' coordinates, so the exported gate has exactly the same 
    outline.

    Going the other way, the file's axis settings often *don't* match
    ``Cytoflow``'s (for instance, a logicle scale with a different width.)  
    Rectangle gates aren't affected.  A polygon is given extra vertices
    so that it follows the original outline, and an ellipse is converted to a 
    128-sided polygon.  You get the right events, but the gate has a lot of 
    vertices.  To avoid this, use the **Axis Transforms** operation to give 
    the channel the same scale settings as the other program *before* you
    import the gates.

FlowJo's biexponential scale
    FlowJo's default "biex" scale is closely related to the logicle scale, but
    FlowJo doesn't publish how it is computed.  ``Cytoflow`` uses a logicle 
    scale with the equivalent settings, which differs from FlowJo's by up to
    about half a percent of the width of the axis.  That can move a few events
    across the edge of a polygon or ellipse gate.  (Rectangle and quadrant
    gates are exact.)

Edges
    Programs differ about events *exactly* on the edge of a gate, and
    ``Cytoflow``'s logicle scale is computed with a (fine) lookup table.  
    Expect an occasional difference of an event or two.


What isn't supported
--------------------

These are skipped, with a warning that says which gates and why:

* Gates with more than two dimensions.

* Gates on ratios of channels (Gating-ML's ``fratio``) and on FlowJo's 
  derived parameters.

* FlowJo's "curly" quadrant gates.

* FlowJo polygon and ellipse gates on axes other than linear, log, 
  biexponential, logicle and arcsinh (for example, hyperlog or Miltenyi 
  scales).  Rectangle gates on those axes are fine.

* FlowJo ellipse gates on a linear axis that doesn't start at zero.

* Anything nested inside, or combined with, a gate that was skipped.

Also:

* Compensation matrices are neither imported nor exported.

* Per-sample differences in FlowJo's gates aren't imported -- you get one
  sample's (or one group's) gates.  ``Cytoflow``'s own "magnetic" gates are
  exported at their un-adjusted position.

* Exporting to a FlowJo workspace isn't supported; FlowJo can import
  Gating-ML.

.. note:: 

   The FlowJo workspace format is not publicly documented.  The importer is
   based on other open-source projects' descriptions of the format, and it has
   been tested against a hand-written workspace file and against FlowKit's 
   importer -- but not (yet) against a wide range of real workspaces.  If a 
   gate comes out wrong, please file a bug report and attach the workspace.
