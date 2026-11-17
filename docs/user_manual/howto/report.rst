.. _user_report:


HOWTO: Export several plots as a PDF or PowerPoint report
=========================================================

**File > Save Plot...** saves the plot that is on the screen.  To put 
*several* plots on a page -- for a lab meeting, a lab notebook or a figure
draft -- use **File > Export report (PDF / PowerPoint)...**

#. Set up the plots you want as usual.  Every view that you have opened on
   an operation is remembered with its settings (channels, scales, gates,
   titles and other plot parameters), and the report draws each one 
   exactly as *Cytoflow* does.

#. Choose **File > Export report (PDF / PowerPoint)...**  The dialog lists 
   every plot in the workflow, named after its operation and the plot's 
   title (or type.)  Untick the ones you don't want.

#. Choose a title (shown on every page), and how many rows and columns of 
   plots to put on each page.

#. Click **OK** and choose a file name ending in ``.pdf`` or ``.pptx``.

The report is made in the background -- the busy indicator shows in the 
corner of the plot -- and a message tells you when it's done.  If a plot 
can't be made (for example, because its operation has an error), the others
are still exported, and the message says which ones are missing.

PDF or PowerPoint?
------------------

* A **PDF with 1 row and 1 column** has one plot per page as a *vector*
  graphic: it stays sharp at any size, and the text can be edited in 
  Illustrator, Affinity Designer or Inkscape.  Use this for figures.
  
* A **PDF with several plots per page** places the plots on the page as 
  images, at the resolution you choose (300 dpi by default.)
  
* A **PowerPoint** file has the same layout on 16:9 slides.  Each plot is a
  picture that you can move, resize or delete, and each caption is an 
  ordinary text box.  The plots are images, not editable charts.

From a script or notebook
-------------------------

The same function is available without the GUI:

.. code-block:: python

    import cytoflow as flow
    
    items = [flow.ReportItem(flow.HistogramView(channel = "FITC-A", scale = "logicle"),
                             ex, caption = "FITC"),
             flow.ReportItem(gate_op.default_view(), ex, caption = "Live cells")]
    
    flow.export_report(items, "report.pptx", title = "Experiment 12", layout = (1, 2))

PowerPoint export needs the ``python-pptx`` package, which is included in the
point-and-click application.
