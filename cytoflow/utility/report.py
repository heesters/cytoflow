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
cytoflow.utility.report
-----------------------

Export several plots as one report: a PDF, or a PowerPoint presentation.

`render_view` -- render one view to an image (or a vector graphic.)

`export_report` -- lay out several plots on pages (or slides) and save them.

`ReportItem` -- one plot of a report.
"""

import io, os, warnings

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure

from .cytoflow_errors import CytoflowError, CytoflowWarning

# width x height, in inches
PAGE_SIZES = {"A4 landscape" : (11.69, 8.27),
              "A4 portrait" : (8.27, 11.69),
              "Letter landscape" : (11.0, 8.5),
              "Letter portrait" : (8.5, 11.0),
              "Slide 16:9" : (13.333, 7.5),
              "Slide 4:3" : (10.0, 7.5)}


class ReportItem(object):
    """
    One plot of a report.
    
    Parameters
    ----------
    view : IView
        The view to plot.
        
    experiment : Experiment
        The experiment to plot it with.
        
    caption : str
        The text under the plot.
        
    plot_kwargs : dict
        Passed to ``view.plot()``.
    """
    
    def __init__(self, view, experiment, caption = "", plot_kwargs = None):
        self.view = view
        self.experiment = experiment
        self.caption = caption
        self.plot_kwargs = dict(plot_kwargs or {})
        
    def plot(self):
        """Plot the view (with pyplot) and return the new figure."""
        before = set(plt.get_fignums())
        self.view.plot(self.experiment, **self.plot_kwargs)
        fig = plt.gcf()
        # close any other figures that the view made along the way
        for num in set(plt.get_fignums()) - before - {fig.number}:
            plt.close(num)
        return fig


def render_view(item, fmt = "png", dpi = 300):
    """
    Render one `ReportItem` and return the image file's contents.
    
    Parameters
    ----------
    item : ReportItem
    
    fmt : {"png", "pdf", "svg"}
        The image format.
        
    dpi : int
        The resolution, for ``png``.
        
    Returns
    -------
    bytes
    """
    fig = item.plot()
    try:
        buf = io.BytesIO()
        with matplotlib.rc_context({'svg.fonttype' : 'none', 'pdf.fonttype' : 42}):
            fig.savefig(buf, format = fmt, dpi = dpi, bbox_inches = 'tight')
        return buf.getvalue()
    finally:
        plt.close(fig)


def _render_all(items, dpi):
    """
    Render the items to PNGs.  Returns a list of (png bytes or None, 
    caption) and a list of the problems.
    """
    images = []
    problems = []
    for i, item in enumerate(items):
        try:
            images.append((render_view(item, "png", dpi), item.caption))
        except Exception as e:
            what = item.caption or "plot {}".format(i + 1)
            msg = e.args[-1] if isinstance(e, CytoflowError) and e.args else str(e)
            problems.append("{}: {}".format(what, msg))
            images.append((None, "{}\n(couldn't be plotted: {})".format(what, msg)))
    return images, problems


def _pages(images, per_page):
    return [images[i : i + per_page] for i in range(0, len(images), per_page)]


def _cells(n_rows, n_cols, x0, y0, width, height, gap):
    """The (left, top, width, height) of each cell of a grid, row by row"""
    w = (width - gap * (n_cols - 1)) / n_cols
    h = (height - gap * (n_rows - 1)) / n_rows
    return [(x0 + c * (w + gap), y0 + r * (h + gap), w, h)
            for r in range(n_rows) for c in range(n_cols)]


def _fit(img_w, img_h, box_w, box_h):
    """Scale an image to fit in a box, keeping its shape.  Returns (w, h)."""
    scale = min(box_w / img_w, box_h / img_h)
    return img_w * scale, img_h * scale


def _png_size(data):
    # the IHDR chunk: width and height are big-endian at bytes 16..24
    return int.from_bytes(data[16:20], 'big'), int.from_bytes(data[20:24], 'big')


def export_report(items, filename, title = "", layout = (2, 2), 
                  page_size = None, dpi = 300):
    """
    Save several plots as a PDF or a PowerPoint presentation.
    
    Parameters
    ----------
    items : list of ReportItem
        The plots, in order.
        
    filename : str
        The file to write.  The extension chooses the format: ``.pdf`` or 
        ``.pptx``.
        
    title : str
        Shown at the top of every page / slide.
        
    layout : (int, int)
        The number of rows and columns of plots on each page.  
        
        For a PDF, ``(1, 1)`` gives one plot per page as a **vector** 
        graphic: sharp at any zoom, with text you can edit in Illustrator,
        Affinity Designer or Inkscape.  With more than one plot per page, 
        the plots are placed on the page as images (at ``dpi``.)
        
        PowerPoint slides always hold the plots as images, each of which 
        you can move and resize, with its caption in a text box.
        
    page_size : (float, float) or str
        The page's width and height in inches, or one of the names in 
        `PAGE_SIZES`.  Default: A4 landscape for a PDF, a 16:9 slide for 
        PowerPoint.
        
    dpi : int
        The resolution of the images.
        
    Returns
    -------
    list of str
        A message for every plot that couldn't be made (it is replaced by 
        the message in the report.)  Empty if all went well.
    """
    items = list(items)
    if not items:
        raise CytoflowError("There are no plots to export")
    
    n_rows, n_cols = int(layout[0]), int(layout[1])
    if n_rows < 1 or n_cols < 1:
        raise CytoflowError("The layout must have at least one row and one column")
    
    ext = os.path.splitext(filename)[1].lower()
    if ext not in (".pdf", ".pptx"):
        raise CytoflowError("Can't export a report to a '{}' file; use .pdf or .pptx"
                            .format(ext))
        
    if page_size is None:
        page_size = "A4 landscape" if ext == ".pdf" else "Slide 16:9"
    if isinstance(page_size, str):
        if page_size not in PAGE_SIZES:
            raise CytoflowError("Unknown page size '{}'; choose from {}"
                                .format(page_size, ", ".join(PAGE_SIZES)))
        page_size = PAGE_SIZES[page_size]
    page_w, page_h = float(page_size[0]), float(page_size[1])
    
    if ext == ".pdf" and (n_rows, n_cols) == (1, 1):
        problems = _vector_pdf(items, filename, title)
    else:
        images, problems = _render_all(items, dpi)
        pages = _pages(images, n_rows * n_cols)
        if ext == ".pdf":
            _grid_pdf(pages, filename, title, n_rows, n_cols, page_w, page_h)
        else:
            _pptx(pages, filename, title, n_rows, n_cols, page_w, page_h)
            
    for p in problems:
        warnings.warn(p, CytoflowWarning)
    return problems


def _vector_pdf(items, filename, title):
    problems = []
    with matplotlib.rc_context({'pdf.fonttype' : 42}), PdfPages(filename) as pdf:
        for i, item in enumerate(items):
            try:
                fig = item.plot()
            except Exception as e:
                what = item.caption or "plot {}".format(i + 1)
                msg = e.args[-1] if isinstance(e, CytoflowError) and e.args else str(e)
                problems.append("{}: {}".format(what, msg))
                fig = Figure(figsize = (6, 3))
                FigureCanvasAgg(fig)
                fig.text(0.5, 0.5, "{}\n(couldn't be plotted: {})".format(what, msg),
                         ha = 'center', va = 'center', wrap = True)
                pdf.savefig(fig)
                continue
            try:
                extra = []
                text = " -- ".join(t for t in (title, item.caption) if t)
                if text:
                    extra.append(fig.text(0.5, 1.0, text, ha = 'center', va = 'bottom',
                                          fontsize = 'small'))
                pdf.savefig(fig, bbox_inches = 'tight', bbox_extra_artists = extra or None)
            finally:
                plt.close(fig)
    return problems


def _grid_pdf(pages, filename, title, n_rows, n_cols, page_w, page_h):
    margin = 0.45
    top = margin + (0.45 if title else 0.0)
    caption_h = 0.3
    gap = 0.2
    
    with matplotlib.rc_context({'pdf.fonttype' : 42}), PdfPages(filename) as pdf:
        for page_num, page in enumerate(pages):
            fig = Figure(figsize = (page_w, page_h))
            FigureCanvasAgg(fig)
            if title:
                fig.text(0.5, 1 - margin / page_h, title, ha = 'center', va = 'top',
                         fontsize = 14, fontweight = 'bold')
            fig.text(1 - margin / page_w, margin / 2 / page_h, 
                     "{} / {}".format(page_num + 1, len(pages)),
                     ha = 'right', va = 'center', fontsize = 8, color = '0.4')
            
            cells = _cells(n_rows, n_cols, margin, top, 
                           page_w - 2 * margin, page_h - top - margin, gap)
            for (png, caption), (x, y, w, h) in zip(page, cells):
                box_h = h - caption_h
                if png is not None:
                    img = plt.imread(io.BytesIO(png), format = 'png')
                    iw, ih = _fit(img.shape[1], img.shape[0], w, box_h)
                    ax = fig.add_axes([(x + (w - iw) / 2) / page_w,
                                       1 - (y + (box_h - ih) / 2 + ih) / page_h,
                                       iw / page_w, ih / page_h])
                    ax.imshow(img, interpolation = 'lanczos', aspect = 'auto')
                    ax.set_axis_off()
                    text, color, va, ty = caption, 'k', 'top', y + box_h + 0.05
                else:
                    text, color, va, ty = caption, '0.4', 'center', y + h / 2
                if text:
                    fig.text((x + w / 2) / page_w, 1 - ty / page_h, text, 
                             ha = 'center', va = va, fontsize = 9, color = color, wrap = True)
            pdf.savefig(fig)


def _pptx(pages, filename, title, n_rows, n_cols, page_w, page_h):
    try:
        from pptx import Presentation
        from pptx.util import Inches, Pt
        from pptx.enum.text import PP_ALIGN
    except ImportError as e:
        raise CytoflowError("Exporting to PowerPoint needs the 'python-pptx' "
                            "package. Install it with 'pip install python-pptx' "
                            "(or 'conda install -c conda-forge python-pptx').") from e
        
    prs = Presentation()
    prs.slide_width = Inches(page_w)
    prs.slide_height = Inches(page_h)
    blank = prs.slide_layouts[6]
    
    margin = 0.4
    top = margin + (0.7 if title else 0.0)
    caption_h = 0.4
    gap = 0.2
    
    def textbox(slide, x, y, w, h, text, size, bold = False, color = None):
        box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
        frame = box.text_frame
        frame.word_wrap = True
        para = frame.paragraphs[0]
        para.alignment = PP_ALIGN.CENTER
        run = para.add_run()
        run.text = text
        run.font.size = Pt(size)
        run.font.bold = bold
        return box
    
    for page in pages:
        slide = prs.slides.add_slide(blank)
        if title:
            textbox(slide, margin, margin * 0.6, page_w - 2 * margin, 0.6, title, 24, bold = True)
        cells = _cells(n_rows, n_cols, margin, top, 
                       page_w - 2 * margin, page_h - top - margin, gap)
        for (png, caption), (x, y, w, h) in zip(page, cells):
            box_h = h - caption_h
            if png is not None:
                pw, ph = _png_size(png)
                iw, ih = _fit(pw, ph, w, box_h)
                slide.shapes.add_picture(io.BytesIO(png), 
                                         Inches(x + (w - iw) / 2), Inches(y + (box_h - ih) / 2),
                                         Inches(iw), Inches(ih))
                if caption:
                    textbox(slide, x, y + box_h, w, caption_h, caption, 12)
            elif caption:
                textbox(slide, x, y + h / 2 - 0.4, w, 0.8, caption, 12)
                
    prs.save(filename)
