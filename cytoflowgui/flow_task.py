#!/usr/bin/env python3.8
# coding: latin-1
from traits.observation._trait_change_event import TraitChangeEvent

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
cytoflowgui.flow_task
---------------------

The main `pyface.tasks.task.Task`, and its associated `pyface.tasks.i_task_pane.ITaskPane`
and `envisage.plugin.Plugin`. 

- `FlowTask` -- the `pyface.tasks.task.Task` that allows the user to analyze flow data.

- `FlowTaskPane` -- the central `pyface.tasks.i_task_pane.ITaskPane` that contains the
  current `IWorkflowView`.
  
- `FlowTaskPlugin` -- the `envisage.plugin.Plugin` that provides the task factory and
  preferences pane.
"""

import os.path, webbrowser, pathlib, sys, warnings

import yaml.parser
from textwrap import dedent

import nbformat as nbf
from yapf.yapflib.yapf_api import FormatCode

from traits.api import Instance, Str, List, on_trait_change, provides, DelegatesTo, observe
from pyface.tasks.api import ITaskPane, TaskPane, Task, TaskLayout, PaneItem, VSplitter  # @UnresolvedImport
from pyface.tasks.action.api import SMenu, SMenuBar, SToolBar, TaskAction, TaskToggleGroup, SGroup
# from pyface.tasks.action.dock_pane_toggle_group import DockPaneToggleGroup, DockPaneToggleAction, ActionItem
from pyface.api import (FileDialog, ImageResource, AboutDialog,  # @UnresolvedImport
                        confirm, OK, YES, ConfirmationDialog, warning,  # @UnresolvedImport
                        error, MessageDialog, information)  # @UnresolvedImport
from pyface.qt import QtGui

from envisage.ui.tasks.api import TaskFactory
from envisage.ui.tasks.action.preferences_action import PreferencesAction
from envisage.api import Plugin

from .cytoflow_application import CytoflowApplication
from .workflow_pane import WorkflowDockPane
from .view_pane import ViewDockPane, PlotParamsPane
from .help_pane import HelpDockPane
from .experiment_pane import ExperimentBrowserDockPane
from .stats_pane import StatsDockPane
from .population_tree_pane import PopulationTreeDockPane
from .matplotlib_backend_local import FigureCanvasQTAggLocal
from .workflow import LocalWorkflow
from .workflow.undo import UndoManager
from .workflow_controller import WorkflowController
from .util import DefaultFileDialog
from .workflow.serialization import save_yaml, load_yaml
from .gate_interchange import import_gates, export_gates, is_flowjo_workspace
from .report_export import list_plots, ReportDialogModel

import cytoflow.utility as util
from cytoflow.operations.gatingml import flowjo_workspace_samples


@provides(ITaskPane)
class FlowTaskPane(TaskPane):
    """
    The center pane for the UI; contains the matplotlib canvas for plotting
    data views.
    """
    
    id = 'cytoflowgui.flow_task_pane'
    name = 'Cytometry Data Viewer'
    
    model = Instance(LocalWorkflow)
    """The shared `LocalWorkflow` model. Set by the task factory"""
    
    handler = Instance(WorkflowController)
    """The shared `WorkflowController`. Set by the task factory"""
    
    layout = Instance(QtGui.QVBoxLayout)                    # @UndefinedVariable
    """The main layout"""
    
    canvas = Instance(FigureCanvasQTAggLocal)
    """The local `matplotlib` canvas, an instance of `FigureCanvasQTAggLocal`"""
        
    def create(self, parent):      
        """Create a layout for the tab widget and the main view"""
        
        self.layout = layout = QtGui.QVBoxLayout()          # @UndefinedVariable
        self.control = QtGui.QWidget()                      # @UndefinedVariable
        self.control.setLayout(layout)
        
        tabs_ui = self.handler.edit_traits(view = 'selected_view_plot_name_view',
                                           context = self.model,
                                           kind = 'subpanel',
                                           parent = parent)
         
        self.layout.addWidget(tabs_ui.control) 
        
        # usually we would add the main plot here -- but a Qt widget
        # can only be part of one layout at a time.  so instead we
        # need to create that layout here, then dynamically add
        # the canvas when the task is activated (see activate(), below)
        
        
    def activate(self):
        """Activate this task by adding the canvas widget to the task pane layout"""
        if self.canvas.layout():
            self.canvas.layout().removeWidget(self.canvas)
        self.layout.addWidget(self.canvas)


class FlowTask(Task):
    """
    This class coordinates all the views and panels on the main model. Thus, you
    can think of it as an MVC controller - this is where all the UI handlers for
    things like the menu and menu bar buttons are.
    """
    
    id = "cytoflowgui.flow_task"
    """This tasks's GUID"""
    
    name = "Cytometry analysis"
    """This task's name"""
    
    # the main workflow instance.
    model = Instance(LocalWorkflow)
    """The shared `LocalWorkflow` model"""
    
    # the handler that connects it to various views
    handler = Instance(WorkflowController)
    """The shared `WorkflowController` controller"""
    
    # the application
    application = Instance(CytoflowApplication)
        
    # side panes
    workflow_pane = Instance(WorkflowDockPane)
    """The workflow dock pane"""
    
    view_pane = Instance(ViewDockPane)
    """The view configuration dock pane"""
    
    browser_pane = Instance(ExperimentBrowserDockPane)
    """The experiment browser dock pane"""
    
    help_pane = Instance(HelpDockPane)
    """The help dock pane"""
    
    plot_params_pane = Instance(PlotParamsPane)
    """The plot parameters dock pane"""

    stats_pane = Instance(StatsDockPane)
    """The live gate statistics dock pane"""

    population_tree_pane = Instance(PopulationTreeDockPane)
    """The population hierarchy tree dock pane"""
    
    undo_manager = Instance(UndoManager)
    """Undo / redo history"""
    
    can_undo = DelegatesTo('undo_manager')
    can_redo = DelegatesTo('undo_manager')
    
    # plugin lists, to setup the interface
    op_plugins = DelegatesTo('handler')
    """The operation plugins"""
    
    view_plugins = DelegatesTo('handler')
    """The view plugins"""
    
    menu_bar = SMenuBar(SMenu(TaskAction(name='New',
                                         method='on_new',
                                         accelerator='Ctrl+N'),
                              TaskAction(name='Open...',
                                         method='on_open',
                                         accelerator='Ctrl+O'),
                              TaskAction(name='Save',
                                         #image='save', 
                                         method='on_save',
                                         accelerator='Ctrl+S'),
                              TaskAction(name='Save As...',
                                         method='on_save_as',
                                         accelerator='Ctrl+Shift+S'),
                              TaskAction(name='Save Plot...',
                                         method='on_export',
                                         accelerator='Ctrl+Shift+E'),
                              TaskAction(name='Export Jupyter notebook...',
                                         method='on_notebook',
                                         accelerator='Ctrl+I'),   
                              TaskAction(name='Import gates...',
                                         method='on_import_gates'),
                              TaskAction(name='Export gates (Gating-ML)...',
                                         method='on_export_gates'),
                              TaskAction(name='Export report (PDF / PowerPoint)...',
                                         method='on_export_report',
                                         accelerator='Ctrl+Shift+R'),
                              TaskAction(name='Save gate as template...',
                                         method='on_save_gate_template',
                                         accelerator='Ctrl+T'),
                              TaskAction(name='Load gate template...',
                                         method='on_load_gate_template',
                                         accelerator='Ctrl+Shift+T'),
                              TaskAction(name='Channel scale defaults...',
                                         method='on_channel_scale_defaults'),
                              PreferencesAction(),
                              id='File', name='&File'),
                        SMenu(TaskAction(name='Undo',
                                         method='on_undo',
                                         accelerator='Ctrl+Z',
                                         enabled_name='can_undo'),
                              TaskAction(name='Redo',
                                         method='on_redo',
                                         accelerator='Ctrl+Shift+Z',
                                         enabled_name='can_redo'),
                              id='Edit', name='&Edit'),
                        SMenu(TaskToggleGroup(),
                              id = 'View', name = '&View'),
                        SMenu(TaskAction(name = 'Online documentation...',
                                         method = 'on_docs'),
                              TaskAction(name='About...',
                                         method='on_about'),
                              id="Help", name ="&Help"))
    """The menu bar schema"""
    
    tool_bars = [ SToolBar(SGroup(TaskAction(method='on_new',
                                             name = "New",
                                             tooltip='New workflow',
                                             image=ImageResource('new')),
                                  TaskAction(method='on_open',
                                             name = "Open",
                                             tooltip='Open a file',
                                             image=ImageResource('open')),
                                  TaskAction(method='on_save',
                                             name = "Save",
                                             tooltip='Save the current file',
                                             image=ImageResource('save')),
                                  TaskAction(method='on_export',
                                             name = "Save Plot",
                                             tooltip='Save the current plot',
                                             image=ImageResource('export')),
                                  TaskAction(method='on_notebook',
                                              name='Notebook',
                                              tooltip="Export to a Jupyter notebook...",
                                              image=ImageResource('jupyter')),
                                  PreferencesAction(image = ImageResource('prefs'))),
                           SGroup(TaskAction(method = 'on_toggle_workflow',
                                             name = 'Workflow',
                                             tooltip = 'Toggle workflow pane',
                                             image = ImageResource('workflow')),
                                  TaskAction(method = 'on_toggle_view',
                                                    name = 'View Properties',
                                                    tooltip = 'Toggle view properties pane',
                                                    image = ImageResource('view_pane')),
                                  TaskAction(method = 'on_toggle_plot_params',
                                                    name = 'Plot Parameters',
                                                    tooltip = 'Toggle plot parameters pane',
                                                    image = ImageResource('plot_params')),
                                  TaskAction(method = 'on_toggle_browser',
                                                    name = 'Browser',
                                                    tooltip = 'Toggle experiment browser pane',
                                                    image = ImageResource('browser')),
                                  TaskAction(method = 'on_toggle_help',
                                                    name = 'Help',
                                                    tooltip = 'Toggle help pane',
                                                    image = ImageResource('help'))))]
    """The tool bar schema"""

    
    filename = Str
    """
    The file to save to if the user clicks "save" and has already clicked
    "open" or "save as".
    """
        
    def initialized(self):
        """
        Called when the task is about to be activated in a TaskWindow for
        the first time.  Called before `activated`. If a filename was given
        on the command line, load it.
        """
        
        if self.filename:
            self.open_file(self.filename)

        
    def activated(self):
        """
        Called after the task has been activated in a TaskWindow.  Initialize the
        model with an import operation and activate the central pane.
        """

        # add the import op
        if not self.model.workflow:
            self.handler.add_operation('cytoflow.operations.import')

        self.window.central_pane.activate()

        self.model.modified = False
        self.undo_manager.clear()

        # Set up progress bar: watch WorkflowItem status changes and
        # any operation-level status messages (e.g. UMAP log messages)
        self.model.observe(self._on_workflow_status_changed,
                           'workflow:items:status',
                           dispatch='ui')
        self.model.observe(self._on_workflow_status_changed,
                           'workflow:items:operation:+status',
                           dispatch='ui')
        # Reflect the current state immediately (in case a file was loaded)
        self._on_workflow_status_changed()

        if self.application.preferences_helper.first_run:
            self.application.preferences_helper.first_run = False
            
            msg = dedent("""
            Welcome to <strong>Cytoflow-GCB!</strong>
            <p>
            Get started: click <strong>Set up experiment</strong> in the Import step
            to load your FCS files. Then add gates from the toolbar on the left.
            <p>
            You can find the manual online at
            <a href="https://cytoflow.readthedocs.io/en/stable/user_manual/user_manual.html">cytoflow.readthedocs.io</a>.
            <p>
            Happy flowing!
            """)

            dialog = MessageDialog(title = "Welcome to Cytoflow-GCB!",
                                   message = msg,
                                   text_format = "rich",
                                   severity = "information")
            dialog.open()
    
    def _default_layout_default(self):
        """
        Returns the default layout for the dock panes.
        """
        
        return TaskLayout(left = VSplitter(PaneItem("cytoflowgui.workflow_pane", width = 280),
                                           PaneItem("cytoflowgui.experiment_pane", width = 280),
                                           PaneItem("cytoflowgui.help_pane", width = 280, height = 280)),
                          right = VSplitter(PaneItem("cytoflowgui.view_traits_pane", width = 340),
                                            PaneItem("cytoflowgui.params_pane", width = 340),
                                            PaneItem("cytoflowgui.stats_pane"),
                                            PaneItem("cytoflowgui.population_tree_pane")),
                          top_left_corner = 'left',
                          bottom_left_corner = 'left',
                          top_right_corner = 'right',
                          bottom_right_corner = 'right')
     
    def create_central_pane(self):  
        """
        Initialize the toolbar image size and return the central `FlowTaskPane`
        """
    
        # this isn't really the right place for this, but it's the only
        # place control passes back to user code before the toolbar
        # is created.
    
        self.tool_bars[0].image_size = (40, 40)
        self.tool_bars[0].show_tool_names = self.application.preferences_helper.show_toolbar_names
    
        return FlowTaskPane(canvas = self.application.canvas,
                            model = self.application.model, 
                            handler = self.handler)
             
    def create_dock_panes(self):
        """
        Create and initialize the dock panes
        """
                
        self.workflow_pane = WorkflowDockPane(model = self.model, 
                                              handler = self.handler,
                                              plugins = self.op_plugins,
                                              task = self)
        
        self.view_pane = ViewDockPane(model = self.model,
                                      handler = self.handler,
                                      plugins = self.view_plugins,
                                      task = self)
        
        self.plot_params_pane = PlotParamsPane(model = self.model,
                                               handler = self.handler)
        
        self.browser_pane = ExperimentBrowserDockPane(model = self.model,
                                                      handler = self.handler,
                                                      task = self)

        self.help_pane = HelpDockPane(model = self.model,
                                      view_plugins = self.view_plugins,
                                      op_plugins = self.op_plugins,
                                      task = self)

        self.stats_pane = StatsDockPane(model = self.model,
                                        handler = self.handler)

        self.population_tree_pane = PopulationTreeDockPane(model = self.model,
                                                           handler = self.handler)

        return [self.workflow_pane, self.view_pane, self.plot_params_pane,
                self.browser_pane, self.help_pane, self.stats_pane,
                self.population_tree_pane]
        
    def on_new(self):
        """
        Create a new workflow when the "New..." button or menu item is clicked
        """
        
        if self.model.modified:
            ret = confirm(parent = None,
                          message = "Are you sure you want to discard the current workflow?",
                          title = "Clear workflow?")
            
            if ret != YES:
                return
            
        self.filename = ""
        self.window.title = "Cytoflow-GCB"
        
        # clear the workflow
        self.model.workflow = []
        
        # add the import op
        self.handler.add_operation('cytoflow.operations.import')
        
        self.model.modified = False
        self.undo_manager.clear()
     
        
    def on_open(self):
        """ 
        Shows a dialog to open a file.
        """
        
        if self.model.modified:
            ret = confirm(parent = None,
                          message = "Are you sure you want to discard the current workflow?",
                          title = "Clear workflow?")
            
            if ret != YES:
                return
        
        dialog = FileDialog(parent = self.window.control, 
                            action = 'open',
                            wildcard = (FileDialog.create_wildcard("Cytoflow workflow", "*.flow") + ';' +  #@UndefinedVariable  
                                        FileDialog.create_wildcard("All files", "*"))) #@UndefinedVariable  
        if dialog.open() == OK:
            self.open_file(dialog.path)
            self.filename = dialog.path
            self.window.title = "Cytoflow-GCB - " + self.filename
            

    def open_file(self, path):
        """
        Load a new workflow from a file
        """
        
        new_workflow = []
        with warnings.catch_warnings(record = True) as w:
            try:

                new_workflow = load_yaml(path)

                if w:
                    warning(None, w[-1].message.__str__())

                # a few things to take care of when reloading.
                # we do this in the try block to catch people who
                # load valid YAML files that aren't from cytoflow.

                for wi_idx, wi in enumerate(new_workflow):

                    # get wi lock
                    wi.lock.acquire()

                    # clear the wi status
                    wi.status = "loading"

                    # re-link the linked list.
                    if wi_idx > 0:
                        wi.previous_wi = new_workflow[wi_idx - 1]

                    if wi_idx < len(new_workflow) - 1:
                        wi.next_wi = new_workflow[wi_idx + 1]

            except yaml.parser.ParserError as e:
                for wi in new_workflow:
                    try: wi.lock.release()
                    except RuntimeError: pass
                error(None,
                      "Parser error loading {} -- is it a Cytoflow-GCB file?\n\n{}"
                      .format(path, str(e)))
                return
            except Exception as e:
                for wi in new_workflow:
                    try: wi.lock.release()
                    except RuntimeError: pass
                error(None,
                      "{} loading {} -- is it a Cytoflow-GCB file?\n\n{}"
                      .format(e.__class__.__name__, path, str(e)))
                return
        
        # are we just running a smoke test?
        if 'startup_test' in new_workflow[0].metadata:
            def quit_app(app):
                app.exit(force = True)
                
            from pyface.timer.api import do_after  # @UnresolvedImport
            do_after(5*1000, quit_app, self.application)
            return
            
        # check that the FCS files are all there
        
        wi = new_workflow[0]
        assert(wi.operation.id == "cytoflow.operations.import")
        missing_tubes = 0
        for tube in wi.operation.tubes:
            file = pathlib.Path(tube.file)
            if not file.exists():
                missing_tubes += 1
                
        if missing_tubes == len(wi.operation.tubes):
            warning(self.window.control,
                    "Cytoflow-GCB couldn't find any of the FCS files from that "
                    "workflow.  If they've been moved, please open one FCS "
                    "file to show Cytoflow-GCB where they've been moved to.")
            
            dialog = FileDialog(parent = self.window.control, 
                                action = 'open',
                                wildcard = (FileDialog.create_wildcard("FCS files", "*.fcs *.lmd")))  # @UndefinedVariable
            
            if dialog.open() == OK:
                # find the "best" file match -- ie, the one with the longest
                # tail match
                fcs_path = pathlib.Path(dialog.path).parts
                best_path_len = -1
                                
                for tube in wi.operation.tubes:
                    tube_path = pathlib.Path(tube.file).parts
                    
                    for i in range(len(fcs_path)):
                        if list(reversed(fcs_path))[:i] == list(reversed(tube_path))[:i] and i > best_path_len:
                            best_path_len = i
                            
                if best_path_len >= 0:
                    for tube in wi.operation.tubes:
                        tube_path = pathlib.Path(tube.file).parts
                        new_path = fcs_path[:-1 * best_path_len] + tube_path[-1 * best_path_len :]
                        tube.file = str(pathlib.Path(*new_path))
                        
        elif missing_tubes > 0:
            warning(self.window.control,
                    "Cytoflow-GCB couldn't find some of the FCS files from that "
                    "workflow.  You'll need to re-load them from the Import "
                    "operation.")

        # replace the current workflow with the one we just loaded
        self.model.workflow = new_workflow
        self.model.modified = False
        self.undo_manager.clear()

        for wi in self.model.workflow:
            wi.lock.release()
            
        self.model.run_all()
        
    def on_save(self):
        """ Save the workflow to the current filename  """
        if self.filename:
            save_yaml(self.model.workflow, self.filename)
            self.model.modified = False
        else:
            self.on_save_as()
            
    def on_save_as(self):
        """ Save the workflow to a different file"""
        dialog = DefaultFileDialog(parent = self.window.control,
                                   action = 'save as', 
                                   default_suffix = "flow",
                                   wildcard = (FileDialog.create_wildcard("Cytoflow workflow", "*.flow") + ';' + #@UndefinedVariable  
                                               FileDialog.create_wildcard("All files", "*")))                    #@UndefinedVariable  
        
        if dialog.open() == OK:
            save_yaml(self.model.workflow, dialog.path)
            self.filename = dialog.path
            self.model.modified = False
            self.window.title = "Cytoflow-GCB - " + self.filename
            
    @on_trait_change('model.modified', post_init = True)
    def _on_model_modified(self, val):
        if val:
            if not self.window.title.endswith("*"):
                self.window.title += "*"
        else:
            if self.window.title.endswith("*"):
                self.window.title = self.window.title[:-1]

    def _on_workflow_status_changed(self, event=None):
        """
        Show or hide the indeterminate progress bar in the window status bar
        whenever any WorkflowItem is in the "applying" or "estimating" state.
        Also displays operation-level progress messages (e.g. from UMAP) when
        available.
        """
        from pyface.qt import QtWidgets
        win = self.window
        if win is None or win.control is None:
            return
        statusbar = win.control.statusBar()

        # Create the progress bar widget once and cache it on the window control
        pb = getattr(win.control, '_gcb_progress_bar', None)
        if pb is None:
            pb = QtWidgets.QProgressBar()
            pb.setMaximumWidth(200)
            pb.setMaximumHeight(16)
            pb.setRange(0, 0)   # indeterminate / busy mode
            pb.setTextVisible(False)
            statusbar.addPermanentWidget(pb)
            win.control._gcb_progress_bar = pb

        busy = any(item.status in ('applying', 'estimating')
                   for item in self.model.workflow)
        pb.setVisible(busy)

        if busy:
            # Show an operation-level progress message when one is available
            # (e.g. UMAP / AutoSpill / PhenoGraph log messages routed back
            # via the UMAPWorkflowOp.status Str trait)
            msg = "Working…"
            for item in self.model.workflow:
                if item.status in ('applying', 'estimating'):
                    op_status = getattr(item.operation, 'status', None)
                    if op_status:
                        msg = op_status
                        break
            statusbar.showMessage(msg)
        else:
            statusbar.clearMessage()

    def _undo_manager_default(self):
        return UndoManager(model = self.model)
    
    def on_undo(self):
        """Undo the last change to the workflow"""
        self.undo_manager.undo()
        self.model.modified = True
        
    def on_redo(self):
        """Redo the last undone change"""
        self.undo_manager.redo()
        self.model.modified = True
        
    def on_export(self):
        """
        Switch to the `ExportTask` task
        """
        
        task = next(x for x in self.window.tasks if x.id == 'cytoflowgui.export_task')
        self.window.activate_task(task)        
        
            
    def on_notebook(self):
        """
        Shows a dialog to export the workflow to an Jupyter notebook
        """

        dialog = DefaultFileDialog(parent = self.window.control,
                                   action = 'save as',
                                   default_suffix = "ipynb",
                                   wildcard = (FileDialog.create_wildcard("Jupyter notebook", "*.ipynb") + ';' + #@UndefinedVariable  
                                               FileDialog.create_wildcard("All files", "*")))  # @UndefinedVariable
        if dialog.open() == OK:
            self.save_notebook(self.model.workflow, dialog.path)

            
    def on_import_gates(self):
        """
        Shows a dialog to import gates from a Gating-ML file or a FlowJo
        workspace, then adds them to the end of the workflow.
        """
        
        dialog = FileDialog(parent = self.window.control, 
                            action = 'open',
                            wildcard = (FileDialog.create_wildcard("Gating-ML or FlowJo workspace", "*.xml *.gml *.gatingml *.wsp") + ';' + #@UndefinedVariable
                                        FileDialog.create_wildcard("All files", "*")))  # @UndefinedVariable
        if dialog.open() != OK:
            return
        
        path = dialog.path
        sample = group = None
        
        try:
            if is_flowjo_workspace(path):
                # FlowJo has a set of gates for each sample; we have one.
                samples, groups = flowjo_workspace_samples(path)
                choices = (["Sample: " + s for s in samples] + 
                           ["Group: " + g for g in groups])
                if not choices:
                    raise util.CytoflowError("There aren't any samples in {}".format(path))
                if len(samples) > 1:
                    choice, ok = QtGui.QInputDialog.getItem(self.window.control,
                                                            "Import gates",
                                                            "Import the gates from which sample or group?",
                                                            choices, 0, False)
                    if not ok:
                        return
                    if choice.startswith("Sample: "):
                        # there may be two samples with the same name
                        sample = samples[choices.index(choice)]
                    else:
                        group = choice[len("Group: "):]
                        
            new_items, warns = import_gates(self.model, path, sample = sample, group = group)
        except Exception as e:
            error(None, "Couldn't import gates from {}:\n\n{}".format(path, str(e)))
            return
        
        if new_items:
            self.model.modified = True
            self.model.selected = new_items[-1]
            if new_items[-1].default_view:
                self.handler.activate_view(new_items[-1].default_view.id)
            
        msg = "Imported {} operation{} from {}".format(len(new_items),
                                                       "" if len(new_items) == 1 else "s",
                                                       path)
        if warns:
            warning(None, msg + "\n\n" + "\n\n".join(warns))
        else:
            information(None, msg)

            
    def on_export_gates(self):
        """
        Shows a dialog to export the workflow's gates to a Gating-ML file
        """
        
        dialog = DefaultFileDialog(parent = self.window.control,
                                   action = 'save as',
                                   default_suffix = "xml",
                                   wildcard = (FileDialog.create_wildcard("Gating-ML", "*.xml") + ';' + #@UndefinedVariable  
                                               FileDialog.create_wildcard("All files", "*")))  # @UndefinedVariable
        if dialog.open() != OK:
            return
        
        try:
            exported, warns = export_gates(self.model, dialog.path)
        except Exception as e:
            error(None, "Couldn't export gates to {}:\n\n{}".format(dialog.path, str(e)))
            return
        
        msg = "Exported {} gate{} to {}".format(len(exported),
                                                "" if len(exported) == 1 else "s",
                                                dialog.path)
        if warns:
            warning(None, msg + "\n\n" + "\n\n".join(warns))
        else:
            information(None, msg)

    
    def on_export_report(self):
        """
        Shows a dialog to export several of the workflow's plots as one 
        report: a PDF, or a PowerPoint presentation.
        """
        
        plots = list_plots(self.model.workflow)
        if not plots:
            information(None, "There are no plots to export yet. Add a view (or "
                              "select an operation with a plot) first.")
            return
        
        labels = [label for _, _, label in plots]
        model = ReportDialogModel(labels = labels, selected = list(labels))
        if not model.edit_traits(kind = 'livemodal').result:
            return
        
        chosen = [(wi, view, label.split(": ", 1)[-1]) 
                  for wi, view, label in plots if label in model.selected]
        if not chosen:
            information(None, "No plots were chosen.")
            return
        
        dialog = DefaultFileDialog(parent = self.window.control,
                                   action = 'save as',
                                   default_suffix = "pdf",
                                   wildcard = (FileDialog.create_wildcard("PDF", "*.pdf") + ';' + #@UndefinedVariable  
                                               FileDialog.create_wildcard("PowerPoint", "*.pptx")))  # @UndefinedVariable
        if dialog.open() != OK:
            return
        
        path = dialog.path
        if not path.lower().endswith((".pdf", ".pptx")):
            path += ".pptx" if dialog.wildcard_index == 1 else ".pdf"
            
        # the report is made in the background; _on_report_result shows the
        # outcome when it's done
        self.model.export_report(chosen, path, **model.options(path))
        
    @observe('model:report_result', dispatch = 'ui')
    def _on_report_result(self, event):
        if not event.new:
            return
        succeeded, message = event.new
        if not succeeded:
            error(None, "Couldn't export the report:\n\n{}".format(message))
        elif "could not be made" in message:
            warning(None, message)
        else:
            information(None, message)

    
    def on_save_gate_template(self):
        """
        Save the current operation's gate parameters to a JSON template file
        so it can be re-applied to future experiments.
        """
        import json, os, pathlib

        if not self.model.selected or not self.model.selected.operation:
            information(None, "Select a gate operation first.")
            return

        op = self.model.selected.operation

        # Collect all traits that are involved in apply/estimate
        apply_traits = [n for n in op.trait_names()
                        if (op.trait(n).apply or op.trait(n).estimate)
                        and not op.trait(n).transient
                        and not op.trait(n).status]

        def _jsonify(v):
            if isinstance(v, (str, int, float, bool)) or v is None:
                return v
            if isinstance(v, (list, tuple)):
                return [_jsonify(x) for x in v]
            if isinstance(v, dict):
                return {str(k): _jsonify(vv) for k, vv in v.items()}
            return str(v)

        payload = {
            '_op_id': op.id,
            '_op_friendly_id': op.friendly_id,
        }
        for name in apply_traits:
            payload[name] = _jsonify(getattr(op, name))

        # Ensure template directory exists
        template_dir = pathlib.Path.home() / '.cytoflowgcb' / 'templates'
        template_dir.mkdir(parents=True, exist_ok=True)

        dialog = DefaultFileDialog(
            parent=self.window.control,
            action='save as',
            default_directory=str(template_dir),
            default_suffix='gatetemplate',
            wildcard=FileDialog.create_wildcard("Gate templates", "*.gatetemplate") + ';' +
                     FileDialog.create_wildcard("All files", "*"))
        if dialog.open() != OK:
            return

        path = dialog.path
        if not path.endswith('.gatetemplate'):
            path += '.gatetemplate'

        with open(path, 'w') as f:
            json.dump(payload, f, indent=2)

        information(None, "Gate template saved to\n{}".format(path))

    def on_load_gate_template(self):
        """
        Load a gate template JSON file and add the corresponding operation
        after the currently selected workflow item.
        """
        import json, pathlib

        template_dir = pathlib.Path.home() / '.cytoflowgcb' / 'templates'
        template_dir.mkdir(parents=True, exist_ok=True)

        dialog = DefaultFileDialog(
            parent=self.window.control,
            action='open',
            default_directory=str(template_dir),
            wildcard=FileDialog.create_wildcard("Gate templates", "*.gatetemplate") + ';' +
                     FileDialog.create_wildcard("All files", "*"))
        if dialog.open() != OK:
            return

        try:
            with open(dialog.path) as f:
                payload = json.load(f)
        except Exception as exc:
            error(None, "Couldn't read template:\n\n{}".format(exc))
            return

        op_id = payload.pop('_op_id', None)
        payload.pop('_op_friendly_id', None)

        # Find the plugin that handles this operation
        plugin = next(
            (p for p in self.op_plugins
             if getattr(p, 'operation_id', None) == op_id),
            None)
        if plugin is None:
            error(None, "No plugin found for operation id '{}'.\n"
                        "Make sure the operation is available in this version."
                        .format(op_id))
            return

        new_op = plugin.get_operation()

        for name, value in payload.items():
            try:
                setattr(new_op, name, value)
            except Exception:
                pass  # skip traits that don't match

        from cytoflowgui.workflow.workflow_item import WorkflowItem
        new_wi = WorkflowItem(operation=new_op)

        wf = self.model.workflow
        if self.model.selected and self.model.selected in wf:
            idx = wf.index(self.model.selected) + 1
        else:
            idx = len(wf)

        wf.insert(idx, new_wi)
        self.model.selected = new_wi

    def on_prefs(self):
        """
        Event handler for the "preferences" button -- currently, do nothing
        """
        pass


    def on_channel_scale_defaults(self):
        """
        Open a simple dialog for editing per-channel default scales.  The
        mapping is stored in user preferences and auto-applied when a
        ScaleSettingsOp is inserted after an Import.
        """
        from pyface.qt import QtGui, QtCore

        prefs = self.application.preferences_helper
        current = prefs.get_channel_scale_defaults()

        # Build a simple two-column table dialog using Qt directly
        SCALES = ['linear', 'log', 'logicle', 'arcsinh', 'hlog']

        dlg = QtGui.QDialog(self.window.control)
        dlg.setWindowTitle("Channel Scale Defaults")
        dlg.setMinimumWidth(420)
        layout = QtGui.QVBoxLayout(dlg)

        label = QtGui.QLabel(
            "Set the default scale for each channel.  These are applied\n"
            "automatically when you open a new experiment.\n"
            "Leave blank or use 'logicle' to keep the automatic default.")
        layout.addWidget(label)

        table = QtGui.QTableWidget(0, 2, dlg)
        table.setHorizontalHeaderLabels(["Channel", "Scale"])
        table.horizontalHeader().setStretchLastSection(True)
        table.setSelectionBehavior(QtGui.QAbstractItemView.SelectionBehavior.SelectRows)

        def _add_row(channel='', scale='logicle'):
            row = table.rowCount()
            table.insertRow(row)
            table.setItem(row, 0, QtGui.QTableWidgetItem(channel))
            combo = QtGui.QComboBox()
            combo.addItems(SCALES)
            idx = combo.findText(scale)
            if idx >= 0:
                combo.setCurrentIndex(idx)
            table.setCellWidget(row, 1, combo)

        for ch, sc in sorted(current.items()):
            _add_row(ch, sc)

        layout.addWidget(table)

        btn_row = QtGui.QHBoxLayout()
        add_btn = QtGui.QPushButton("Add row")
        del_btn = QtGui.QPushButton("Remove selected")
        add_btn.clicked.connect(lambda: _add_row())
        del_btn.clicked.connect(lambda: table.removeRow(table.currentRow()))
        btn_row.addWidget(add_btn)
        btn_row.addWidget(del_btn)
        btn_row.addStretch()
        layout.addLayout(btn_row)

        buttons = QtGui.QDialogButtonBox(
            QtGui.QDialogButtonBox.StandardButton.Ok |
            QtGui.QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        layout.addWidget(buttons)

        if dlg.exec() == QtGui.QDialog.DialogCode.Accepted:
            new_map = {}
            for row in range(table.rowCount()):
                ch_item = table.item(row, 0)
                combo = table.cellWidget(row, 1)
                if ch_item and ch_item.text().strip() and combo:
                    new_map[ch_item.text().strip()] = combo.currentText()
            import json
            prefs.channel_scale_defaults_json = json.dumps(new_map)

    def on_docs(self):
        """
        Event handler for the "Online documentation...." menu item.
        Opens a webbrowser to the online manual at readthedocs.io
        """

        webbrowser.open_new_tab("https://cytoflow.readthedocs.io/en/stable/user_manual/user_manual.html")

    
    def on_toggle_workflow(self):
        """Toggle the visibility of the workflow pane"""
        self.workflow_pane.visible = not self.workflow_pane.visible
        
        
    def on_toggle_view(self):
        """Toggle the visibility of the view properties pane"""
        self.view_pane.visible = not self.view_pane.visible
        
        
    def on_toggle_plot_params(self):
        """Toggle the visibility of the plot parameters pane"""
        self.plot_params_pane.visible = not self.plot_params_pane.visible
        
        
    def on_toggle_browser(self):
        """Toggle the visiblity of the experiment browser pane"""
        self.browser_pane.visible = not self.browser_pane.visible
        
        
    def on_toggle_help(self):
        """Toggle the visiblity of the help pane"""
        self.help_pane.visible = not self.help_pane.visible
        
    
    def _get_package_versions(self):    
        
        # if you add anything here, you should also make sure the package is
        # added to hook-yapf.py
        
        from importlib.metadata import version as get_version
        cf_version = get_version('cytoflow')
        from fcsparser.fcsparser import __version__ as fcs_version
        np_version = get_version('numpy')
        pd_version = get_version('pandas')
        mpl_version = get_version('matplotlib')
        btl_version = get_version('bottleneck')
        nxp_version = get_version('numexpr')
        scipy_version = get_version('scipy')
        skl_version = get_version('scikit-learn')
        sns_version = get_version('seaborn')
        stats_version = get_version('statsmodels')
        natsort_version = get_version('natsort')
        numba_version = get_version('numba')
        trt_version = get_version('traits')
        trt_ui_version = get_version('traitsui')
        pyf_version = get_version('pyface')
        env_version = get_version('envisage')
        nb_version = get_version('nbformat')
        yapf_version = get_version('yapf')
        from yaml import __version__ as yaml_version
        
        return {"python" : sys.version,
                "cytoflow" : cf_version,
                "fcsparser" : fcs_version,
                "numpy" : np_version,
                "pandas" : pd_version,
                "matplotlib" : mpl_version,
                "bottleneck" : btl_version,
                "numexpr" : nxp_version,
                "scipy" : scipy_version,
                "scikit-learn" : skl_version,
                "seaborn" : sns_version,
                "statsmodels" : stats_version,
                "natsort" : natsort_version,
                "numba" : numba_version,
                "traits" : trt_version,
                "traitsui" : trt_ui_version,
                "pyface" : pyf_version,
                "envisage" : env_version,
                "nbformat" : nb_version,
                "yapf" : yapf_version,
                "yaml" : yaml_version}
        
        
    def on_about(self):
        """
        Event handler for "About..." menu item.  Shows an about dialog box.
        """
        
        versions = self._get_package_versions()
        text = ["<b>Cytoflow-GCB {0}</b>".format(versions['cytoflow']),
                "<p>Based on <a href=https://github.com/cytoflow/cytoflow>Cytoflow</a> "
                "by Brian Teague. Extended for the GCB lab.</p>",
                "<p>"]

        ver_text = ["{0} {1}".format(key, value) for key, value in versions.items()]
        text.extend(ver_text)

        text.extend(["<p>Icons from the <a href=http://tango.freedesktop.org>Tango Desktop Project</a>",
                "<a href=http://www.freepik.com/free-photos-vectors/background>App icon from Starline - Freepik.com</a>",
                "Cuvette image from Wikimedia Commons user <a href=http://commons.wikimedia.org/wiki/File:Hellma_Large_cone_cytometry_cell.JPG>HellmaUSA</a></p>"])

        dialog = AboutDialog(text = text,
                             parent = self.window.control,
                             title = "About Cytoflow-GCB",
                             image = ImageResource('cuvette'),
                             additions = text)
        dialog.open()

    # Jupyter notebook serialization
    def save_notebook(self, workflow, path):
        """
        Saves the workflow to a Jupyter notebook.
        """
        
        nb = nbf.v4.new_notebook()
        
        # todo serialize here
        header = dedent("""\
            import cytoflow as flow""")
        nb['cells'].append(nbf.v4.new_code_cell(header))
            
        for i, wi in enumerate(workflow):
            try:
                code = wi.operation.get_notebook_code(i)
                code = FormatCode(code, style_config = 'pep8')[0]
            except Exception as e:
                error(parent = None,
                      message = "Had trouble serializing the {} operation:\n{}"
                                .format(wi.operation.friendly_id,
                                        repr(e)))
            
            nb['cells'].append(nbf.v4.new_code_cell(code))
                        
            for view in wi.views:
                try:
                    code = view.get_notebook_code(i)
                    code = FormatCode(code, style_config = 'pep8')[0]
                except Exception as e:
                    error(parent = None,
                          message = "Had trouble serializing the {} view of the {} operation:\n{}"
                                     .format(view.friendly_id, 
                                             wi.operation.friendly_id,
                                             repr(e)))
                
                nb['cells'].append(nbf.v4.new_code_cell(code))
                
        with open(path, 'w') as f:
            nbf.write(nb, f)
            
    @observe('application.preferences_helper.show_toolbar_names', post_init = True)
    def show_toolbar_names(self, _):
        self.window.tool_bar_managers[0].show_tool_names = self.application.preferences_helper.show_toolbar_names
        self.window._update_tool_bar_managers(TraitChangeEvent(object = None,
                                                               name = None,
                                                               old = None,
                                                               new = True))

        
class FlowTaskPlugin(Plugin):
    """
    An Envisage plugin wrapping FlowTask
    """

    # Extension point IDs.
    PREFERENCES       = 'envisage.preferences'
    PREFERENCES_PANES = 'envisage.ui.tasks.preferences_panes'
    TASKS             = 'envisage.ui.tasks.tasks' 
    
    #### 'IPlugin' interface ##################################################

    # The plugin's unique identifier.
    id = 'cytoflow'

    # The plugin's name (suitable for displaying to the user).
    name = 'Cytoflow-GCB'

    ###########################################################################
    # Protected interface.
    ###########################################################################

    preferences = List(contributes_to = PREFERENCES)
    def _preferences_default(self):
        filename = os.path.join(os.path.dirname(__file__), 'preferences.ini')
        return [ 'file://' + filename ]
    
    preferences_panes = List(contributes_to = PREFERENCES_PANES)
    def _preferences_panes_default(self):
        from .preferences import CytoflowPreferencesPane
        return [CytoflowPreferencesPane]

    tasks = List(contributes_to = TASKS)
    def _tasks_default(self):
        return [TaskFactory(id = 'cytoflowgui.flow_task',
                            name = 'Cytometry analysis',
                            factory = lambda **x: FlowTask(application = self.application,
                                                           model = self.application.model,
                                                           handler = self.application.controller,
                                                           filename = self.application.filename,
                                                           **x))]
