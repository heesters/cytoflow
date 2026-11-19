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
cytoflowgui.run
---------------

The entry-point for the GUI -- sets up and starts the remote process, configures
logging, loads the Envisage plugins, and starts the GUI loop.
"""

try:
    import faulthandler       # @UnresolvedImport
    faulthandler.enable()     # @UndefinedVariable
except:
    # if there's no console, this fails
    pass

import sys, multiprocessing, logging, traceback, threading, argparse

logger = logging.getLogger(__name__)

def log_notification_handler(_, trait_name, old, new):
    """
    Exception handler for traits notifications
    """
    
    (exc_type, exc_value, tb) = sys.exc_info()
    logger.debug('Exception occurred in traits notification '
                  'handler for object: %s, trait: %s, old value: %s, '
                  'new value: %s.\n%s\n' % ( object, trait_name, old, new,
                  ''.join( traceback.format_exception(exc_type, exc_value, tb) ) ) )

    err_string = traceback.format_exception_only(exc_type, exc_value)[0]
    err_loc = traceback.format_tb(tb)[-1]
    err_ctx = threading.current_thread().name
    
    logging.error("Error: {0}\nLocation: {1}Thread: {2}" \
                  .format(err_string, err_loc, err_ctx) )
    
def log_excepthook(typ, val, tb):
    """Exception handler for global exceptions"""
    tb_str = "".join(traceback.format_tb(tb))
    logger.debug("Global exception: {0}\n{1}: {2}".format(tb_str, typ, val))
    
    tb_str = traceback.format_tb(tb)[-1]
    logging.error("Error: {0}: {1}\nLocation: {2}Thread: Main"
                  .format(typ, val, tb_str))
    
                         
_MATPLOTLIB_DARK_PARAMS = {
    'figure.facecolor':  '#1e1e1e',
    'axes.facecolor':    '#1e1e1e',
    'savefig.facecolor': '#1e1e1e',
    'text.color':        '#dcdcdc',
    'axes.labelcolor':   '#dcdcdc',
    'xtick.color':       '#dcdcdc',
    'ytick.color':       '#dcdcdc',
    'axes.edgecolor':    '#555555',
    'grid.color':        '#444444',
    'legend.facecolor':  '#2d2d2d',
    'legend.edgecolor':  '#555555',
}


def _detect_dark_mode_early():
    """Detect OS dark mode before QApplication is created.

    Uses platform-native APIs so the remote process can be started with
    the correct matplotlib style before Qt is set up.
    """
    try:
        if sys.platform == 'darwin':
            import subprocess
            r = subprocess.run(['defaults', 'read', '-g', 'AppleInterfaceStyle'],
                               capture_output=True, text=True, timeout=1)
            return r.returncode == 0 and r.stdout.strip().lower() == 'dark'
        elif sys.platform == 'win32':
            import winreg
            key = winreg.OpenKey(
                winreg.HKEY_CURRENT_USER,
                r'Software\Microsoft\Windows\CurrentVersion\Themes\Personalize')
            val, _ = winreg.QueryValueEx(key, 'AppsUseLightTheme')
            return val == 0
        else:  # Linux / other
            import subprocess
            r = subprocess.run(
                ['gsettings', 'get', 'org.gnome.desktop.interface', 'color-scheme'],
                capture_output=True, text=True, timeout=1)
            return 'dark' in r.stdout.lower()
    except Exception:
        pass
    return False


def _apply_theme(qapp):
    """
    Apply a light or dark Fusion palette matching the OS system theme, and
    set matplotlib's default style to match.  Called once after QApplication
    is configured.
    """
    from pyface.qt import QtGui, QtCore
    import matplotlib

    # Detect dark mode: check if the window background is darker than 50%.
    bg = qapp.palette().color(QtGui.QPalette.ColorRole.Window)
    is_dark = bg.lightness() < 128

    if is_dark:
        # Dark Fusion palette
        p = QtGui.QPalette()
        c = QtGui.QColor
        p.setColor(QtGui.QPalette.ColorRole.Window,          c(45,  45,  45))
        p.setColor(QtGui.QPalette.ColorRole.WindowText,      c(220, 220, 220))
        p.setColor(QtGui.QPalette.ColorRole.Base,            c(30,  30,  30))
        p.setColor(QtGui.QPalette.ColorRole.AlternateBase,   c(53,  53,  53))
        p.setColor(QtGui.QPalette.ColorRole.ToolTipBase,     c(255, 255, 220))
        p.setColor(QtGui.QPalette.ColorRole.ToolTipText,     c(0,   0,   0))
        p.setColor(QtGui.QPalette.ColorRole.Text,            c(220, 220, 220))
        p.setColor(QtGui.QPalette.ColorRole.Button,          c(53,  53,  53))
        p.setColor(QtGui.QPalette.ColorRole.ButtonText,      c(220, 220, 220))
        p.setColor(QtGui.QPalette.ColorRole.BrightText,      c(255,  80,  80))
        p.setColor(QtGui.QPalette.ColorRole.Link,            c(42, 130, 218))
        p.setColor(QtGui.QPalette.ColorRole.Highlight,       c(42, 130, 218))
        p.setColor(QtGui.QPalette.ColorRole.HighlightedText, c(0,   0,   0))
        # Disabled colours
        p.setColor(QtGui.QPalette.ColorGroup.Disabled,
                   QtGui.QPalette.ColorRole.ButtonText, c(127, 127, 127))
        p.setColor(QtGui.QPalette.ColorGroup.Disabled,
                   QtGui.QPalette.ColorRole.WindowText, c(127, 127, 127))
        p.setColor(QtGui.QPalette.ColorGroup.Disabled,
                   QtGui.QPalette.ColorRole.Text,       c(127, 127, 127))
        qapp.setPalette(p)

        # Matplotlib dark style (main process — for export_task etc.)
        matplotlib.rcParams.update(_MATPLOTLIB_DARK_PARAMS)
    else:
        # Keep matplotlib defaults (white background)
        pass


def run_gui():
    """Run the GUI!"""
    import os
    
    try:
        # if we're running as a one-click from a MacOS app,
        # we need to reset the working directory
        os.chdir(sys._MEIPASS)  # @UndefinedVariable
    except:
        # if we're not running as a one-click, fail gracefully
        pass
    
    # disable OpenGL in the qt web engine
    os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = "--disable-gpu"
    os.environ["QMLSCENE_DEVICE"] = "softwarecontext"

    # disable OpenGL on raster surfaces on linux / XCB
    os.environ["QT_XCB_GL_INTEGRATION"] = "none"

    # on macOS: opt into Retina/HiDPI rendering for the web engine
    if sys.platform == 'darwin':
        os.environ.setdefault("QT_ENABLE_HIGHDPI_SCALING", "1")

    from pyface.qt import QtGui, QtCore

    # take care of the 3 places in the cytoflow module that
    # need different behavior in a GUI
    import cytoflow
    cytoflow.RUNNING_IN_GUI = True
            
    # check that we're using the right Qt API
    from pyface.qt import qt_api

    cmd_line = " ".join(sys.argv)
    
    # parse args
    parser = argparse.ArgumentParser(description = 'Cytoflow GUI')
    parser.add_argument("--debug", action = 'store_true')
    parser.add_argument("filename", nargs='?', default = "")
    
    args = parser.parse_args()

    # detect dark mode before QApplication exists so the remote process
    # can be started with the right matplotlib style
    _early_dark = _detect_dark_mode_early()

    # start the remote process

    remote_process, remote_workflow_connection, remote_canvas_connection, queue_listener = start_remote_process(is_dark=_early_dark)
    
    # getting real tired of the matplotlib deprecation warnings
    import warnings
    warnings.filterwarnings('ignore', '.*is deprecated and replaced with.*')
    
    # if we're frozen, add _MEIPASS to the pyface search path for icons etc
    if getattr(sys, 'frozen', False):
        from pyface.resource_manager import resource_manager
        resource_manager.extra_paths.append(sys._MEIPASS)    # @UndefinedVariable
        
    # these three lines stop pkg_resources from trying to load resources
    # from the __main__ module, which is frozen (and thus not loadable.)
    from pyface.image_resource import ImageResource
    icon = ImageResource('icon')
    icon.search_path = []
    
    # monkey patch the resource manager to use SVGs for icons
    import pyface.resource.resource_manager
    pyface.resource.resource_manager.ResourceManager.IMAGE_EXTENSIONS.append('.svg')
    
    # monkey patch checklist editor to stop lowercasing
    import traitsui.qt.check_list_editor  # @UnusedImport
    traitsui.qt.check_list_editor.capitalize = lambda s: s
    
    # monkey patch ApplicationWindow to fix toolbar bug
    from traits.api import observe
    @observe("tool_bar_managers.items")
    def _update_tool_bar_managers(self, event):
        if self.control is not None:
            # Remove the old toolbars.
            for child in self.control.children():
                if isinstance(child, QtGui.QToolBar):
                    self.control.removeToolBar(child)
                    child.deleteLater()

            # Add the new toolbars.
            if event.new is not None:
                self._create_tool_bar(self.control)
                
    from pyface.ui.qt.application_window import ApplicationWindow
    ApplicationWindow._update_tool_bar_managers = _update_tool_bar_managers

    if sys.platform == 'darwin':
        # Use the native macOS unified title-and-toolbar so the toolbar merges
        # into the title bar, saving vertical space on a laptop.
        _orig_create = ApplicationWindow._create_control
        def _mac_create_control(self, parent):
            ctrl = _orig_create(self, parent)
            try:
                ctrl.setUnifiedTitleAndToolBarOnMac(True)
            except Exception:
                pass
            return ctrl
        ApplicationWindow._create_control = _mac_create_control

    # -------------------------------------------------------------------------
    # Drag-and-drop FCS files onto the main window
    # -------------------------------------------------------------------------

    def _handle_fcs_drop(window, paths):
        """Add dropped FCS files to the import operation (creates one if needed).
        Deferred to the next event-loop tick so the window is fully set up."""
        from pyface.qt import QtCore
        QtCore.QTimer.singleShot(0, lambda: _do_fcs_drop(window, paths))

    def _do_fcs_drop(window, paths):
        app = window.application   # CytoflowApplication
        if app is None:
            return

        workflow = getattr(app, 'model', None)
        if workflow is None:
            return

        from cytoflowgui.workflow.operations.import_op import ImportWorkflowOp
        from cytoflow.operations.import_op import Tube
        from cytoflowgui.workflow.workflow_item import WorkflowItem

        # find existing ImportWorkflowOp or create one
        import_op = None
        for item in workflow.workflow:
            if isinstance(item.operation, ImportWorkflowOp):
                import_op = item.operation
                break

        if import_op is None:
            import_op = ImportWorkflowOp()
            for path in paths:
                import_op.tubes.append(Tube(file=path, conditions={}))
            wi = WorkflowItem(operation=import_op)
            workflow.workflow.append(wi)
        else:
            existing = {t.file for t in import_op.tubes}
            for path in paths:
                if path not in existing:
                    import_op.tubes.append(Tube(file=path, conditions={}))
                    existing.add(path)

    # Patch ApplicationWindow._create_control to install drop support.
    # We always do this (all platforms), chaining from whatever is currently set.
    _orig_create_for_drop = ApplicationWindow._create_control

    def _drop_create_control(self, parent):
        ctrl = _orig_create_for_drop(self, parent)
        ctrl.setAcceptDrops(True)

        def dragEnterEvent(event):
            urls = event.mimeData().urls() if event.mimeData().hasUrls() else []
            if any(u.toLocalFile().lower().endswith('.fcs') for u in urls):
                event.acceptProposedAction()
            else:
                event.ignore()

        def dropEvent(event):
            paths = [u.toLocalFile() for u in event.mimeData().urls()
                     if u.toLocalFile().lower().endswith('.fcs')]
            if paths:
                _handle_fcs_drop(self, paths)
            event.acceptProposedAction()

        ctrl.dragEnterEvent = dragEnterEvent
        ctrl.dropEvent = dropEvent
        return ctrl

    ApplicationWindow._create_control = _drop_create_control

    # macOS: also handle "Open With" / Dock-drag via QFileOpenEvent
    if sys.platform == 'darwin':
        from pyface.qt import QtGui as _QtGui, QtCore as _QtCore

        _qapp_inst = _QtGui.QApplication.instance()
        if _qapp_inst is not None:
            _orig_qapp_event = _qapp_inst.__class__.event

            def _qapp_event(self, e):
                if e.type() == _QtCore.QEvent.FileOpen:
                    path = e.file()
                    if path.lower().endswith('.fcs'):
                        # The window may not exist yet; use a generous delay.
                        def _try_drop(p=path):
                            from pyface.qt import QtGui as _Qg
                            app_obj = _Qg.QApplication.instance()
                            # Walk top-level widgets looking for an ApplicationWindow
                            for widget in app_obj.topLevelWidgets():
                                win = getattr(widget, '_pyface_window', None)
                                if win is not None:
                                    _handle_fcs_drop(win, [p])
                                    return
                        _QtCore.QTimer.singleShot(500, _try_drop)
                return _orig_qapp_event(self, e)

            _qapp_inst.__class__.event = _qapp_event

    # -------------------------------------------------------------------------
    # End drag-and-drop setup
    # -------------------------------------------------------------------------

    # monkey patch pyface.image.image.ImageVolume to stop overwriting the
    # standard image libraries, which breaks Mac code signing
    
    def empty_save(self):
        pass
    
    from pyface.image.image import ImageVolume
    ImageVolume.save = empty_save
    
    # define and install a message handler for Qt errors
    from traits.api import push_exception_handler  # @UnresolvedImport
                             
    def QtMsgHandler(msg_type, msg_context, msg_string):
        # Convert Qt msg type to logging level
        log_level = {QtCore.QtMsgType.QtDebugMsg : logging.DEBUG,
                     QtCore.QtMsgType.QtInfoMsg : logging.INFO,
                     QtCore.QtMsgType.QtWarningMsg : logging.WARN,
                     QtCore.QtMsgType.QtCriticalMsg : logging.ERROR,
                     QtCore.QtMsgType.QtFatalMsg : logging.FATAL}
        logging.log(log_level[msg_type], 'Qt message: ' + msg_string)
        
    from pyface.qt.QtCore import qInstallMessageHandler  # @UnresolvedImport
    qInstallMessageHandler(QtMsgHandler)
    
    # install a global (gui) error handler for traits notifications
    push_exception_handler(handler = log_notification_handler,
                           reraise_exceptions = False, 
                           main = True)
    
    sys.excepthook = log_excepthook
    
    # Import, then load, the envisage plugins
    
    from envisage.core_plugin import CorePlugin
    from envisage.ui.tasks.tasks_plugin import TasksPlugin
    
    from cytoflowgui.flow_task import FlowTaskPlugin
    from cytoflowgui.export_task import ExportFigurePlugin
    from cytoflowgui.cytoflow_application import CytoflowApplication

    from cytoflowgui.op_plugins import (OpPluginManager,
                                        ImportPlugin,
                                        SampleLabelsPlugin,
                                        ThresholdPlugin,
                                        RangePlugin, 
                                        QuadPlugin,
                                        Range2DPlugin, 
                                        PolygonPlugin, 
                                        EllipsePlugin,
                                        FreehandPlugin,
                                        HierarchyPlugin,
                                        CategoryPlugin,
                                        BooleanGatePlugin,
                                        BinningPlugin,
                                        GaussianMixture1DPlugin, 
                                        GaussianMixture2DPlugin,
                                        DensityGatePlugin, 
                                        BleedthroughLinearPlugin,
                                        CompensationPlugin,
                                        AutoSpillPlugin,
                                        CrossEntropyTestPlugin,
                                        BeadCalibrationPlugin, 
                                        AutofluorescencePlugin,
                                        ColorTranslationPlugin, 
                                        TasbePlugin, 
                                        ChannelStatisticPlugin,
                                        MultiChannelStatisticPlugin,
                                        TransformStatisticPlugin,
                                        RatioPlugin,
                                        ScaleSettingsPlugin,
                                        StatisticalTestPlugin,
                                        DerivedParameterPlugin,
                                        FlowPeaksPlugin,
                                        KMeansPlugin,
                                        PhenographPlugin,
                                        PopulationTreePlugin,
                                        ProliferationPlugin,
                                        CellCyclePlugin,
                                        PCAPlugin,
                                        FlowCleanPlugin, 
                                        tSNEPlugin,
                                        UMAPPlugin,
                                        SOMPlugin,
                                        MSTPlugin as MSTOpPlugin,
                                        RegistrationPlugin)
    
    from cytoflowgui.view_plugins import (ViewPluginManager,
                                          HistogramPlugin, 
                                          Histogram2DPlugin, 
                                          ScatterplotPlugin,
                                          BarChartPlugin,
                                          Stats1DPlugin, 
                                          Stats2DPlugin,
                                          MatrixPlugin,
                                          MSTPlugin as MSTViewPlugin,
                                          Kde1DPlugin, 
                                          Kde2DPlugin,
                                          ViolinPlotPlugin,
                                          TablePlugin, 
                                          LongTablePlugin,
                                          DensityPlugin,
                                          ParallelCoordinatesPlugin,
                                          RadvizPlugin,
                                          ExportFCSPlugin)

    plugins = [CorePlugin(), TasksPlugin(), FlowTaskPlugin(), ExportFigurePlugin()]

    # ordered as we want them to show up in the toolbar    
    view_plugins = [ViewPluginManager(),
                    HistogramPlugin(),
                    ScatterplotPlugin(),
                    Histogram2DPlugin(),
                    DensityPlugin(),
                    Kde1DPlugin(),
                    Kde2DPlugin(),
                    RadvizPlugin(),
                    ParallelCoordinatesPlugin(),
                    ViolinPlotPlugin(),
                    BarChartPlugin(),
                    Stats1DPlugin(),
                    Stats2DPlugin(),
                    MatrixPlugin(),
                    MSTViewPlugin(),
                    TablePlugin(),
                    LongTablePlugin(),
                    ExportFCSPlugin()]
    
    plugins.extend(view_plugins)
    
    op_plugins = [OpPluginManager(),
                  ImportPlugin(),
                  SampleLabelsPlugin(),
                  ThresholdPlugin(),
                  RangePlugin(),
                  QuadPlugin(),
                  Range2DPlugin(),
                  PolygonPlugin(),
                  EllipsePlugin(),
                  FreehandPlugin(),
                  HierarchyPlugin(),
                  CategoryPlugin(),
                  BooleanGatePlugin(),
                  ChannelStatisticPlugin(),
                  MultiChannelStatisticPlugin(),
                  TransformStatisticPlugin(),
                  #MergeStatisticsPlugin(),
                  RatioPlugin(),
                  ScaleSettingsPlugin(),
                  StatisticalTestPlugin(),
                  DerivedParameterPlugin(),
                  FlowCleanPlugin(),
                  RegistrationPlugin(),
                  BinningPlugin(),
                  GaussianMixture1DPlugin(),
                  GaussianMixture2DPlugin(),
                  DensityGatePlugin(),
                  KMeansPlugin(),
                  PhenographPlugin(),
                  PopulationTreePlugin(),
                  ProliferationPlugin(),
                  CellCyclePlugin(),
                  SOMPlugin(),
                  FlowPeaksPlugin(),
                  PCAPlugin(),
                  tSNEPlugin(),
                  UMAPPlugin(),
                  MSTOpPlugin(),
                  CompensationPlugin(),
                  AutoSpillPlugin(),
                  CrossEntropyTestPlugin(),
                  AutofluorescencePlugin(),
                  BleedthroughLinearPlugin(),
                  BeadCalibrationPlugin(),
                  ColorTranslationPlugin(),
                  TasbePlugin()]

    plugins.extend(op_plugins)
    
    # start the app

    app = CytoflowApplication(id = 'cytoflow',
                              plugins = plugins,
                              icon = icon,
                              remote_process = remote_process,
                              remote_workflow_connection = remote_workflow_connection,
                              remote_canvas_connection = remote_canvas_connection,
                              filename = args.filename,
                              debug = args.debug)

    qapp = QtGui.QApplication.instance()
    qapp.setStyle(QtGui.QStyleFactory.create('Fusion'))

    # Detect system dark mode (works on macOS and Windows 10+).
    # Qt6 / PyQt5 may auto-detect; we do it explicitly so matplotlib also
    # follows the same palette.
    _apply_theme(qapp)

    app.run()

    remote_process.join()
    queue_listener.stop()
    logging.shutdown()
    

        
def monitor_remote_process(proc):
    """The main method for the (local) thread that monitors the remote process"""
    proc.join()
    if proc.exitcode:
        logging.error("Remote process exited with {}".format(proc.exitcode))
    

def start_remote_process(is_dark=False):
    """
    Start the remote process.  Creates pipes and synchronization primitives,
    sets up logging, and starts the remote process and the monitoring thread.
    """

    # communications channels

    parent_workflow_conn, child_workflow_conn = multiprocessing.Pipe()
    parent_mpl_conn, child_matplotlib_conn = multiprocessing.Pipe()
    running_event = multiprocessing.Event()

    # logging
    from logging.handlers import QueueListener
    from cytoflowgui.utility import CallbackHandler
    log_q = multiprocessing.Queue()
    def handle(record):
        logger = logging.getLogger(record.name)
        if logger.isEnabledFor(record.levelno):
            logger.handle(record)

    handler = CallbackHandler(handle)
    queue_listener = QueueListener(log_q, handler)
    queue_listener.start()

    remote_process = multiprocessing.Process(target = remote_main,
                                             name = "remote process",
                                             args = [parent_workflow_conn,
                                                     parent_mpl_conn,
                                                     log_q,
                                                     running_event,
                                                     is_dark])
    
    remote_process.daemon = True
    remote_process.start() 
    running_event.wait()
    
    remote_process_thread = threading.Thread(target = monitor_remote_process,
                                             name = "monitor remote process",
                                             args = [remote_process])
    remote_process_thread.daemon = True
    remote_process_thread.start()
    
    return (remote_process, child_workflow_conn, child_matplotlib_conn, queue_listener)
 

def remote_main(parent_workflow_conn, parent_mpl_conn, log_q, running_event,
                is_dark=False):
    """
    The main method for the remote process. Configures logging, sets up
    the matplotlib backend, and instantiates the `RemoteWorkflow`.
    """
    
    # this should only ever be main method after a spawn() call 
    # (not fork). So we should have a fresh logger to set up.
        
    # messages that end up at the root logger to go (only to) log_q
    logging.getLogger().handlers.clear()
    from logging.handlers import QueueHandler
    h = QueueHandler(log_q) 
    logging.getLogger().addHandler(h)
    
    # make sure the root logger has a level of DEBUG -- we'll sort out what
    # to show or not on the local logger
    logging.getLogger().setLevel(logging.DEBUG)
    
    # capture everything sent to stdout in the log too
    class StreamToLogger(object):
        """
        Fake file-like stream object that redirects writes to a logger instance.
        """
        def __init__(self, logger, level):
            self.logger = logger
            self.level = level
            self.linebuf = ''
    
        def write(self, buf):
            for line in buf.rstrip().splitlines():
                self.logger.log(self.level, line.rstrip())
    
        def flush(self):
            pass
        
    sys.stdout = StreamToLogger(logging.getLogger(),logging.INFO)
    
    # We want matplotlib to use our backend .... in both the GUI and the
    # remote process.  Must be called BEFORE cytoflow is imported
    
    import matplotlib
    matplotlib.use('module://cytoflowgui.matplotlib_backend_remote')

    if is_dark:
        matplotlib.rcParams.update(_MATPLOTLIB_DARK_PARAMS)

    # install a global (gui) error handler for traits notifications
        
    from traits.api import push_exception_handler        # @UnresolvedImport
    push_exception_handler(handler = log_notification_handler,
                           reraise_exceptions = False,
                           main = True)
    
    sys.excepthook = log_excepthook
    
    # take care of the 3 places in the cytoflow module that
    # need different behavior in a GUI
    import cytoflow
    cytoflow.RUNNING_IN_GUI = True
    
    running_event.set()
    
    from cytoflowgui.workflow import RemoteWorkflow
    RemoteWorkflow().run(parent_workflow_conn, parent_mpl_conn)

if __name__ == '__main__':
    multiprocessing.freeze_support()
    multiprocessing.set_start_method('spawn', force = True)
    run_gui()
