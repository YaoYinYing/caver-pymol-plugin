"""
Advanced UI for Caver, originally written by Yinying for REvoDesign Project.
"""


import math
import os
import time
import warnings
from collections.abc import Iterable
from contextlib import contextmanager
from typing import TYPE_CHECKING, Any, Callable, NoReturn, Optional, TypeVar, Union, overload

from ..caver_pymol import ROOT_LOGGER

if TYPE_CHECKING:
    from PyQt5 import QtCore, QtGui, QtWidgets
else:
    from pymol.Qt import QtCore, QtGui, QtWidgets

logging=ROOT_LOGGER.getChild('UI_Tape')

@overload
def set_widget_value(widget: QtWidgets.QStackedWidget, value: list): ...


@overload
def set_widget_value(widget: Union[QtWidgets.QProgressBar, QtWidgets.QSlider], value: Union[int, list[int], tuple[int, int]]): ...


@overload
def set_widget_value(
    widget: Union[QtWidgets.QDoubleSpinBox, QtWidgets.QSpinBox], value: Union[int, float, list[str], tuple[str, str]]
): ...


@overload
def set_widget_value(widget: QtWidgets.QComboBox, value: Union[list, tuple, dict, str, int, float, bool]): ...


@overload
def set_widget_value(widget: QtWidgets.QGridLayout, value: str): ...


@overload
def set_widget_value(
    widget: Union[
        QtWidgets.QLineEdit, QtWidgets.QTextEdit, QtWidgets.QLCDNumber, QtWidgets.QCheckBox, QtWidgets.QLabel
    ],
    value: Any,
): ...


def set_widget_value(widget, value):
    """
    Sets the value of a PyQt5 widget based on the provided value.

    Args:
    - widget: The PyQt5 widget whose value needs to be set.
    - value: The value to be set on the widget.

    Supported Widgets and Value Types:
    - QDoubleSpinBox: Supports int, float, list or tuple (for setting range).
    - QSpinBox: Supports int, float, list or tuple (for setting range).
    - QComboBox: Supports str, list, tuple, dict.
    - QLineEdit: Supports str.
    - QLabel: Supports str.
    - QProgressBar: Supports int, list or tuple (for setting range).
    - QLCDNumber: Supports any value (converted to string for display).
    - QCheckBox: Supports bool.
    - QStackedWidget: Supports list of image paths (adds ImageWidget widgets).
    - QGridLayout: Supports a string (image path) to add an ImageWidget widget.
    """

    def set_value_error(widget: QtWidgets.QWidget, value: Any):
        logging.warning(f"FIX ME: Value {value} is not currently supported on widget {type(widget).__name__}")

    # Preprocess values according to types
    if callable(value):
        value = value()  # Call the function to get the value if value is callable

    if isinstance(value, Iterable) and not isinstance(value, (str, list, tuple, dict)):
        value = list(value)  # Convert iterable (excluding strings, lists, tuples, dicts) to list

    # Setting values
    if isinstance(widget, QtWidgets.QDoubleSpinBox):
        if isinstance(value, (int, float)):
            widget.setValue(float(value))
        elif isinstance(value, (list, tuple)) and len(value) > 1:
            widget.setRange(float(value[0]), float(value[1]))
    elif isinstance(widget, QtWidgets.QSpinBox):
        if isinstance(value, (int, float)):
            widget.setValue(int(value))
        elif isinstance(value, (list, tuple)) and len(value) > 1:
            widget.setRange(int(value[0]), int(value[1]))

    elif isinstance(widget, QtWidgets.QComboBox):
        if isinstance(value, (list, tuple)):
            widget.clear()
            widget.addItems(map(str, value))
        elif isinstance(value, dict):
            widget.clear()
            for k, v in value.items():
                widget.addItem(v, k)
        else:
            widget.setCurrentText(str(value))
    elif isinstance(widget, QtWidgets.QLineEdit) or isinstance(widget, QtWidgets.QLabel):
        widget.setText(str(value))
    elif isinstance(widget, QtWidgets.QTextEdit):
        widget.setPlainText(str(value))
    elif isinstance(widget, QtWidgets.QProgressBar, QtWidgets.QSlider):
        if isinstance(value, int):
            widget.setValue(value)
        elif isinstance(value, (list, tuple)) and len(value) == 2:
            widget.setRange(*value)
    elif isinstance(widget, QtWidgets.QLCDNumber):
        widget.display(str(value))
    elif isinstance(widget, QtWidgets.QCheckBox):
        widget.setChecked(bool(value))

    else:
        set_value_error(widget, value)


@overload
def get_widget_value(widget: QtWidgets.QCheckBox) -> bool: ...  # type: ignore


@overload
def get_widget_value(
    widget: Union[QtWidgets.QComboBox, QtWidgets.QLineEdit, QtWidgets.QTextEdit],  # type: ignore
) -> str: ...


@overload
def get_widget_value(widget: Union[QtWidgets.QDoubleSpinBox, QtWidgets.QLCDNumber]) -> float: ...  # type: ignore


@overload
def get_widget_value(widget: Union[QtWidgets.QSpinBox, QtWidgets.QProgressBar, QtWidgets.QSlider]) -> int: ...  # type: ignore


def get_widget_value(widget: QtWidgets.QWidget) -> Any:
    """
    Retrieves the value from a PyQt5 widget.

    Args:
    - widget: The PyQt5 widget from which the value needs to be retrieved.

    Returns:
    The current value of the widget.

    Supported Widgets:
    - QDoubleSpinBox, QSpinBox: Returns the current value as float or int.
    - QComboBox: Returns the current text or the userData of the current item if any.
    - QLineEdit: Returns the current text as str.
    - QProgressBar: Returns the current value as int.
    - QLCDNumber: Returns the current value as float.
    - QCheckBox: Returns the checked state as bool.

    Raises:
    - ValueError: If the widget type is not supported for value retrieval.
    """
    if isinstance(widget, QtWidgets.QDoubleSpinBox) or isinstance(widget, QtWidgets.QSpinBox):
        return widget.value()
    if isinstance(widget, QtWidgets.QComboBox):
        return widget.currentText()
    if isinstance(widget, QtWidgets.QLineEdit):
        return widget.text()
    if isinstance(widget, QtWidgets.QProgressBar, QtWidgets.QSlider):
        return widget.value()
    if isinstance(widget, QtWidgets.QLCDNumber):
        return float(widget.value())
    if isinstance(widget, QtWidgets.QCheckBox):
        return widget.isChecked()
    if isinstance(widget, QtWidgets.QTextEdit):
        return widget.toPlainText()

    raise ValueError(f"Widget type {type(widget).__name__} is not supported for value retrieval.")


def widget_signal_tape(widget: QtWidgets.QWidget, event):
    """
    Connects the appropriate signal of a QWidget to the specified event handler.

    This function connects specific signals from different types of QWidgets to a unified event handler.
    It handles several common Qt widget types such as QDoubleSpinBox, QSpinBox, QProgressBar,
    QComboBox, QLineEdit, and QCheckBox, binding their respective signals to the provided event handler.

    Parameters:
    - widget (QtWidgets.QWidget): The widget instance whose signal will be connected.
    - event (callable): The event handler function that will be called when the widget's signal is emitted.

    Raises:
    - NotImplementedError: If the widget type is not supported by this function.
    """

    # Handle numeric input widgets and progress bar
    if isinstance(
        widget,
        (
            QtWidgets.QDoubleSpinBox,
            QtWidgets.QSpinBox,
            QtWidgets.QProgressBar,
        ),
    ):
        widget.valueChanged.connect(event)

    # Handle combo box widgets with text change signals
    elif isinstance(widget, QtWidgets.QComboBox):
        widget.currentTextChanged.connect(event)
        widget.editTextChanged.connect(event)

    # Handle line edit widgets with text change signals
    elif isinstance(widget, QtWidgets.QLineEdit):
        widget.textChanged.connect(event)
        widget.textEdited.connect(event)

    # Handle checkbox widgets with state change signals
    elif isinstance(widget, QtWidgets.QCheckBox):
        widget.stateChanged.connect(event)

    # Raise an error for unsupported widget types
    else:
        raise NotImplementedError(f"{widget} {type(widget)} is not supported yet")


def getExistingDirectory():
    return QtWidgets.QFileDialog.getExistingDirectory(  # type: ignore
        None,
        "Open Directory",
        os.path.expanduser("~"),
        QtWidgets.QFileDialog.ShowDirsOnly | QtWidgets.QFileDialog.DontResolveSymlinks,  # type: ignore
    )


def refresh_window():
    """
    Refresh the application window by processing all pending events.
    This function is copied from `REvoDesign/tools/customized_widgets.py`.

    No parameters are required for this function.

    Returns:
        None
    """
    QtWidgets.QApplication.processEvents()


@overload
def notify_box(
    message: str = "", error_type: Union[None, type[Warning]] = None, details: Optional[str] = None
) -> None: ...


# Overload #2: Exception => NoReturn


@overload
def notify_box(message: str, error_type: type[Exception], details: Optional[str] = None) -> NoReturn: ...


def notify_box(
    message: str = "", error_type: Optional[type[Exception]] = None, details: Optional[str] = None
) -> Union[None, NoReturn]:
    """
    Display a notification message box.

    Parameters:
    - message: str, the content of the message box.
    - error_type: Optional[Union[Exception, Warning]], the type of error or warning, can be None.

    If `error_type` is None or a Warning, returns bool.
    If `error_type` is an Exception (and not a Warning), raises => NoReturn.
    """
    refresh_window()
    # Create an information message box
    msg = QtWidgets.QMessageBox()

    if error_type is None:
        msg.setIcon(QtWidgets.QMessageBox.Information)
    elif issubclass(error_type, Warning):
        msg.setIcon(QtWidgets.QMessageBox.Warning)
    elif issubclass(error_type, Exception):
        msg.setIcon(QtWidgets.QMessageBox.Critical)

    msg.setText(message)
    if details is not None:
        msg.setDetailedText(details)

    msg.setStandardButtons(QtWidgets.QMessageBox.Ok)
    # Display the message box
    msg.exec_()
    # If error_type is None, end the function execution
    if error_type is None:
        return

    # error_type is a Warning => also bool
    if issubclass(error_type, Warning):
        warnings.warn(error_type(message))
        return

    # Otherwise, raise => NoReturn
    raise_error(error_type, message)


def raise_error(error_type: type[Exception], message: str) -> NoReturn:
    """
    Raises an error of the specified type with the given message.

    Args:
    - error_type: Type[Exception], the type of error to raise.
    - message: str, the error message.
    """
    raise error_type(message)


class CheckableListView(QtWidgets.QWidget):
    """
    Checkable list view widget, allowing users to check items in the list.

    Attributes:
        list_view: The QListView instance this widget operates on.
        model: The data model instance used by the list view.
    """

    # Custom signal, emits (item_text: str, new_state: Qt.CheckState)
    checkStateChanged = QtCore.pyqtSignal(list)

    def __init__(self, list_view, items: dict[str, str] = {}, parent=None):
        """
        Initializes the CheckableListView instance.

        Parameters:
            listView: The QListView instance to use.
            items: Optional list of item texts to add to the list.
            separators: Optional list of separator texts, used to categorize items.
            parent: The parent widget, defaults to None.
        """
        super().__init__(parent)

        # Use the existing list view
        self.list_view = list_view

        # Set up the model (use existing one if set, otherwise create a new one)
        if self.list_view.model() is None:
            self.model = QtGui.QStandardItemModel(self.list_view)
            self.list_view.setModel(self.model)
        else:
            self.model = self.list_view.model()

        # Connect to model's itemChanged signal to detect check state changes
        self.model.itemChanged.connect(self._on_item_changed)

        # Clear the model before adding new items
        self.model.clear()

        # Add items to the model with optional separators
        if not items:
            return

        self.items = {}

        self.update(items)

    def _on_item_changed(self):
        """
        Internal slot: triggered when an item's state changes.
        Emits the checkStateChanged signal only for checkable items.
        """
        self.checkStateChanged.emit(self.get_checked_items())

    def _get_items_by_check_state(self, check_state):
        """
        Helper function to get items based on their check state.

        Args:
            check_state (int): The check state to filter items by (e.g., QtCore.Qt.Checked).

        Returns:
            A list of strings representing the texts of items with the specified check state.
        """
        items = []
        for row in range(self.model.rowCount()):
            item = self.model.item(row)
            if item.isCheckable() and item.checkState() == check_state:
                items.append(self.items.get(item.text(), None))
        return items

    def get_checked_items(self):
        """
        Returns a list of all checked items' text.

        Returns:
            A list of strings representing the texts of all checked items.
        """
        checked_items = self._get_items_by_check_state(QtCore.Qt.Checked)
        logging.debug(f"Checked: {checked_items}")
        return checked_items

    def get_unchecked_items(self):
        """
        Returns a list of all unchecked items' text.

        Returns:
            A list of strings representing the texts of all unchecked items.
        """
        return self._get_items_by_check_state(QtCore.Qt.Unchecked)

    def check_all(self):
        """
        Check all items in the list, excluding separators.
        """
        for row in range(self.model.rowCount()):
            item = self.model.item(row)
            if item.isCheckable():
                item.setCheckState(QtCore.Qt.Checked)

    def uncheck_all(self):
        """
        Uncheck all items in the list, excluding separators.
        """
        for row in range(self.model.rowCount()):
            item = self.model.item(row)
            if item.isCheckable():
                item.setCheckState(QtCore.Qt.Unchecked)

    def reverse_check(self):
        """
        Reverse the check state of all items in the list, excluding separators.
        """
        for row in range(self.model.rowCount()):
            item = self.model.item(row)
            if item.isCheckable():
                if item.checkState() == QtCore.Qt.Checked:
                    item.setCheckState(QtCore.Qt.Unchecked)
                else:
                    item.setCheckState(QtCore.Qt.Checked)

    def check_these(self, required_items: list[str], clear_before_check: bool = True):
        """
        Selects the items in the list.
        """
        if not required_items:
            return
        if clear_before_check:
            self.uncheck_all()

        for row in range(self.model.rowCount()):
            item = self.model.item(row)
            if not (item.isCheckable() and item.text() in required_items):
                continue
            item.setCheckState(QtCore.Qt.Checked)

    def update(self, new_items: dict):
        for k, v in new_items.items():
            if not v:
                # Add as a separator
                separator_item = QtGui.QStandardItem(k)
                separator_item.setEnabled(False)  # Non-interactive
                separator_item.setSelectable(False)  # Non-selectable
                separator_item.setCheckable(False)  # Non-checkable
                separator_item.setForeground(QtGui.QBrush(QtCore.Qt.yellow))
                separator_item.setBackground(QtGui.QBrush(QtCore.Qt.blue))  # Different background
                separator_item.setFont(QtGui.QFont("Arial", weight=QtGui.QFont.Bold))  # Bold text
                self.model.appendRow(separator_item)
            else:
                # Add as a regular checkable item
                item = QtGui.QStandardItem(k)
                item.setCheckable(True)
                item.setCheckState(QtCore.Qt.Unchecked)  # Default unchecked
                self.model.appendRow(item)

        self.items.update(new_items)


def getOpenFileNameWithExt(*args, **kwargs):
    """
    Return a file name, append extension from filter if no extension provided.
    """
    import re

    fname, filter = QtWidgets.QFileDialog.getOpenFileName(*args, **kwargs)  # type: ignore

    if not fname:
        return ""

    if "." not in os.path.split(fname)[-1]:
        m = re.search(r"\*(\.[\w\.]+)", filter)
        if m:
            # append first extension from filter
            fname += m.group(1)

    return fname


@contextmanager
def hold_trigger_button(
    buttons: Union[tuple[QtWidgets.QPushButton, ...], QtWidgets.QPushButton],
    animation_duration: int = 1000,  # Duration of the breathing cycle (in milliseconds)
):
    """
    A context manager for holding and releasing trigger buttons with a breathing effect
    using the system's accent color.

    Args:
        buttons: One or more QPushButton objects.
        animation_duration: Duration of the breathing animation cycle (in milliseconds).
    """
    if not isinstance(buttons, (tuple, list, set)):
        buttons = (buttons,)

    timers = []

    def get_accent_color():
        color = QtGui.QColor(76, 217, 100)
        return color

    def start_breathing_animation(button: QtWidgets.QPushButton):
        accent_color = get_accent_color()
        base_color = accent_color.lighter(150)  # Start with a lighter shade
        darker_color = accent_color.darker(150)  # Use a darker shade for the trough

        timer = QtCore.QTimer(button)
        timer.setInterval(30)  # Update every 30 milliseconds
        elapsed = 0

        def update_stylesheet():
            nonlocal elapsed
            elapsed += timer.interval()
            t = (elapsed % animation_duration) / animation_duration  # Normalized time [0, 1]
            # Calculate intermediate intensity using sine wave
            factor = (1 + math.sin(2 * math.pi * t)) / 2  # Normalized to [0, 1]
            r = int(base_color.red() * factor + darker_color.red() * (1 - factor))
            g = int(base_color.green() * factor + darker_color.green() * (1 - factor))
            b = int(base_color.blue() * factor + darker_color.blue() * (1 - factor))
            button.setStyleSheet(f"background-color: rgb({r}, {g}, {b});")

        timer.timeout.connect(update_stylesheet)
        timer.start()
        timers.append(timer)

    def stop_breathing_animation(button: QtWidgets.QPushButton):
        # Stop all timers associated with this button
        for timer in timers:
            if timer.parent() == button:
                timer.stop()
                timers.remove(timer)
        button.setStyleSheet("")  # Reset the button's style

    try:
        for b in buttons:
            b.setEnabled(False)
            b.setProperty("held", True)  # Mark the button as held
            b.setProperty("original_style", b.styleSheet() if b.styleSheet() else "")
            start_breathing_animation(b)
            logging.debug(f"Held button: {b.text()}: ({b.objectName()})")
        yield
    finally:
        for b in buttons:
            b.setProperty("held", False)  # Remove the held mark
            stop_breathing_animation(b)
            b.setStyleSheet(b.property("original_style") if b.property("original_style") else "")
            b.setEnabled(True)  # Re-enable the button
            logging.debug(f"Released button: {b.text()}: ({b.objectName()})")


R = TypeVar("R")


class WorkerThread(QtCore.QThread):
    """
    Custom worker thread for executing a function in a separate thread.

    Attributes:
    - result_signal (QtCore.pyqtSignal): Signal emitted when the result is available.
    - finished_signal (QtCore.pyqtSignal): Signal emitted when the thread finishes its execution.
    - interrupt_signal (QtCore.pyqtSignal): Signal to interrupt the thread.

    Methods:
    - __init__: Initializes the WorkerThread object.
    - run: Executes the specified function with arguments and emits the result through signals.
    - handle_result: Returns the result obtained after the thread execution.
    - interrupt: Interrupts the execution of the thread.

    Example Usage:
    ```python
    def some_function(x, y):
        return x + y

    worker = WorkerThread(func=some_function, args=(10, 20))
    worker.result_signal.connect(handle_result_function)
    worker.finished_signal.connect(handle_finished_function)
    worker.interrupt_signal.connect(handle_interrupt_function)
    worker.start()
    # To interrupt the execution:
    # worker.interrupt()
    ```
    """

    result_signal = QtCore.pyqtSignal(list)
    finished_signal = QtCore.pyqtSignal()
    interrupt_signal = QtCore.pyqtSignal()

    def __init__(self, func, args=None, kwargs=None):
        super().__init__()
        self.func = func
        self.args = args or ()
        self.kwargs = kwargs or {}
        self.results = None  # Define the results attribute

    def run(self):
        """
        Executes the task and handles the results.

        This function checks if an interruption has been requested. If not, it runs the specified function with
        given arguments and keyword arguments.
        The result is then emitted through a signal if available, and a completion signal is emitted at the end.

        Parameters:
        - self: The instance of the class containing this method. It should have the following attributes:
            - func: The function to be executed.
            - args: A tuple of positional arguments for the function.
            - kwargs: A dictionary of keyword arguments for the function.
            - result_signal: A signal to emit the results.
            - finished_signal: A signal to indicate the task has finished.
            - isInterruptionRequested: A method that returns True if an interruption has been requested,
            otherwise False.
        """
        # Check if an interruption has been requested
        if not self.isInterruptionRequested():
            # Execute the function with provided arguments and store the result
            self.results = [self.func(*self.args, **self.kwargs)]

            # Emit the result if it exists
            if self.results:
                self.result_signal.emit(self.results)

            # Emit the finished signal
            self.finished_signal.emit()

    def handle_result(self):
        """
        Retrieves the results from the current instance.

        This method returns the 'results' attribute of the current instance.
        It is used to obtain the result data within other methods of the class.
        """
        return self.results

    def interrupt(self):
        """
        Emit an interrupt signal.

        This function triggers an interrupt signal.
        """
        self.interrupt_signal.emit()


@overload
def run_worker_thread_with_progress(worker_function: Callable[..., R], *args, **kwargs) -> R: ...


def run_worker_thread_with_progress(worker_function: Callable[..., Optional[R]], *args, **kwargs) -> Optional[R]:
    """
    Runs a worker function in a separate thread and optionally updates a progress bar.

    This function is designed to execute a given task (worker_function) in a separate thread,
    allowing the main thread to remain responsive, such as updating a progress bar.
    After the task is completed, it restores the progress bar's state and returns the result of the task.

    Parameters:
    - worker_function: The function to execute in a separate thread.
    - *args, **kwargs: Additional arguments and keyword arguments to pass to the worker function.

    Returns:
    - The result of the worker function or None if no result is available.
    """

    # Create and start a worker thread with the given function and parameters
    work_thread = WorkerThread(worker_function, args=args, kwargs=kwargs)
    work_thread.start()

    # Keep the main thread running until the worker thread finishes
    while not work_thread.isFinished():
        refresh_window()
        time.sleep(0.01)

    # Obtain and return the result of the worker function
    result = work_thread.handle_result()

    return result[0] if result else None
