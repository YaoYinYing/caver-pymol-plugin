import logging
import os
from typing import TYPE_CHECKING, Any, Dict, Iterable, List, NoReturn, Optional, Type, Union, overload
import warnings
if TYPE_CHECKING:
    from PyQt5 import QtCore, QtGui, QtWidgets
else:
    from pymol.Qt import QtCore, QtGui, QtWidgets



@overload
def set_widget_value(widget: QtWidgets.QStackedWidget, value: list): ...


@overload
def set_widget_value(widget: QtWidgets.QProgressBar, value: Union[int, List[int], tuple[int, int]]): ...


@overload
def set_widget_value(widget: Union[
    QtWidgets.QDoubleSpinBox,
    QtWidgets.QSpinBox],
    value: Union[int, float, list[str], tuple[str, str]]): ...




@overload
def set_widget_value(widget: QtWidgets.QComboBox, value: Union[list, tuple, dict, str, int, float, bool]): ...


@overload
def set_widget_value(widget: QtWidgets.QGridLayout, value: str): ...


@overload
def set_widget_value(widget: Union[
    QtWidgets.QLineEdit,
    QtWidgets.QLCDNumber,
    QtWidgets.QCheckBox
], value: Any): ...


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
    elif isinstance(widget, QtWidgets.QLineEdit):
        widget.setText(str(value))
    elif isinstance(widget, QtWidgets.QProgressBar):
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
def get_widget_value(widget: Union[  # type: ignore
    QtWidgets.QComboBox,
    QtWidgets.QLineEdit]) -> str: ...


@overload
def get_widget_value(widget: Union[  # type: ignore
    QtWidgets.QDoubleSpinBox,
    QtWidgets.QLCDNumber
]) -> float: ...


@overload
def get_widget_value(widget: Union[  # type: ignore
    QtWidgets.QSpinBox,
    QtWidgets.QProgressBar]) -> int: ...



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
    elif isinstance(widget, QtWidgets.QComboBox):
        return widget.currentText()
    elif isinstance(widget, QtWidgets.QLineEdit):
        return widget.text()
    elif isinstance(widget, QtWidgets.QProgressBar):
        return widget.value()
    elif isinstance(widget, QtWidgets.QLCDNumber):
        return float(widget.value())
    elif isinstance(widget, QtWidgets.QCheckBox):
        return widget.isChecked()

    else:
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
        raise NotImplementedError(
            f"{widget} {type(widget)} is not supported yet"
        )

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
    message: str = "",
    error_type: Union[None, Type[Warning]] = None,
    details: Optional[str] = None
) -> None:
    ...

# Overload #2: Exception => NoReturn


@overload
def notify_box(
    message: str,
    error_type: Type[Exception],
    details: Optional[str] = None
) -> NoReturn:
    ...


def notify_box(
    message: str = "",
    error_type: Optional[Type[Exception]] = None,
    details: Optional[str] = None
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


def raise_error(error_type: Type[Exception], message: str) -> NoReturn:
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

    def __init__(self, list_view, items: Dict[str, str] = {}, parent=None):
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
        logging.debug(f'Checked: {checked_items}')
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
            if item.isCheckable() and item.text() != 'Test':
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

    def check_these(self, items):
        """
        Selects the items in the list.
        """
        self.uncheck_all()
        for item in items:
            if not item in self.items.items():
                logging.warning("Item %s not found in list." % item)
                continue
            self.model.itemFromIndex(item).setCheckState(QtCore.Qt.Checked)

    def update(self, new_items: Dict):
        for k, v in new_items.items():
            if not v:
                # Add as a separator
                separator_item = QtGui.QStandardItem(k)
                separator_item.setEnabled(False)  # Non-interactive
                separator_item.setSelectable(False)  # Non-selectable
                separator_item.setCheckable(False)  # Non-checkable
                separator_item.setForeground(QtGui.QBrush(QtCore.Qt.yellow))
                separator_item.setBackground(QtGui.QBrush(QtCore.Qt.blue))   # Different background
                separator_item.setFont(QtGui.QFont("Arial", weight=QtGui.QFont.Bold))  # Bold text
                self.model.appendRow(separator_item)
            else:
                # Add as a regular checkable item
                item = QtGui.QStandardItem(k)
                item.setCheckable(True)
                item.setCheckState(QtCore.Qt.Unchecked)   # Default unchecked
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

