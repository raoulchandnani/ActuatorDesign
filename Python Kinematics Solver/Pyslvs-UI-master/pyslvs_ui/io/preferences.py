# -*- coding: utf-8 -*-

"""This module contains the preferences dialog."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from typing import Optional
from dataclasses import fields, Field
from pyslvs_ui.widgets import MainWindowBase
from qtpy.QtCore import Slot
from qtpy.QtWidgets import (
    QDialog,
    QLineEdit,
    QSpinBox,
    QDoubleSpinBox,
    QCheckBox,
    QComboBox,
    QDialogButtonBox,
    QMessageBox,
)
from qtpy.QtGui import QCloseEvent
from pyslvs_ui.info import kernel_list
from pyslvs_ui.widgets import Preferences
from .format_editor import PROJECT_FORMAT
from .preference_ui import Ui_Dialog


class PreferencesDialog(QDialog, Ui_Dialog):
    """Preference dialog."""

    def __init__(self, parent: MainWindowBase):
        super(PreferencesDialog, self).__init__(parent)
        self.setupUi(self)
        self.file_type_option.addItems(PROJECT_FORMAT)
        self.planar_solver_option.addItems(kernel_list)
        self.path_preview_option.addItems(kernel_list + ("Same as solver kernel",))
        self.prefer = parent.prefer.copy()
        self.prefer_applied = parent.prefer.copy()

        self.default_zoom_option.setMaximum(parent.zoom_bar.maximum())
        self.default_zoom_option.setMinimum(parent.zoom_bar.minimum())
        self.accepted.connect(self.__get_settings)
        self.button_box.button(QDialogButtonBox.Apply).clicked.connect(self.__get_settings)
        self.button_box.button(QDialogButtonBox.Cancel).clicked.connect(self.__cancel)
        self.button_box.button(QDialogButtonBox.RestoreDefaults).clicked.connect(self.__reset)
        self.__load_settings()

    @Slot()
    def __reset(self) -> None:
        """Reset user options."""
        self.prefer.reset()
        self.__load_settings()

    @Slot()
    def __load_settings(self) -> None:
        """Load settings on UI."""
        for field in fields(self.prefer):  # type: Field
            widget = getattr(self, field.name)
            value = getattr(self.prefer, field.name)
            if type(widget) is QSpinBox or type(widget) is QDoubleSpinBox:
                widget.setValue(value)
            elif type(widget) is QLineEdit:
                widget.setText(value)
            elif type(widget) is QCheckBox:
                widget.setChecked(value)
            elif type(widget) is QComboBox:
                widget.setCurrentIndex(value)

    @Slot()
    def __get_settings(self, prefer: Optional[Preferences] = None) -> None:
        """Save settings after clicked apply."""
        if prefer is None:
            prefer = self.prefer_applied
        for field in fields(prefer):  # type: Field
            widget = getattr(self, field.name)
            if type(widget) is QSpinBox or type(widget) is QDoubleSpinBox:
                setattr(prefer, field.name, widget.value())
            elif type(widget) is QLineEdit:
                setattr(prefer, field.name, widget.text())
            elif type(widget) is QCheckBox:
                setattr(prefer, field.name, widget.isChecked())
            elif type(widget) is QComboBox:
                setattr(prefer, field.name, widget.currentIndex())

    @Slot()
    def __cancel(self) -> None:
        """Cancel button clicked."""
        if self.__cancel_check():
            self.reject()

    def closeEvent(self, event: QCloseEvent):
        if self.__cancel_check():
            event.accept()
        else:
            event.ignore()

    def __cancel_check(self) -> bool:
        """Ask for saving options."""
        self.__get_settings(self.prefer)
        if sum(1 for _ in self.prefer_applied.diff(self.prefer)) < 1:
            return True
        reply = QMessageBox.question(
            self,
            "Option changed",
            "Do you want to save the changes?",
            QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel
        )
        if reply == QMessageBox.Save:
            self.prefer_applied = self.prefer
            return True
        elif reply == QMessageBox.Discard:
            return True
        else:
            return False
