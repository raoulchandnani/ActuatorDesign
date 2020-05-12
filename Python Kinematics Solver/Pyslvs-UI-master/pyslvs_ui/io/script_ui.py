# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'pyslvs_ui/io/script.ui'
#
# Created by: PyQt5 UI code generator 5.14.1
#
# WARNING! All changes made in this file will be lost!


from qtpy import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.setEnabled(True)
        Dialog.resize(533, 564)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/icons/script.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        Dialog.setWindowIcon(icon)
        Dialog.setAutoFillBackground(True)
        Dialog.setSizeGripEnabled(True)
        Dialog.setModal(True)
        self.main_layout = QtWidgets.QVBoxLayout(Dialog)
        self.main_layout.setObjectName("main_layout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.show_qrcode = QtWidgets.QPushButton(Dialog)
        self.show_qrcode.setAutoDefault(False)
        self.show_qrcode.setObjectName("show_qrcode")
        self.horizontalLayout.addWidget(self.show_qrcode)
        self.copy = QtWidgets.QPushButton(Dialog)
        self.copy.setContextMenuPolicy(QtCore.Qt.NoContextMenu)
        self.copy.setAutoDefault(False)
        self.copy.setObjectName("copy")
        self.horizontalLayout.addWidget(self.copy)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.button_box = QtWidgets.QDialogButtonBox(Dialog)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.button_box.sizePolicy().hasHeightForWidth())
        self.button_box.setSizePolicy(sizePolicy)
        self.button_box.setOrientation(QtCore.Qt.Horizontal)
        self.button_box.setStandardButtons(QtWidgets.QDialogButtonBox.Close)
        self.button_box.setObjectName("button_box")
        self.horizontalLayout.addWidget(self.button_box)
        self.main_layout.addLayout(self.horizontalLayout)

        self.retranslateUi(Dialog)
        self.button_box.rejected.connect(Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        self.show_qrcode.setText(_translate("Dialog", "Show QR code"))
        self.copy.setText(_translate("Dialog", "Copy"))
        self.button_box.setWhatsThis(_translate("Dialog", "Click to exit"))
from pyslvs_ui import icons_rc
