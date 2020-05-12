# -*- coding: utf-8 -*-

"""'widgets' module contains the custom widgets
that design without Qt designer.
"""

__all__ = [
    'MainWindowBase',
    'Preferences',
    'PointArgs',
    'LinkArgs',
    'PointTableWidget',
    'LinkTableWidget',
    'QRotatableView',
    'AddTable',
    'AddPath',
    'AddStorage',
    'AddStorageName',
    'AddInput',
    'ClearStorageName',
    'DeletePath',
    'DeleteStorage',
    'DeleteTable',
    'DeleteInput',
    'EditPointTable',
    'EditLinkTable',
    'FixSequenceNumber',
]
__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from .main_base import MainWindowBase, Preferences
from .tables import PointArgs, LinkArgs, PointTableWidget, LinkTableWidget
from .rotatable import QRotatableView
from .undo_redo import (
    AddTable,
    AddPath,
    AddStorage,
    AddStorageName,
    AddInput,
    ClearStorageName,
    DeletePath,
    DeleteStorage,
    DeleteTable,
    DeleteInput,
    EditPointTable,
    EditLinkTable,
    FixSequenceNumber,
)
