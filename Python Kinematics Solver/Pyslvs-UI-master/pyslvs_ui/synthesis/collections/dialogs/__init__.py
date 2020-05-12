# -*- coding: utf-8 -*-

"""'dialogs' module contains
the dialog of this tab.
"""

__all__ = [
    'CollectionsDialog',
    'CustomsDialog',
    'TargetsDialog',
    'list_texts',
]
__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from .collections import CollectionsDialog
from .customs import CustomsDialog
from .targets import TargetsDialog, list_texts
