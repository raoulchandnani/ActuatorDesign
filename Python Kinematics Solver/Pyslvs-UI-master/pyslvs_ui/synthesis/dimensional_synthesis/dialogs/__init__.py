# -*- coding: utf-8 -*-

"""'dialogs' module contains contains the dialog of this tab."""

__all__ = [
    'AlgorithmOptionDialog',
    'EditPathDialog',
    'ProgressDialog',
    'PreviewDialog',
    'ChartDialog'
]
__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from .options import AlgorithmOptionDialog
from .edit_path import EditPathDialog
from .progress import ProgressDialog
from .preview import PreviewDialog
from .chart import ChartDialog
