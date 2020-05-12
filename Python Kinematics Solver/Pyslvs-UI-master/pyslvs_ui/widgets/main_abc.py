# -*- coding: utf-8 -*-

"""Predefined methods of main window."""

from __future__ import annotations

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from typing import Tuple, List, Sequence, Dict, Optional, Any
from abc import abstractmethod
from qtpy.QtCore import Qt, QPoint
from qtpy.QtWidgets import QMainWindow
from pyslvs import VPoint
from pyslvs.graph import Graph
from pyslvs_ui.qt_patch import QABCMeta
from .main_ui import Ui_MainWindow

_Coord = Tuple[float, float]


class MainWindowABC(QMainWindow, Ui_MainWindow, metaclass=QABCMeta):
    """Main window abstract class."""

    @abstractmethod
    def __init__(self):
        super(MainWindowABC, self).__init__()
        self.setupUi(self)
        self.setAttribute(Qt.WA_DeleteOnClose)

    @abstractmethod
    def command_reload(self, index: int) -> None:
        ...

    @abstractmethod
    def new_point(self) -> None:
        ...

    @abstractmethod
    def add_normal_point(self) -> None:
        ...

    @abstractmethod
    def add_fixed_point(self) -> None:
        ...

    @abstractmethod
    def edit_point(self) -> None:
        ...

    @abstractmethod
    def delete_selected_points(self) -> None:
        ...

    @abstractmethod
    def lock_points(self) -> None:
        ...

    @abstractmethod
    def new_link(self) -> None:
        ...

    @abstractmethod
    def edit_link(self) -> None:
        ...

    @abstractmethod
    def delete_selected_links(self) -> None:
        ...

    @abstractmethod
    def delete_empty_links(self) -> None:
        ...

    @abstractmethod
    def constrain_link(self) -> None:
        ...

    @abstractmethod
    def release_ground(self) -> None:
        ...

    @abstractmethod
    def add_target_point(self) -> None:
        ...

    @abstractmethod
    def set_free_move(self, args: Sequence[Tuple[int, Tuple[float, float, float]]]) -> None:
        ...

    @abstractmethod
    def add_point_by_pos(self, x: float, y: float) -> None:
        ...

    @abstractmethod
    def set_mouse_pos(self, x: float, y: float) -> None:
        ...

    @abstractmethod
    def get_back_position(self) -> None:
        ...

    @abstractmethod
    def solve(self) -> None:
        ...

    @abstractmethod
    def resolve(self) -> None:
        ...

    @abstractmethod
    def enable_mechanism_actions(self) -> None:
        ...

    @abstractmethod
    def clone_point(self) -> None:
        ...

    @abstractmethod
    def copy_coord(self) -> None:
        ...

    @abstractmethod
    def copy_points_table(self) -> None:
        ...

    @abstractmethod
    def copy_links_table(self) -> None:
        ...

    @abstractmethod
    def canvas_context_menu(self, point: QPoint) -> None:
        ...

    @abstractmethod
    def point_context_menu(self, point: QPoint) -> None:
        ...

    @abstractmethod
    def link_context_menu(self, point: QPoint) -> None:
        ...

    @abstractmethod
    def preview_path(
        self,
        auto_preview: List[List[_Coord]],
        slider_auto_preview: Dict[int, List[_Coord]],
        vpoints: Sequence[VPoint]
    ) -> None:
        ...

    @abstractmethod
    def reload_canvas(self) -> None:
        ...

    @abstractmethod
    def output_to(self, format_name: str, format_choose: Sequence[str]) -> str:
        ...

    @abstractmethod
    def right_input(self) -> bool:
        ...

    @abstractmethod
    def set_coords_as_current(self) -> None:
        ...

    @abstractmethod
    def dof(self) -> int:
        ...

    @abstractmethod
    def save_reply_box(self, title: str, file_name: str) -> None:
        ...

    @abstractmethod
    def input_from(
        self,
        format_name: str,
        format_choose: Sequence[str]
    ) -> str:
        ...

    @abstractmethod
    def input_from_multiple(
        self,
        format_name: str,
        format_choose: Sequence[str]
    ) -> List[str]:
        ...

    @abstractmethod
    def get_graph(self) -> Tuple[
        Graph,
        List[int],
        List[Tuple[int, int]],
        Dict[int, _Coord],
        Dict[int, int],
        Dict[int, int]
    ]:
        ...

    @abstractmethod
    def get_configure(self) -> Dict[str, Any]:
        ...

    @abstractmethod
    def project_no_save(self) -> None:
        ...

    @abstractmethod
    def project_saved(self) -> None:
        ...

    @abstractmethod
    def merge_result(self, expr: str, path: Sequence[Sequence[_Coord]]) -> None:
        ...

    @abstractmethod
    def check_file_changed(self) -> bool:
        ...

    @abstractmethod
    def get_storage(self) -> Dict[str, str]:
        ...

    @abstractmethod
    def add_empty_links(self, link_color: Dict[str, str]) -> None:
        ...

    @abstractmethod
    def parse_expression(self, expr: str) -> None:
        ...

    @abstractmethod
    def add_multiple_storage(self, exprs: Dict[str, str]) -> None:
        ...

    @abstractmethod
    def clear(self) -> None:
        ...

    @abstractmethod
    def add_points(
        self,
        p_attr: Sequence[Tuple[float, float, str, str, int, float]]
    ) -> None:
        ...

    @abstractmethod
    def set_window_title_full_path(self) -> None:
        ...

    @abstractmethod
    def import_pmks_url(self) -> None:
        ...

    @abstractmethod
    def save_picture_clipboard(self) -> None:
        ...

    @abstractmethod
    def show_expr(self) -> None:
        ...

    @abstractmethod
    def py_script(self) -> None:
        ...

    @abstractmethod
    def export_dxf(self) -> None:
        ...

    @abstractmethod
    def export_slvs(self) -> None:
        ...

    @abstractmethod
    def save_pmks(self) -> None:
        ...

    @abstractmethod
    def export_image(self) -> None:
        ...

    @abstractmethod
    def show_overview(self) -> None:
        ...

    @staticmethod
    @abstractmethod
    def new() -> MainWindowABC:
        ...

    @abstractmethod
    def point_alignment(self) -> None:
        ...

    @abstractmethod
    def get_expression(
        self,
        points: Optional[Sequence[VPoint]] = None,
        indent: int = -1
    ) -> str:
        ...
