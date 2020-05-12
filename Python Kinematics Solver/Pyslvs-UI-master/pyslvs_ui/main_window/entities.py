# -*- coding: utf-8 -*-

"""This module contains the functions that main window needed."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2020"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from typing import (cast, Tuple, Sequence, Set, FrozenSet, Dict,
                    Counter as Counter_t, Union, Optional)
from abc import abstractmethod, ABC
from collections import Counter
from qtpy.QtCore import Slot
from qtpy.QtWidgets import (
    QDialogButtonBox,
    QDialog,
    QDoubleSpinBox,
    QLabel,
    QHBoxLayout,
    QVBoxLayout,
    QComboBox,
    QMessageBox,
    QInputDialog,
)
from pyslvs import VJoint, VLink, edges_view, SolverSystem, PointArgs, LinkArgs
from pyslvs.graph import Graph
from pyslvs_ui.entities import EditPointDialog, EditLinkDialog
from pyslvs_ui.widgets import (
    AddTable,
    DeleteTable,
    EditPointTable,
    EditLinkTable,
    FixSequenceNumber,
)
from pyslvs_ui.widgets import MainWindowBase

_Coord = Tuple[float, float]


class _ScaleDialog(QDialog):
    """Scale mechanism dialog."""

    def __init__(self, parent: MainWindowBase) -> None:
        super(_ScaleDialog, self).__init__(parent)
        self.setWindowTitle("Scale Mechanism")
        self.main_layout = QVBoxLayout(self)
        self.enlarge = QDoubleSpinBox(self)
        self.shrink = QDoubleSpinBox(self)
        self.__add_option("Enlarge", self.enlarge)
        self.__add_option("Shrink", self.shrink)

        button_box = QDialogButtonBox(self)
        button_box.setStandardButtons(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        button_box.button(QDialogButtonBox.Ok).setEnabled(bool(parent.vpoint_list))
        self.main_layout.addWidget(button_box)

    def __add_option(self, name: str, option: QDoubleSpinBox) -> None:
        """Add widgets for option."""
        layout = QHBoxLayout()
        label = QLabel(name, self)
        option.setValue(1)
        option.setMaximum(10000)
        option.setMinimum(0.01)
        layout.addWidget(label)
        layout.addWidget(option)
        self.main_layout.addLayout(layout)

    def factor(self) -> float:
        """Return scale value."""
        return self.enlarge.value() / self.shrink.value()


class _LinkLengthDialog(QDialog):
    """Link length dialog."""

    def __init__(self, parent: MainWindowBase) -> None:
        super(_LinkLengthDialog, self).__init__(parent)
        self.setWindowTitle("Set Link Length")
        self.main_layout = QVBoxLayout(self)
        layout = QHBoxLayout()
        self.leader = QComboBox(self)
        self.follower = QComboBox(self)
        self.length = QDoubleSpinBox(self)
        layout.addWidget(self.leader)
        layout.addWidget(self.follower)
        layout.addWidget(self.length)
        self.main_layout.addLayout(layout)

        self.vpoints = parent.vpoint_list
        self.vlinks: Dict[str, FrozenSet[int]] = {
            vlink.name: frozenset(vlink.points) for vlink in parent.vlink_list
        }
        self.leader.currentTextChanged.connect(self.__set_follower)
        self.follower.currentTextChanged.connect(self.__set_length)
        self.leader.addItems([f"P{i}" for i in range(len(self.vpoints))])
        self.leader.setCurrentIndex(0)
        self.length.setMaximum(100000)

        button_box = QDialogButtonBox(self)
        button_box.setStandardButtons(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        button_box.button(QDialogButtonBox.Ok).setEnabled(bool(parent.vpoint_list))
        self.main_layout.addWidget(button_box)

    @Slot(str)
    def __set_follower(self, leader: str) -> None:
        """Set follower options."""
        self.follower.clear()
        n = int(leader.replace('P', ''))
        options: Set[int] = set()
        for name, points in self.vlinks.items():
            if name == VLink.FRAME:
                continue
            if n in points:
                options.update(points)
        options.discard(n)
        self.follower.addItems([f"P{i}" for i in options])

    @Slot(str)
    def __set_length(self, follower: str) -> None:
        """Set the current length of two points."""
        if not follower:
            return
        n1 = self.get_leader()
        n2 = int(follower.replace('P', ''))
        self.length.setValue(self.vpoints[n1].distance(self.vpoints[n2]))

    def get_leader(self) -> int:
        """Get current leader."""
        return int(self.leader.currentText().replace('P', ''))

    def get_follower(self) -> int:
        """Get current follower."""
        return int(self.follower.currentText().replace('P', ''))

    def get_length(self) -> float:
        """Get current length."""
        return self.length.value()


class EntitiesMethodInterface(MainWindowBase, ABC):
    """Abstract class for entities methods."""

    @abstractmethod
    def __init__(self):
        """Defined mouse position value on main canvas."""
        super(EntitiesMethodInterface, self).__init__()
        self.mouse_pos_x = 0.
        self.mouse_pos_y = 0.

    def __edit_point(self, row: Union[int, bool] = False) -> None:
        """Edit point function."""
        dlg = EditPointDialog(self.vpoint_list, self.vlink_list, row, self)
        dlg.show()
        if not dlg.exec_():
            dlg.deleteLater()
            return

        row_count = self.entities_point.rowCount()
        type_str = dlg.type_box.currentText().split()[0]
        if type_str != 'R':
            type_str += f":{dlg.angle_box.value() % 360}"
        args = PointArgs(
            ','.join(
                dlg.selected.item(link).text()
                for link in range(dlg.selected.count())
            ),
            type_str,
            dlg.color_box.currentText(),
            dlg.x_box.value(),
            dlg.y_box.value()
        )
        if row is False:
            self.command_stack.beginMacro(f"Add {{Point{row_count}}}")
            self.command_stack.push(AddTable(self.vpoint_list, self.entities_point))
            row = row_count
        else:
            row = dlg.name_box.currentIndex()
            self.command_stack.beginMacro(f"Edit {{Point{row}}}")

        dlg.deleteLater()

        self.command_stack.push(EditPointTable(
            row,
            self.vpoint_list,
            self.vlink_list,
            self.entities_point,
            self.entities_link,
            args
        ))
        self.command_stack.endMacro()

    def __edit_link(self, row: Union[int, bool] = False) -> None:
        """Edit link function."""
        dlg = EditLinkDialog(self.vpoint_list, self.vlink_list, row, self)
        dlg.show()
        if not dlg.exec_():
            dlg.deleteLater()
            return

        name = dlg.name_edit.text()
        args = LinkArgs(name, dlg.color_box.currentText(), ','.join(
            dlg.selected.item(point).text()
            for point in range(dlg.selected.count())
        ))
        if row is False:
            self.command_stack.beginMacro(f"Add {{Link: {name}}}")
            self.command_stack.push(AddTable(self.vlink_list, self.entities_link))
            row = self.entities_link.rowCount() - 1
        else:
            row = dlg.name_box.currentIndex()
            self.command_stack.beginMacro(f"Edit {{Link: {name}}}")

        dlg.deleteLater()

        self.command_stack.push(EditLinkTable(
            row,
            self.vpoint_list,
            self.vlink_list,
            self.entities_point,
            self.entities_link,
            args
        ))
        self.command_stack.endMacro()

    def __get_link_serial_number(self) -> str:
        """Return a new serial number name of link."""
        names = {vlink.name for vlink in self.vlink_list}
        i = 1
        while f"link_{i}" in names:
            i += 1
        return f"link_{i}"

    @Slot(name='on_action_delete_point_triggered')
    def delete_point(self, row: Optional[int] = None) -> None:
        """Push delete point command to stack."""
        if row is None:
            row = self.entities_point.currentRow()
        if row < 0:
            return
        args = self.entities_point.row_data(row)
        args.links = ''
        self.command_stack.beginMacro(f"Delete {{Point{row}}}")
        for i in reversed([
            i for i, (b, d, a) in enumerate(self.inputs_widget.input_pairs())
            if row in {b, d}
        ]):
            self.inputs_widget.remove_var(i)
        self.command_stack.push(EditPointTable(
            row,
            self.vpoint_list,
            self.vlink_list,
            self.entities_point,
            self.entities_link,
            args
        ))
        for i in range(self.entities_link.rowCount()):
            self.command_stack.push(FixSequenceNumber(
                self.vlink_list,
                self.entities_link,
                i,
                row
            ))
        self.command_stack.push(DeleteTable(
            row,
            self.vpoint_list,
            self.entities_point,
            is_rename=True
        ))
        self.inputs_widget.variable_excluding(row)
        self.command_stack.endMacro()
        if self.prefer.auto_remove_link_option:
            self.delete_empty_links()

    @Slot(name='on_action_delete_link_triggered')
    def delete_link(self, row: Optional[int] = None) -> None:
        """Push delete link command to stack.

        Remove link will not remove the points.
        """
        if row is None:
            row = self.entities_link.currentRow()
        if row < 1:
            return
        args = self.entities_link.row_data(row)
        args.points = ''
        self.command_stack.beginMacro(f"Delete {{Link: {self.vlink_list[row].name}}}")
        self.command_stack.push(EditLinkTable(
            row,
            self.vpoint_list,
            self.vlink_list,
            self.entities_point,
            self.entities_link,
            args
        ))
        self.command_stack.push(DeleteTable(
            row,
            self.vlink_list,
            self.entities_link,
            is_rename=False
        ))
        self.command_stack.endMacro()

    def delete_points(self, points: Sequence[int]) -> None:
        """Delete multiple points."""
        if not points:
            return
        self.command_stack.beginMacro(f"Delete points: {sorted(points)}")
        for row in sorted(points, reverse=True):
            self.delete_point(row)
        self.command_stack.endMacro()

    def delete_links(self, links: Sequence[int]) -> None:
        """Delete multiple links."""
        if not links:
            return
        names = ", ".join(self.vlink_list[i].name for i in sorted(links))
        self.command_stack.beginMacro(f"Delete links: [{names}]")
        for row in sorted(links, reverse=True):
            if row == 0:
                continue
            self.delete_link(row)
        self.command_stack.endMacro()

    @Slot(float, float)
    def add_point_by_pos(self, x: float, y: float) -> None:
        """Add point group using alt key."""
        if (
            self.main_panel.currentWidget() is self.synthesis_tab
            and self.synthesis_tab_widget.currentWidget() is self.dimensional_synthesis
        ):
            self.add_target_point()
        else:
            self.add_point(x, y)

    def add_normal_point(self) -> None:
        """Add a point (not fixed)."""
        self.add_point(self.mouse_pos_x, self.mouse_pos_y)

    def add_fixed_point(self) -> None:
        """Add a point (fixed)."""
        self.add_point(self.mouse_pos_x, self.mouse_pos_y, VLink.FRAME, 'Blue')

    def add_point(
        self,
        x: float,
        y: float,
        links: str = "",
        color: str = 'Green',
        type_num: Union[int, VJoint] = VJoint.R,
        angle: float = 0.
    ) -> int:
        """Add an ordinary point. Return the row count of new point."""
        row_count = self.entities_point.rowCount()
        self.command_stack.beginMacro(f"Add {{Point{row_count}}}")
        self.command_stack.push(AddTable(self.vpoint_list, self.entities_point))
        if type_num == VJoint.R:
            type_str = 'R'
        elif type_num == VJoint.P:
            type_str = f'P:{angle}'
        else:
            type_str = f'RP:{angle}'
        self.command_stack.push(EditPointTable(
            row_count,
            self.vpoint_list,
            self.vlink_list,
            self.entities_point,
            self.entities_link,
            PointArgs(links, type_str, color, x, y)
        ))
        self.command_stack.endMacro()
        return row_count

    def add_points(
        self,
        p_attr: Sequence[Tuple[float, float, str, str, int, float]]
    ):
        """Add multiple points."""
        for attr in p_attr:
            self.add_point(*attr)

    def add_points_by_graph(
        self,
        graph: Graph,
        pos: Dict[int, Tuple[float, float]],
        ground_link: Optional[int]
    ):
        """Add points by NetworkX graph and position dict."""
        base_count = self.entities_point.rowCount()
        self.command_stack.beginMacro(
            "Merge mechanism kit from {Number and Type Synthesis}"
        )

        for i in range(len(pos)):
            x, y = pos[i]
            self.add_point(x, y)

        ground: Optional[int] = None
        for link in graph.vertices:
            self.add_link(self.__get_link_serial_number(), 'Blue', [
                base_count + n for n, edge in edges_view(graph) if link in edge
            ])
            if link == ground_link:
                ground = self.entities_link.rowCount() - 1
        self.command_stack.endMacro()
        if ground_link is not None:
            self.constrain_link(ground)

    @Slot(list)
    def add_normal_link(self, points: Sequence[int]) -> None:
        """Add a link."""
        self.add_link(self.__get_link_serial_number(), 'Blue', points)

    def add_link(
        self,
        name: str,
        color: str,
        points: Optional[Sequence[int]] = None
    ) -> None:
        """Push a new link command to stack."""
        if points is None:
            points = []
        args = LinkArgs(name, color, ','.join(f'Point{i}' for i in points))
        self.command_stack.beginMacro(f"Add {{Link: {name}}}")
        self.command_stack.push(AddTable(self.vlink_list, self.entities_link))
        self.command_stack.push(EditLinkTable(
            self.entities_link.rowCount() - 1,
            self.vpoint_list,
            self.vlink_list,
            self.entities_point,
            self.entities_link,
            args
        ))
        self.command_stack.endMacro()

    @Slot(name='on_action_new_point_triggered')
    def new_point(self) -> None:
        """Create a point with arguments."""
        self.__edit_point()

    @Slot(name='on_action_edit_point_triggered')
    def edit_point(self) -> None:
        """Edit a point with arguments."""
        row = self.entities_point.currentRow()
        self.__edit_point(row if row > -1 else 0)

    def lock_points(self) -> None:
        """Turn a group of points to fixed on ground or not."""
        to_fixed = self.action_p_lock.isChecked()
        selected_rows = self.entities_point.selected_rows()
        self.command_stack.beginMacro(
            f"{'Grounded' if to_fixed else 'Ungrounded'} {sorted(selected_rows)}"
        )
        for row in selected_rows:
            new_links = list(self.vpoint_list[row].links)
            if to_fixed:
                if VLink.FRAME not in new_links:
                    new_links.append(VLink.FRAME)
            elif VLink.FRAME in new_links:
                new_links.remove(VLink.FRAME)
            args = self.entities_point.row_data(row)
            args.links = ','.join(s for s in new_links if s)
            self.command_stack.push(EditPointTable(
                row,
                self.vpoint_list,
                self.vlink_list,
                self.entities_point,
                self.entities_link,
                args
            ))
        self.command_stack.endMacro()

    def clone_point(self) -> None:
        """Clone a point (with orange color)."""
        row = self.entities_point.currentRow()
        args = self.entities_point.row_data(row)
        args.color = 'Orange'
        row_count = self.entities_point.rowCount()
        self.command_stack.beginMacro(f"Clone {{Point{row}}} as {{Point{row_count}}}")
        self.command_stack.push(AddTable(self.vpoint_list, self.entities_point))
        self.command_stack.push(EditPointTable(
            row_count,
            self.vpoint_list,
            self.vlink_list,
            self.entities_point,
            self.entities_link,
            args
        ))
        self.command_stack.endMacro()

    @Slot(name="on_action_scale_points_triggered")
    def __set_scale(self) -> None:
        """Scale the mechanism."""
        dlg = _ScaleDialog(self)
        if not dlg.exec_():
            dlg.deleteLater()
            return

        factor = dlg.factor()
        dlg.deleteLater()

        self.command_stack.beginMacro(f"Scale mechanism: {factor}")
        for row in range(self.entities_point.rowCount()):
            args = self.entities_point.row_data(row)
            args.x *= factor
            args.y *= factor
            self.command_stack.push(EditPointTable(
                row,
                self.vpoint_list,
                self.vlink_list,
                self.entities_point,
                self.entities_link,
                args
            ))
        self.command_stack.endMacro()

    @Slot(name='on_action_set_link_length_triggered')
    def __set_link_length(self) -> None:
        """Set link length."""
        dlg = _LinkLengthDialog(self)
        dlg.show()
        if not dlg.exec_():
            return
        data = {(dlg.get_leader(), dlg.get_follower()): dlg.get_length()}
        dlg.deleteLater()
        system = SolverSystem(
            self.vpoint_list,
            {(b, d): a for b, d, a in self.inputs_widget.input_pairs()}
        )
        system.set_data(data)
        try:
            result = system.solve()
        except ValueError:
            QMessageBox.warning(self, "Solved error", "The condition is not valid.")
            return
        self.command_stack.beginMacro(f"Set link length:{set(data)}")
        for row, c in enumerate(result):
            args = self.entities_point.row_data(row)
            if type(c[0]) is float:
                args.x, args.y = cast(_Coord, c)
            else:
                (args.x, args.y), _ = cast(Tuple[_Coord, _Coord], c)
            self.command_stack.push(EditPointTable(
                row,
                self.vpoint_list,
                self.vlink_list,
                self.entities_point,
                self.entities_link,
                args
            ))
        self.command_stack.endMacro()

    @Slot(tuple)
    def set_free_move(
        self,
        args: Sequence[Tuple[int, Tuple[float, float, float]]]
    ):
        """Free move function."""
        points_text = ", ".join(f"Point{c[0]}" for c in args)
        self.command_stack.beginMacro(f"Moved {{{points_text}}}")
        for row, (x, y, angle) in args:
            arg = self.entities_point.row_data(row)
            arg.x = x
            arg.y = y
            if arg.type != 'R':
                angle_tag = arg.type.split(':')[0]
                arg.type = f"{angle_tag}:{angle:.02f}"
            self.command_stack.push(EditPointTable(
                row,
                self.vpoint_list,
                self.vlink_list,
                self.entities_point,
                self.entities_link,
                arg
            ))
        self.command_stack.endMacro()

    @Slot(name='on_action_new_link_triggered')
    def new_link(self) -> None:
        """Create a link with arguments.

        + Last than one point:
            + Create a new link (Dialog)
        + Search method:
            + Find the intersection between points that was including any link.
            + Add the points that is not in the intersection to the link.
        + If no, just create a new link by selected points.
        """
        rows = self.entities_point.selected_rows()
        if not len(rows) > 1:
            self.__edit_link()
            return
        inter: Counter_t[str] = Counter()
        for p in rows:
            inter.update(self.vpoint_list[p].links)
        name = max(inter, key=inter.get) if inter else ""
        if inter[name] < 2:
            self.add_normal_link(rows)
            return
        row = self.entities_link.find_name(name)
        self.command_stack.beginMacro(f"Edit {{Link: {name}}}")
        args = self.entities_link.row_data(row)
        points = set(self.entities_link.get_points(row))
        points.update(rows)
        args.points = ','.join(f'Point{p}' for p in points)
        self.command_stack.push(EditLinkTable(
            row,
            self.vpoint_list,
            self.vlink_list,
            self.entities_point,
            self.entities_link,
            args
        ))
        self.command_stack.endMacro()

    @Slot(name='on_action_edit_link_triggered')
    def edit_link(self) -> None:
        """Edit a link with arguments."""
        self.__edit_link(self.entities_link.currentRow())

    @Slot()
    def release_ground(self) -> None:
        """Clone ground to a new link, then make ground no points."""
        name = self.__get_link_serial_number()
        args = LinkArgs(name, 'Blue', self.entities_link.item(0, 2).text())
        self.command_stack.beginMacro(f"Release ground to {{Link: {name}}}")
        # Free all points.
        self.command_stack.push(EditLinkTable(
            0,
            self.vpoint_list,
            self.vlink_list,
            self.entities_point,
            self.entities_link,
            LinkArgs(VLink.FRAME, 'White', '')
        ))
        # Create new link.
        self.command_stack.push(AddTable(self.vlink_list, self.entities_link))
        self.command_stack.push(EditLinkTable(
            self.entities_link.rowCount() - 1,
            self.vpoint_list,
            self.vlink_list,
            self.entities_point,
            self.entities_link,
            args
        ))
        self.command_stack.endMacro()

    @Slot()
    def constrain_link(self, row1: Optional[int] = None, row2: int = 0) -> None:
        """Turn a link to ground, then delete this link."""
        if row1 is None:
            row1 = self.entities_link.currentRow()
        vlink1 = self.vlink_list[row1]
        link_args = self.entities_link.row_data(row1)
        link_args.points = ''
        new_points = sorted(set(self.vlink_list[0].points) | set(vlink1.points))
        base_args = self.entities_link.row_data(row2)
        base_args.points = ','.join(f"Point{e}" for e in new_points if e)
        self.command_stack.beginMacro(f"Constrain {{Link: {vlink1.name}}} to ground")
        # Turn to ground.
        self.command_stack.push(EditLinkTable(
            row2,
            self.vpoint_list,
            self.vlink_list,
            self.entities_point,
            self.entities_link,
            base_args
        ))
        # Free all points and delete the link.
        self.command_stack.push(EditLinkTable(
            row1,
            self.vpoint_list,
            self.vlink_list,
            self.entities_point,
            self.entities_link,
            link_args
        ))
        self.command_stack.push(DeleteTable(
            row1,
            self.vlink_list,
            self.entities_link,
            is_rename=False
        ))
        self.command_stack.endMacro()

    @Slot()
    def delete_selected_points(self) -> None:
        """Delete the selected points."""
        self.delete_points(self.entities_point.selected_rows())

    @Slot()
    def delete_selected_links(self) -> None:
        """Delete the selected links."""
        self.delete_links(self.entities_link.selected_rows())

    @Slot()
    def delete_empty_links(self) -> None:
        """Delete empty link names."""
        self.delete_links([
            i for i, vlink in enumerate(self.vlink_list)
            if vlink.name != VLink.FRAME and len(vlink.points) < 2
        ])

    def set_coords_as_current(self) -> None:
        """Update points position as current coordinate."""
        self.set_free_move(tuple(
            (row, (vpoint.cx, vpoint.cy, vpoint.angle))
            for row, vpoint in enumerate(self.vpoint_list)
        ))

    @Slot()
    def point_alignment(self) -> None:
        """Alignment function."""
        selected_rows = self.entities_point.selected_rows()
        if not selected_rows:
            QMessageBox.warning(
                self,
                "Points alignment",
                "No selected points with this operation."
            )
            return
        if self.alignment_mode == 0:
            axis = "x"
        elif self.alignment_mode == 1:
            axis = "y"
        else:
            raise ValueError("no such alignment option")
        value, ok = QInputDialog.getDouble(
            self,
            f"Set {axis} axis",
            f"Align the selected points into {axis} axis:",
            0, -9999, 9999, 4)
        if not ok:
            return
        self.command_stack.beginMacro(f"Align points with {axis}")
        for row in selected_rows:
            args = self.entities_point.row_data(row)
            if self.alignment_mode == 0:
                args.x = value
            elif self.alignment_mode == 1:
                args.y = value
            else:
                raise ValueError("no such alignment option")
            self.command_stack.push(EditPointTable(
                row,
                self.vpoint_list,
                self.vlink_list,
                self.entities_point,
                self.entities_link,
                args
            ))
        self.command_stack.endMacro()
