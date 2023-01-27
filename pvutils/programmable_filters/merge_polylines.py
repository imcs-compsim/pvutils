# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# PVUTILS: Utility tools for the ParaView python interface
#
# Copyright 2023 Ivo Steinbrecher & Matthias Mayr
#                Institute for Mathematics and Computer-based Simulation
#                University of the Bundeswehr Munich
#                https://www.unibw.de/imcs-en
#
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# -----------------------------------------------------------------------------
"""
This script creates a single poly line for a beam filament. For this to work,
the connectivity of the nodes has to be resolved (CleanToGrid filter). Also it
is advisable to use the cell to point data filter before applying this one.

Use with programmable filter in paraview and execute with:
>> exec(open("<path to file>").read())
"""

# Import paraview moudles.
import paraview
from paraview import vtk
import numpy as np

# Get "keyword arguments".
# Default value, if the angle of a corner is more than 45 degrees, the polyline
# will be split there.
max_angle = 0.5 * np.pi
if hasattr(paraview, "programmable_filter_kwargs"):
    if "max_angle" in paraview.programmable_filter_kwargs[kwargs_id].keys():
        tmp = paraview.programmable_filter_kwargs[kwargs_id]["max_angle"]
        if tmp is not None:
            max_angle = tmp

# Get input and output objects.
pdi = self.GetInput()
pdo = self.GetOutput()

# Reset output.
pdo.Initialize()

# Add the points to the output.
pdo.SetPoints(pdi.GetPoints())
for i in range(pdi.GetPointData().GetNumberOfArrays()):
    pdo.GetPointData().AddArray(pdi.GetPointData().GetArray(i))

# Check that all cells are poly lines.
n_cells = pdi.GetNumberOfCells()
for i in range(n_cells):
    if not pdi.GetCellType(i) == 4:
        raise ValueError("Only poly lines (vtk type 4) are supported")


def find_connected_cells(pdi, cell_id):
    """
    Start with poly line "cell_id" and search all poly lines connected to that
    one. Also return all points of those lines in order for the combined poly
    line.
    """

    global max_angle

    def vtk_id_to_list(vtk_id_list):
        return [
            int(vtk_id_list.GetId(i_id)) for i_id in range(vtk_id_list.GetNumberOfIds())
        ]

    def add_cell_recursive(connected_cell_points, old_cells, initial_point_id):
        import numpy as np

        pdi.GetPointCells(initial_point_id, id_list)
        cell_connectivity = vtk_id_to_list(id_list)
        new_cell_ids = [
            cell_id for cell_id in cell_connectivity if cell_id not in old_cells
        ]
        if len(cell_connectivity) > 2 or len(new_cell_ids) == 0:
            # In this case we either have a bifurcation point or are at the end
            # the poly line.
            return connected_cell_points, old_cells

        # Add this cell and its points (in correct order).
        new_cell_id = new_cell_ids[0]
        new_cell_point_ids = vtk_id_to_list(pdi.GetCell(new_cell_id).GetPointIds())
        if new_cell_point_ids[0] == connected_cell_points[-1]:
            # First point of this cell is added to the last point of the last
            # cell.
            extend = True
        elif new_cell_point_ids[-1] == connected_cell_points[-1]:
            # Last point of this cell is added to the last point of the last
            # cell.
            extend = True
            new_cell_point_ids.reverse()
        elif new_cell_point_ids[0] == connected_cell_points[0]:
            # First point of this cell is added to the first point of the last
            # cell.
            extend = False
            new_cell_point_ids.reverse()
        elif new_cell_point_ids[-1] == connected_cell_points[0]:
            # Last point of this cell is added to the first point of the last
            # cell.
            extend = False
        else:
            raise ValueError("This should not happen")

        # Check the angel of the corner and decide whether or not to merge
        # this.
        if extend:
            corner_point_ids = connected_cell_points[-2:] + [new_cell_point_ids[1]]
        else:
            corner_point_ids = [new_cell_point_ids[-2]] + connected_cell_points[:2]
        points = [np.array(pdi.GetPoint(j)) for j in corner_point_ids]
        vec_1 = points[1] - points[0]
        vec_1 = vec_1 / np.linalg.norm(vec_1)
        vec_2 = points[2] - points[1]
        vec_2 = vec_2 / np.linalg.norm(vec_2)
        dot = np.dot(vec_1, vec_2)
        if 1.0 <= dot and dot < 1.0 + 1e-12:
            dot = 1.0
        angle = np.arccos(dot)
        if np.abs(angle) > max_angle:
            # Angle between beam elements is more than the maximum angle,
            # which indicates an edge and not a continuous polyline.
            return connected_cell_points, old_cells

        # Extend the merged poly line points.
        if extend:
            connected_cell_points.extend(new_cell_point_ids[1:])
            next_start_index = -1
        else:
            connected_cell_points = new_cell_point_ids[:-1] + connected_cell_points
            next_start_index = 0

        old_cells.append(new_cell_id)
        connected_cell_points, old_cells = add_cell_recursive(
            connected_cell_points, old_cells, new_cell_point_ids[next_start_index]
        )

        return connected_cell_points, old_cells

    id_list = vtk.vtkIdList()
    connected_cell_points = vtk_id_to_list(pdi.GetCell(cell_id).GetPointIds())
    old_cells = [cell_id]

    start_id = connected_cell_points[0]
    end_id = connected_cell_points[-1]
    connected_cell_points, old_cells = add_cell_recursive(
        connected_cell_points, old_cells, start_id
    )
    connected_cell_points, old_cells = add_cell_recursive(
        connected_cell_points, old_cells, end_id
    )

    return old_cells, connected_cell_points


# Start with the first poly line and search all poly lines connected to that
# line and so on. Then do the same with the first poly line that was not found
# and so on.
cell_ids = list(range(n_cells))
new_cells = []
while True:
    # Take the next available cell and look for all connected cells.
    for i in cell_ids:
        if i is not None:
            next_id = i
            break
    else:
        break
    old_cells, connected_cell_points = find_connected_cells(pdi, next_id)

    # Mark the found cell IDs.
    for found_cell_id in old_cells:
        cell_ids[found_cell_id] = None

    # Create the poly line for this beam.
    new_cell = vtk.vtkPolyLine()
    new_cell.GetPointIds().SetNumberOfIds(len(connected_cell_points))
    for i, index in enumerate(connected_cell_points):
        new_cell.GetPointIds().SetId(i, index)
    new_cells.append(new_cell)


# Add all new cells to the output data.
pdo.Allocate(len(new_cells), 1)
for cell in new_cells:
    pdo.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
