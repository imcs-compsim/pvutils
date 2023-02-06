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
This script filters all cells for polylines and extracts the first and last
nodes from that cell (nodal points for beams). It also keeps all point data
values.

Use with programmable filter in paraview and execute with:
>> exec(open("<path to file>").read())
"""

# Import paraview moudles.
from paraview import vtk

# Get input and output objects.
pdi = self.GetInput()
pdo = self.GetOutput()

# Get empty point data arrays.
n_point_data_arrays = pdi.GetPointData().GetNumberOfArrays()
point_data_arrays = []
for i in range(n_point_data_arrays):
    name = pdi.GetPointData().GetArray(i).GetName()
    dim = pdi.GetPointData().GetArray(i).GetNumberOfComponents()

    array = vtk.vtkDoubleArray()
    array.SetName(name)
    array.SetNumberOfComponents(dim)
    point_data_arrays.append(array)

# Loop over cells and get the start and endpoints.
n_cells = pdi.GetNumberOfCells()
vtk_points = vtk.vtkPoints()
for i in range(n_cells):
    if isinstance(pdi.GetCell(i), vtk.vtkPolyLine):
        # Get the first and last node indices on the polyline.
        n_cell_nodes = pdi.GetCell(i).GetNumberOfPoints()
        id1 = pdi.GetCell(i).GetPointId(0)
        id2 = pdi.GetCell(i).GetPointId(n_cell_nodes - 1)

        # Add the points to the output.
        for index in [id1, id2]:
            coord = pdi.GetPoint(index)
            x, y, z = coord[:3]
            vtk_points.InsertNextPoint(x, y, z)

            # Add all point data.
            for j, array in enumerate(point_data_arrays):
                array.InsertNextTuple(pdi.GetPointData().GetArray(j).GetTuple(index))

# Reset output.
pdo.Initialize()

# Add the nodal points to the output.
pdo.SetPoints(vtk_points)

# Add point data to the output.
for array in point_data_arrays:
    pdo.GetPointData().AddArray(array)
