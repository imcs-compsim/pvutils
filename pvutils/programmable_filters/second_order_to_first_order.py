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
This converts second order cells to first order cells.

Use with programmable filter in paraview and execute with:
>> exec(open("<path to file>").read())
"""

# Import paraview moudles.
from paraview import vtk

# Get input and output objects.
pdi = self.GetInput()
pdo = self.GetOutput()

# Add old fields.
for i in range(pdi.GetCellData().GetNumberOfArrays()):
    pdo.GetCellData().AddArray(pdi.GetCellData().GetArray(i))
for i in range(pdi.GetPointData().GetNumberOfArrays()):
    pdo.GetPointData().AddArray(pdi.GetPointData().GetArray(i))

# Loop over cells and convert them.
n_cells = pdi.GetNumberOfCells()
new_cells = []
for i in range(n_cells):

    old_cell = pdi.GetCell(i)
    n_points_old_cell = old_cell.GetNumberOfPoints()

    if (
        type(old_cell) == vtk.vtkQuadraticHexahedron
        or type(old_cell) == vtk.vtkTriQuadraticHexahedron
    ):
        n_points_new_cell = 8
        new_cell = vtk.vtkHexahedron()
    elif type(old_cell) == vtk.vtkQuadraticTetra:
        n_points_new_cell = 4
        new_cell = vtk.vtkTetra()
    else:
        n_points_new_cell = n_points_old_cell
        new_cell = type(old_cell)()

    new_cell.GetPointIds().SetNumberOfIds(n_points_new_cell)
    for i_point in range(n_points_new_cell):
        new_cell.GetPointIds().SetId(i_point, old_cell.GetPointId(i_point))
    new_cells.append(new_cell)

# Add all new cells to the output data.
pdo.Allocate(len(new_cells), 1)
for cell in new_cells:
    pdo.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
