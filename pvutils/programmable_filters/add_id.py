# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# PVUTILS: Utility tools for the ParaView python interface
#
# Copyright 2023 Ivo Steinbrecher & Matthias Mayr
#                 Institute for Mathematics and Computer-based Simulation
#                 University of the Bundeswehr Munich
#                 https://www.unibw.de/imcs-en
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
This script adds a field variable with the node and cell ids respectively.

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

# Loop over cells and nodes and get ids.
id_items = [
    [
        "cell_id",
        vtk.vtkDoubleArray(),
        pdi.GetNumberOfCells(),
        pdo.GetCellData().AddArray,
    ],
    [
        "node_id",
        vtk.vtkDoubleArray(),
        pdi.GetNumberOfPoints(),
        pdo.GetPointData().AddArray,
    ],
]
for name, array, n_items, add_function in id_items:
    array.SetName(name)
    for i in range(n_items):
        array.InsertNextTuple1(i)
    add_function(array)
