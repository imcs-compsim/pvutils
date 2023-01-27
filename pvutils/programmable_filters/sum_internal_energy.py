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
This script sum the cell value element_internal_energy for all cells.

Use with programmable filter in paraview and execute with:
>> exec(open("<path to file>").read())
"""

# Import paraview moudles.
from paraview import vtk

# Get input and output objects.
pdi = self.GetInput()
pdo = self.GetOutput()

# Calculate the global internal energy.
global_internal_energy = 0.0
energy_array = pdi.GetCellData().GetArray("element_internal_energy")
for i_cell in range(energy_array.GetNumberOfTuples()):
    global_internal_energy += energy_array.GetTuple1(i_cell)

# Create array for the global energy data.
array = vtk.vtkDoubleArray()
array.SetName("global_internal_energy")
array.SetNumberOfComponents(1)
array.InsertNextTuple1(global_internal_energy)

# Add to output object.
pdo.GetFieldData().AddArray(array)
