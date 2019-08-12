"""
This script sum the cell value element_internal_energy for all cells.

Use with programmable filter in paraview and execute with:
>> execfile('path to file')
"""

# Import paraview moudles.
from paraview import vtk

# Get input and output objects.
pdi = self.GetInput()
pdo = self.GetOutput()

# Calculate the global internal energy.
global_internal_energy = 0.
energy_array = pdi.GetCellData().GetArray('element_internal_energy')
for i_cell in range(energy_array.GetNumberOfTuples()):
    global_internal_energy += energy_array.GetTuple1(i_cell)

# Create array for the global energy data.
array = vtk.vtkDoubleArray()
array.SetName('global_internal_energy')
array.SetNumberOfComponents(1)
array.InsertNextTuple1(global_internal_energy)

# Add to output object.
pdo.GetFieldData().AddArray(array)
