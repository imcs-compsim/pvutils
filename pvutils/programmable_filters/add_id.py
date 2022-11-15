"""
This script adds a field variable with the node and cell ids respectively.

Use with programmable filter in paraview and execute with:
>> exec(open('path to file').read())
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
