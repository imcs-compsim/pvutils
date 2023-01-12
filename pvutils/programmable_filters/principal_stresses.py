"""
This script adds a field variable with the node and cell ids respectively.

Use with programmable filter in paraview and execute with:
>> exec(open('path to file').read())
"""

# Import modules
from paraview import vtk
import numpy as np

# Get input and output objects
pdi = self.GetInput()
pdo = self.GetOutput()

# Add old fields
for i in range(pdi.GetCellData().GetNumberOfArrays()):
    pdo.GetCellData().AddArray(pdi.GetCellData().GetArray(i))
for i in range(pdi.GetPointData().GetNumberOfArrays()):
    pdo.GetPointData().AddArray(pdi.GetPointData().GetArray(i))

# Loop over nodes and perform the eigenvalue calculations
for i in range(pdi.GetPointData().GetNumberOfArrays()):
    stress_array = pdi.GetPointData().GetArray(i)
    if stress_array.GetName() == "nodal_cauchy_stresses_xyz":
        break
else:
    stress_array = None

# Loop over points and do the eigenvalue calculation
principal_stresses = vtk.vtkDoubleArray()
principal_stresses.SetName("principal_stresses")
principal_stresses.SetNumberOfComponents(3)
directions = [vtk.vtkDoubleArray() for i in range(3)]
for i, direction in enumerate(directions):
    direction.SetName("principal_direction_{}".format(i + 1))
    direction.SetNumberOfComponents(3)
for i in range(pdi.GetNumberOfPoints()):
    point_stress = stress_array.GetTuple(i)
    stress_tensor = np.array(
        [
            [point_stress[0], point_stress[3], point_stress[5]],
            [point_stress[3], point_stress[1], point_stress[4]],
            [point_stress[5], point_stress[4], point_stress[2]],
        ]
    )
    eigen_values, eigen_vectors = np.linalg.eig(stress_tensor)
    principal_stresses.InsertNextTuple3(*eigen_values)
    for i, direction in enumerate(directions):
        direction.InsertNextTuple3(*(eigen_vectors[i]))

pdo.GetPointData().AddArray(principal_stresses)
for direction in directions:
    pdo.GetPointData().AddArray(direction)
