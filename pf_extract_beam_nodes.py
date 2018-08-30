"""
This script filters all cells for polylines and extracts the first and last
nodes from that cell (nodal points for beams). It also keeps all point data
values.

Use with programmable filter in paraview and execute with:
>> execfile('path to file')
"""

# Import paraview moudles.
from paraview import vtk

# Get input and output objects (for now output is equal to input).
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
                if array.GetNumberOfComponents() == 1:
                    array.InsertNextTuple1(pdi.GetPointData().GetArray(j).GetTuple1(index))
                elif array.GetNumberOfComponents() == 3:
                    array.InsertNextTuple3(*pdi.GetPointData().GetArray(j).GetTuple3(index)[:3])
                else:
                    raise ValueError('Only 1 and 3 tuples allowed.')

# Reset output.
pdo.Initialize()

# Add the nodal points to the output.
pdo.SetPoints(vtk_points)

# Add point data to the output.
for array in point_data_arrays:
    pdo.GetPointData().AddArray(array)
