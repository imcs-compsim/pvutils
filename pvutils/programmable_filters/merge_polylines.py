"""
This script creates a single poly line for a beam filament. For this to work,
the connectivity of the nodes has to be resolved (CleanToGrid filter). Also it
is advisable to use the cell to point data filter before applying this one.

Use with programmable filter in paraview and execute with:
>> execfile('path to file')
"""

# Import paraview moudles.
from paraview import vtk

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
        raise ValueError('Only polylines (vtk type 4) are supported')

# Check that no point connects more than one cell.
n_points = pdi.GetNumberOfPoints()
id_list = vtk.vtkIdList()
for i in range(n_points):
    pdi.GetPointCells(i, id_list)
    if id_list.GetNumberOfIds() > 2:
        raise ValueError('Bifurcations are not supported')


def find_connected_cells(pdi, cell_id):
    """
    Start with poly line "cell_id" and search all poly lines connected to that
    one. Also return all points of those lines in order for the combined poly
    line.
    """

    def vtk_id_to_list(vtk_id_list):
        return [int(vtk_id_list.GetId(i_id)) for i_id in
            range(vtk_id_list.GetNumberOfIds())]

    def add_cell_recursive(connected_cell_points, old_cells, initial_point_id):
        pdi.GetPointCells(initial_point_id, id_list)
        cell_connectivity = vtk_id_to_list(id_list)
        new_cell_ids = [cell_id for cell_id in cell_connectivity if
            cell_id not in old_cells]
        if len(new_cell_ids) > 1:
            raise ValueError('This should not happen')
        elif len(new_cell_ids) == 0:
            return connected_cell_points, old_cells

        # Add this cell and its points (in correct order).
        new_cell_id = new_cell_ids[0]
        old_cells.append(new_cell_id)
        new_cell_point_ids = vtk_id_to_list(
            pdi.GetCell(new_cell_id).GetPointIds())
        if new_cell_point_ids[0] == connected_cell_points[-1]:
            # First point of this cell is added to the last point of the last
            # cell.
            connected_cell_points.extend(new_cell_point_ids[1:])
            connected_cell_points, old_cells = add_cell_recursive(
                connected_cell_points, old_cells, new_cell_point_ids[-1])
        elif new_cell_point_ids[-1] == connected_cell_points[-1]:
            # Last point of this cell is added to the last point of the last
            # cell.
            new_cell_point_ids.reverse()
            connected_cell_points.extend(new_cell_point_ids[1:])
            connected_cell_points, old_cells = add_cell_recursive(
                connected_cell_points, old_cells, new_cell_point_ids[-1])
        elif new_cell_point_ids[0] == connected_cell_points[0]:
            # First point of this cell is added to the first point of the last
            # cell.
            new_cell_point_ids.reverse()
            connected_cell_points = (new_cell_point_ids[:-1] +
                connected_cell_points)
            connected_cell_points, old_cells = add_cell_recursive(
                connected_cell_points, old_cells, new_cell_point_ids[0])
        elif new_cell_point_ids[-1] == connected_cell_points[0]:
            # Last point of this cell is added to the first point of the last
            # cell.
            connected_cell_points = (new_cell_point_ids[:-1] +
                connected_cell_points)
            connected_cell_points, old_cells = add_cell_recursive(
                connected_cell_points, old_cells, new_cell_point_ids[0])
        else:
            raise ValueError('This should not happen')
        return connected_cell_points, old_cells

    id_list = vtk.vtkIdList()
    connected_cell_points = vtk_id_to_list(pdi.GetCell(cell_id).GetPointIds())
    old_cells = [cell_id]

    start_id = connected_cell_points[0]
    end_id = connected_cell_points[-1]
    connected_cell_points, old_cells = add_cell_recursive(
        connected_cell_points, old_cells, start_id)
    connected_cell_points, old_cells = add_cell_recursive(
        connected_cell_points, old_cells, end_id)

    return old_cells, connected_cell_points


# Start with the first poly line and search all poly lines connected to that
# line and so on. Then do the same with the first poly line that was not found
# and so on.
cell_ids = list(range(n_cells))
new_cells = []
while True:
    # Get the next available cell id.
    for i in cell_ids:
        if i is not None:
            next_id = i
            break
    else:
        break

    # Take the next available cell and look for all connected cells.
    old_cells, connected_cell_points = find_connected_cells(pdi, next_id)

    # Mark the found cell IDs since each one can only occur once in the
    # output.
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
