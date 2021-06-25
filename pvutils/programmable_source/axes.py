"""
This script generates a point with vector data. Can be used to visualize
coordinate axes.

Use with programmable source filter in paraview and execute with:
>> execfile('path to file')
"""

from paraview import vtk


def get_keyword_argument(name, kwargs_id):
    """
    Return the wanted keyword argument.
    """
    import paraview
    if hasattr(paraview, 'programmable_source_kwargs'):
        if name in paraview.programmable_source_kwargs[kwargs_id].keys():
            return paraview.programmable_source_kwargs[kwargs_id][name]
    raise ValueError('Keyword argument "{}" not found'.format(name))


# Get keyword arguments.
origin = get_keyword_argument('origin', kwargs_id)
basis = get_keyword_argument('basis', kwargs_id)

# Add the origin.
output = self.GetOutput()
newPts = vtk.vtkPoints()
newPts.InsertPoint(0, origin[0], origin[1], origin[2])
output.SetPoints(newPts)

# Add the basis vectors.
for i, base in enumerate(basis):
    base_data = vtk.vtkFloatArray();
    base_data.SetNumberOfComponents(3);
    base_data.SetName('base_{}'.format(i + 1));
    base_data.InsertNextTuple3(base[0], base[1], base[2]);
    output.GetPointData().AddArray(base_data)
