# -*- coding: utf-8 -*-
"""
Functions to simplify the use of scripts in paraview.
"""


# Import python modules.
import os

# Import paraview module.
from paraview.simple import *


def load_file(path):
    """
    Load a file according to its extension.
    """

    if not os.path.isfile(path):
        raise ValueError('The file "{}" does not exist!'.format(path))

    split_path = os.path.basename(path).split('.')
    extension = split_path[-1].lower()
    base_name = ''.join(split_path[:-1])

    if extension == 'pvd':
        data = PVDReader(FileName=path)
    elif extension == 'exo':
        data = ExodusIIReader(FileName=[path])
    else:
        raise ValueError('Extension "{}" not defined!'.format(extension))

    RenameSource(base_name, data)
    return data


def display(data, line_width=None, line_color=None, solid_color=None,
        representation=None, nonlinear_subdividison=None):
    """
    Set the display options for the paraview object data.
    """

    view = GetActiveViewOrCreate('RenderView')
    display = Show(data, view)

    if representation is not None:
        display.Representation = representation
    if line_color is not None:
        display.EdgeColor = line_color
    if line_width is not None:
        display.LineWidth = line_width
    if solid_color is not None:
        display.DiffuseColor = solid_color
    if nonlinear_subdividison is not None:
        display.NonlinearSubdivisionLevel = nonlinear_subdividison


def warp(data, field='displacement', scale_factor=1):
    """
    Wrap the data by the displacement vector.
    """

    check_data(data, field)
    warp = WarpByVector(Input=data)
    warp.Vectors = ['POINTS', field]
    warp.ScaleFactor = scale_factor
    return warp


def check_data(data, name, data_type='point', dimension=3,
        fail_on_error=True):
    """
    Check if data with the given name and dimension exists.
    """

    vtk_data = servermanager.Fetch(data)
    point_data = vtk_data.GetPointData()
    field_data = point_data.GetArray(name)

    if field_data is None:
        if fail_on_error:
            names = [point_data.GetArrayName(i)
                for i in range(point_data.GetNumberOfArrays())]
            raise ValueError(('Could not find {} data with the name {}! '
                + 'Available names: {}').format(data_type, name, names))
        return False

    if not field_data.GetNumberOfComponents() == dimension:
        if fail_on_error:
            raise ValueError(
                'The field {} has {} instead of {} dimensions!'.format(
                    name, field_data.GetNumberOfComponents(), dimension)
                )
        return False

    return True


def beam_tube(data):
    """
    Apply the tube filter to a beam.
    """
    pass
