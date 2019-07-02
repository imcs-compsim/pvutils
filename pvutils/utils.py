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

    _dummy, extension = os.path.splitext(path)
    extension = extension.split('.')[-1].lower()

    if extension == 'pvd':
        return PVDReader(FileName=path)
    else:
        raise ValueError('Extension "{}" not defined!'.format(extension))


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
