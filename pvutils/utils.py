# -*- coding: utf-8 -*-
"""
Functions to simplify the use of scripts in paraview.
"""


# Import python modules.
import os
import sys
import numpy as np

# Import paraview module.
import paraview
import paraview.simple as pa


def _print_attibutes(obj, attributes, variable_name):
    """
    Print attributes of an object, given by a list of strings.
    """

    for att in attributes:
        attribute = getattr(obj, att)
        if isinstance(attribute, str):
            att_str = '\'{}\''.format(attribute)
            att_str = att_str.replace('\\', '\\\\')
        elif isinstance(attribute, float):
            att_str = '{0:.6g}'.format(attribute)
        elif type(attribute) == paraview.servermanager.VectorProperty:
            att_str_list = []
            for val in attribute:
                att_str = att_str_list.append('{0:.6g}'.format(val))
            att_str = '[{}]'.format(', '.join(att_str_list))
        else:
            att_str = str(attribute)
        print('{}.{} = {}'.format(variable_name, att, att_str))


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
        data = pa.PVDReader(FileName=path)
    elif extension == 'vtu':
        data = pa.XMLUnstructuredGridReader(FileName=[path])
    elif extension == 'exo':
        data = pa.ExodusIIReader(FileName=[path])
    elif extension == 'case':
        # In this case we get a multiblock structure, we do not want this.
        data = pa.EnSightReader(CaseFileName=path)
        data = pa.MergeBlocks(Input=data)
    else:
        raise ValueError('Extension "{}" not defined!'.format(extension))

    pa.RenameSource(base_name, data)
    return data


def display(data, line_width=None, line_color=None, solid_color=None,
        representation=None, nonlinear_subdividison=None, opacity=None):
    """
    Set the display options for the paraview object data.

    Colors are given as arrays with [R,G,B] values.
    """

    view = pa.GetActiveViewOrCreate('RenderView')
    display = pa.Show(data, view)

    if representation is not None:
        display.Representation = representation
    if line_color is not None:
        display.EdgeColor = line_color
    if line_width is not None:
        display.LineWidth = line_width
    if solid_color is not None:
        display.DiffuseColor = solid_color
        display.AmbientColor = solid_color
    if nonlinear_subdividison is not None:
        display.NonlinearSubdivisionLevel = nonlinear_subdividison
    if opacity is not None:
        display.Opacity = opacity


def contour(data, field='displacement', data_type='POINTS',
        vector_type='Magnitude'):
    """
    Set the contour options for a data item.

    Args
    ----
    data: ParaView data object
        ParaView item the warp filter will be applied to.
    field: String
        Name of field to be used for coloring/contouring.
    data_type: String
        Type/topology of data, encoded using ParaView-style namings.
    vector_type: String
        Component (if field is a vector field)
    """

    check_data(data, field, data_type=data_type)
    view = pa.GetActiveViewOrCreate('RenderView')
    display = pa.GetDisplayProperties(data, view=view)
    pa.ColorBy(display, (data_type, field, vector_type))


def warp(data, field='displacement', scale_factor=1.0):
    """
    Warp the data by a vector field.

    Args
    ----
    data: ParaView data object
        ParaView item the warp filter will be applied to.
    field: String
        Name of vector field to be used as "displacement" of the warp.
    scale_factor: Float
        Scaling factor.

    Return
    ------
    warp: ParaView data object
        The resulting ParaView item after applying the filter.
    """

    check_data(data, field, dimension=3)
    warp = pa.WarpByVector(Input=data)
    warp.Vectors = ['POINTS', field]
    warp.ScaleFactor = scale_factor
    return warp


def transform(data, translate=None, rotate=None, scale=None):
    """
    Apply a 'Transform' filter to a ParaView item.

    Args
    ----
    data: ParaView data object
        ParaView item the transpose filter will be applied to.
    transform: 3D-array
        Translation increment vector.
    rotate: 3D array
        Rotation increment vector.
    scale:
        Scaling in each Cartesian direction.

    Return
    ------
    transform: ParaView data object
        The resulting ParaView item after applying the filter.
    """

    if translate is None:
        translate = [0.0, 0.0, 0.0]
    if rotate is None:
        rotate = [0.0, 0.0, 0.0]
    if scale is None:
        scale = [1.0, 1.0, 1.0]

    transform = pa.Transform(Input=data)
    transform.Transform.Translate = translate
    transform.Transform.Rotate = rotate
    transform.Transform.Scale = scale
    return transform


def threshold(data, field='displacement', data_type='POINTS',
        threshold_range=None):
    """
    Apply a 'Threshold' filter to a ParaView item.

    Args
    ----
    data: ParaView data object
        ParaView item the threshold filter will be applied to.
    field: String
        Name of field whose values will be subject to thresholding.
    data_type: String
        Type/topology of data, encoded using ParaView-style namings.
    threshold_range: Array with 2 entries
        Array with min and max values of valid range of values.

    Return
    ------
    threshold: ParaView data object
        The resulting ParaView item after applying the filter.
    """

    if threshold_range is None:
        threshold_range = [-1.0e+12, 1.0e+12]

    threshold = pa.Threshold(Input=data)
    threshold.Scalars = [data_type, field]
    threshold.ThresholdRange = threshold_range
    return threshold


def check_data(data, name, data_type='POINTS', dimension=None,
        fail_on_error=True):
    """
    Check if data with the given name and dimension exists.

    Args
    ----
    data: ParaView data object
        ParaView item to be checked.
    name: String
        Name of field whose existence will be checked for.
    data_type: String
        Type of data, encoded using ParaView-style namings.
    dimension:
        Spatial dimension of requested data field.
    """

    visualization_data = pa.servermanager.Fetch(data)

    if data_type == 'POINTS':
        data = visualization_data.GetPointData()
    elif data_type == 'CELLS':
        data = visualization_data.GetCellData()
    else:
        raise ValueError(('Can\'t access data of type {} so far. '
            + 'Needs to be implemented.').format(data_type))
    field_data = data.GetArray(name)

    if field_data is None:
        if fail_on_error:
            names = [data.GetArrayName(i)
                for i in range(data.GetNumberOfArrays())]
            raise ValueError(('Could not find {} data with the name {}! '
                + 'Available names: {}').format(data_type, name, names))
        return False

    if dimension is not None:
        if not field_data.GetNumberOfComponents() == dimension:
            if fail_on_error:
                raise ValueError(
                    'The field {} has {} instead of {} dimensions!'.format(
                        name, field_data.GetNumberOfComponents(), dimension)
                    )
            return False

    return True


def tube(data, slices=8):
    """
    Apply the tube filter to a beam.
    """

    check_data(data, 'cross_section_radius', dimension=1)
    tube = pa.Tube(Input=data)
    tube.Scalars = [None, 'cross_section_radius']
    tube.VaryRadius = 'By Absolute Scalar'
    tube.NumberofSides = slices
    return tube


def get_base(data):
    """
    Return the root item of a given item, e.g. the base geometry object.
    """

    if 'Input' in dir(data):
        return get_base(data.Input)
    else:
        return data


def reset_paraview():
    """
    Delete all data in paraview.
    https://stackoverflow.com/questions/48580653/paraview-programmatically-reset-session
    """

    pa.Disconnect()
    pa.Connect()


def programmable_filter(source, name):
    """
    Apply a programmable filter from this git repository.
    """

    filter_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        'filters',
        '{}.py'.format(name))

    pv_filter = pa.ProgrammableFilter(Input=source)
    pv_filter.Script = 'execfile("{}")'.format(filter_path)
    return pv_filter


def setup_view(view, *args, **kwargs):
    """
    Allow the user to setup and return the relevant values.
    """

    # When python3 is used this construct will be obsolete.
    if sys.version_info >= (3, 0):
        raise ValueError('The keyword management in setup_view should be '
            + 'adapted to python3.')

    # Default keyword arguments.
    kwargs_default = {
        # If the size of the view should be fixed, i.e. preview mode should be
        # used. If pvpython is used, this does not have an effect. The size
        # will be taken from view.
        'fixed_size': False,
        }

    # Set the keyword arguments.
    for key in kwargs_default.keys():
        if key in kwargs:
            value = kwargs[key]
            del(kwargs[key])
        else:
            value = kwargs_default[key]
        exec(key + ' = value')
    # Check that no keyword arguments remain.
    if len(kwargs) > 0:
        raise ValueError('Unsupported keyword arguments {} given.'.format(
            kwargs.keys()))

    # Check which paraview interpreter is used and setup the view accordingly.
    pa.Render(view)
    if is_pvpython():
        # Stop for user modifications.
        pa.Interact(view)
    else:
        # Get the desired aspect ratio.
        if fixed_size:
            # Enter preview mode.
            layout = pa.GetLayout()
            layout.PreviewMode = view.ViewSize

    print_view_state(view, *args)


def print_view_state(view, *args):
    """
    Print the relevant view state information, so the user can set the view in
    the GUI and then copy the state to a script.
    """

    # Display the view attributes.
    attributes = [
        'CameraPosition',
        'CameraFocalPoint',
        'CameraViewUp',
        'CameraViewAngle',
        'CameraParallelScale',
        'OrientationAxesVisibility',
        'CameraParallelProjection',
        'ViewSize'
        ]
    _print_attibutes(view, attributes, 'view')
    print('')

    # If additional items are given to this function, print their properties.
    for arg in args:

        item_string = str(arg)
        if 'ScalarBarWidgetRepresentation' in item_string:
            # Display the color bar attributes.
            attributes = [
                'Title',
                'ComponentTitle',
                'WindowLocation',
                'Orientation',
                'ScalarBarLength',
                'ScalarBarThickness',
                'Position'
                ]
            print('')
            _print_attibutes(arg, attributes, 'color_bar')


def get_size_pixel(size, dpi):
    """
    Convert a output size in cm to a size in pixel.
    """

    inch = 2.54
    return (np.array(size) / inch * dpi).astype(int)


def is_pvpython():
    """
    Check if the current script is executed by paraview or pvpthon.
    """
    return not sys.argv[0] == ''


def set_colorbar_font(color_bar, font_size, dpi, font=None):
    """
    Set the font options for the a color bar.
    This only makes sense if the screenshots are exported with the option:
        FontScaling='Do not scale fonts'
    """

    # 72 points are in an inch. Calculate the needed pixels for the font size.
    dpi_font = 72.
    font_size_pixel = (np.array(font_size) / dpi_font * dpi).astype(int)

    color_bar.TitleFontSize = font_size_pixel[0]
    color_bar.LabelFontSize = font_size_pixel[1]

    if font == 'TeX':
        dirname = os.path.dirname(__file__)
        font_file = os.path.join(dirname, '..', 'utilities',
            'latin-modern-regular.ttf')
        color_bar.TitleFontFamily = 'File'
        color_bar.TitleFontFile = font_file
        color_bar.LabelFontFamily = 'File'
        color_bar.LabelFontFile = font_file
