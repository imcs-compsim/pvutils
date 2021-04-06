# -*- coding: utf-8 -*-
"""
Functions to simplify the use of scripts in ParaView.
"""


# Import python modules.
import os
import sys
import numpy as np

# Import ParaView module.
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
    elif extension == 'pvtu':
        data = pa.XMLPartitionedUnstructuredGridReader(FileName=[path])
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
    Set the display options for the ParaView object data.

    Colors are given as arrays with [R,G,B] values.
    """

    view = get_view()
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
    display = get_display(data)
    pa.ColorBy(display, (data_type, field, vector_type))


def get_source_name(item):
    """
    Return the name of a source item in the ParaView tree. The name does not
    have to be unique.
    """

    # Get all sources in the session.
    sources = pa.GetSources()
    for key, value in sources.items():
        if value == item:
            return key[0]
    else:
        raise ValueError('The item "{}" was not found in sources!'.format(
            item))


def get_field_names(item):
    """
    Return a dictionary with the available field data names in data.

    Return
    ----
    names: dict
        The return dictionary has the following structure:
        {
            'FIELD': [
                ('field_name_1', number_of_components),
                ('field_name_2', number_of_components)
                ],
            'CELLS': [
                ('cell_field_name_1', number_of_components),
                ('cell_field_name_2', number_of_components)
                ],
            'POINTS': [
                ('point_field_name_1', number_of_components),
                ('point_field_name_2', number_of_components)
                ]
        }
    """

    # Loop over global, cell and point data to get the available names.
    return_dict = {}
    visualization_data = pa.servermanager.Fetch(item)
    for funct, data_type in [
            (visualization_data.GetFieldData, 'FIELD'),
            (visualization_data.GetCellData, 'CELLS'),
            (visualization_data.GetPointData, 'POINTS')
            ]:

        # Create the entry in the return dictionary.
        return_dict[data_type] = []

        # Get all names for this data type.
        data = funct()
        names = [data.GetArrayName(i)
                 for i in range(data.GetNumberOfArrays())]

        # Add the number of components for the fields.
        for name in names:
            field_data = data.GetArray(name)
            return_dict[data_type].append(
                (name, field_data.GetNumberOfComponents())
                )

    return return_dict


def check_data(item, name, data_type='POINTS', dimension=None,
        fail_on_error=True):
    """
    Check if data with the given name and dimension exists.

    Args
    ----
    item: ParaView data object
        ParaView item to be checked.
    name: String
        Name of field whose existence will be checked for.
    data_type: String
        Type of data, encoded using ParaView-style namings.
    dimension:
        Spatial dimension of requested data field.
    """

    # Check if the field exists with the correct dimension.
    field_names = get_field_names(item)[data_type]
    for field_name, field_dimension in field_names:
        if field_name == name:
            if (dimension is not None) and (not field_dimension == dimension):
                if fail_on_error:
                    raise ValueError((
                        'The field {} has {} instead of {} dimensions in the '
                        + 'source "{}"!').format(
                            name, field_dimension, dimension,
                            get_source_name(item))
                        )
                return False
            return True
    else:
        # No match was found in the field names.
        if fail_on_error:
            raise ValueError(('Could not find {} data with the name {} in the'
                + 'source "{}"! Available names: {}').format(data_type, name,
                    get_source_name(item), field_names))
        return False


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
    Delete all data in ParaView.
    https://stackoverflow.com/questions/48580653/paraview-programmatically-reset-session
    """

    pa.Disconnect()
    pa.Connect()


def programmable_filter(source, name, **kwargs):
    """
    Apply a programmable filter from this git repository.
    """

    filter_path = os.path.join(
        os.path.dirname(__file__),
        'programmable_filters',
        '{}.py'.format(name))

    # Store the kwargs in a variable in the namespace of the ParaView python
    # module.
    paraview.programmable_filter_kwargs = kwargs

    pv_filter = pa.ProgrammableFilter(Input=source)
    pv_filter.Script = 'execfile("{}")'.format(filter_path)
    return pv_filter


def get_display(data, view=None):
    """Return the display object from ParaView."""
    if view is None:
        view = get_view()
    return pa.GetDisplayProperties(data, view=view)


def get_view():
    """Return the view object from ParaView."""
    return pa.GetActiveViewOrCreate('RenderView')


def setup_view(*args, **kwargs):
    """
    Allow the user to setup and return the relevant values.
    """

    # When python3 is used this construct will be obsolete.
    if sys.version_info >= (3, 0):
        raise ValueError('The keyword management in setup_view should be '
            + 'adapted to python3.')

    # Default keyword arguments.
    kwargs_default = {
        # The current view object. If none is given, the default one is taken.
        'view': None,
        # If the size of the view should be fixed, i.e. preview mode should be
        # used. If pvpython is used, this does not have an effect. The size
        # will be taken from view.
        'fixed_size': False
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

    # Get the view object.
    if view is None:
        view = get_view()

    # Check which ParaView interpreter is used and setup the view accordingly.
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
        'ViewSize',
        'InteractionMode'
        ]
    _print_attibutes(view, attributes, 'view')
    print('')

    # If additional items are given to this function, print their properties.
    args = list(args)
    if hasattr(paraview, 'print_view_state_scalar_bar'):
        args.extend(paraview.print_view_state_scalar_bar)
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


def print_view_state_set_scalar_bar(*args):
    """
    Set a scalar bar that will be included in print_view_state.
    """

    paraview.print_view_state_scalar_bar = args


def get_size_pixel(size, dpi):
    """
    Convert a output size in cm to a size in pixel.
    """

    inch = 2.54
    return (np.array(size) / inch * dpi).astype(int)


def is_pvpython():
    """
    Check if the current script is executed by ParaView or pvpthon.
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


def get_available_timesteps():
    """
    Return a list with all available time steps in the current session.
    """

    scene = pa.GetAnimationScene()
    times = scene.TimeKeeper.TimestepValues

    # In the case of a single time step the function
    # scene.TimeKeeper.TimestepValues returns a float instead of a list with
    # one entry. Here we check that the return value of this function is
    # always a list.
    if isinstance(times, float):
        return [times]
    else:
        return times


def set_timestep(time, fail_on_not_available_time=True):
    """
    Set the time step in the current session.
    This works fine, BUT ParaView will use a wrong time (first available time),
    even tough a different time step is displayed.

    Args
    ----
    time: scalar
        Time that should be set.
    fail_on_not_available_time: bool
        If this is true and the given time does not exist in ParaView an error
        will be thrown.
    """

    if fail_on_not_available_time:
        # Check that the given time is a time step in the current state.
        times = get_available_timesteps()
        if min(np.abs(np.array(times) - time)) > 1e-10:
            raise ValueError(('The given time {} is not a valid time step in '
                + 'ParaView. The valid times are: {}').format(time, times))

    scene = pa.GetAnimationScene()
    scene.TimeKeeper.Time = time

    # This render command is disabled as it resulted in failing pipelines in
    # GitLab.
    # pa.Render()


def get_parents(source):
    """
    Return a list with all direct parents of an object in ParaView, including
    the object itself.
    """

    parents = [source]
    # Only add the parent if the object has a method 'get_parents'.
    if 'Input' in dir(source):
        # In some cases input returns a list, e.g. for programmable filters. In
        # this case this function only works if there is exactly one element in
        # that list.
        pv_input = source.Input
        is_pv_list = ('paraview.servermanager.InputProperty'
            in str(type(pv_input)))
        if is_pv_list and len(pv_input) > 1:
            raise ValueError(('"Input" for the object "{}" returned the '
                + 'following list with more than 1 entry: {}'
                ).format(source, pv_input))
        elif is_pv_list:
            # For some reason we can not access the __getitem__ method with
            # pv_input[index]. Therefore, we call the method directly.
            pv_input = pa.servermanager.InputProperty.__getitem__(pv_input, 0)
        parents.extend(get_parents(pv_input))
    return parents