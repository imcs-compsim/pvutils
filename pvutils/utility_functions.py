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
import pvutils
from vtk.util import numpy_support as VN


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

    reset_function_arguments()
    pa.Disconnect()
    pa.Connect()


def reset_layout():
    """
    Closes the current layout and opens a new clean one. All the sources stay
    loaded.
    """

    # Get the old layout and its name.
    old_layout = pa.GetLayout()
    old_view = get_view()
    layouts = pa.GetLayouts()
    for key in layouts.keys():
        name, _id = key
        if layouts[key] == old_layout:
            old_name = name
            break
    else:
        raise ValueError('Could not get the layout name!')

    # Create the new layout.
    new_layout = pa.CreateLayout('temp layout')
    new_view = pa.CreateView('RenderView')
    new_layout.AssignView(0, new_view)

    # Delete the old stuff.
    pa.Delete(old_view)
    del old_view
    pa.RemoveLayout(old_layout)

    # Rename the new layout.
    pa.RenameLayout(old_name, new_layout)


def programmable_filter(source, name, **kwargs):
    """
    Apply a programmable filter from this git repository.
    """

    filter_path = os.path.join(
        os.path.dirname(__file__),
        'programmable_filters',
        '{}.py'.format(name))

    # Store the kwargs in a variable in the namespace of the ParaView python
    # module. This variable can be retrieved by the programmable filter. The
    # kwargs_id stores which keyword arguments belong to this programmable
    # filter.
    if hasattr(paraview, 'programmable_filter_kwargs'):
        paraview.programmable_filter_kwargs.append(kwargs)
    else:
        paraview.programmable_filter_kwargs = [kwargs]
    kwargs_id = len(paraview.programmable_filter_kwargs) - 1

    pv_filter = pa.ProgrammableFilter(Input=source)
    pv_filter.Script = (
        'kwargs_id = {}\n'.format(kwargs_id) +
        'execfile("{}")'.format(filter_path))
    pa.RenameSource(name, pv_filter)
    return pv_filter


def programmable_source(name, output_data='vtkUnstructuredGrid', **kwargs):
    """
    Apply a programmable source filter from this git repository.
    """

    filter_path = os.path.join(
        os.path.dirname(__file__),
        'programmable_source',
        '{}.py'.format(name))

    # Store the kwargs in a variable in the namespace of the ParaView python
    # module. This variable can be retrieved by the programmable source. The
    # kwargs_id stores which keyword arguments belong to this programmable
    # source.
    if hasattr(paraview, 'programmable_source_kwargs'):
        paraview.programmable_source_kwargs.append(kwargs)
    else:
        paraview.programmable_source_kwargs = [kwargs]
    kwargs_id = len(paraview.programmable_source_kwargs) - 1

    pv_source = pa.ProgrammableSource()
    pv_source.OutputDataSetType = output_data
    pv_source.Script = (
        'kwargs_id = {}\n'.format(kwargs_id) +
        'execfile("{}")'.format(filter_path)
        )
    pa.RenameSource(name, pv_source)
    return pv_source


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
        # Name of the view object in the print out. You can add spaces here to
        # directly copy the resulting print output to the code.
        'view_name': None,
        # Name of the variable that holds the view size in pixel.
        'size_pixel_name': None,
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

    if view_name is not None:
        set_function_arguments('print_view_state', 'view_name', view_name)
    if size_pixel_name is not None:
        set_function_arguments('print_view_state', 'size_pixel_name',
            size_pixel_name)
    print_view_state(view, *args)


def print_view_state(view, *args):
    """
    Print the relevant view state information, so the user can set the view in
    the GUI and then copy the state to a script.
    """

    view_name = get_function_arguments('view_name', default_value='view')
    size_pixel_name = get_function_arguments('size_pixel_name',
        default_value='size_pixel')

    # Display the view attributes.
    attributes = [
        'CameraPosition',
        'CameraFocalPoint',
        'CameraViewUp',
        'CameraViewAngle',
        'CameraParallelScale',
        'OrientationAxesVisibility',
        'CameraParallelProjection',
        'InteractionMode'
        ]
    print('{} = {}\n\n'.format(size_pixel_name, view.ViewSize))
    print('{}.ViewSize = {}'.format(view_name, size_pixel_name))
    _print_attibutes(view, attributes, view_name)
    print('')

    # If additional items are given to this function, print their properties.
    args = list(args)
    extra_args = get_function_arguments('extra_args', default_value=[])
    args.extend(extra_args)
    for name, arg in args:

        item_string = str(arg)
        if 'ScalarBarWidgetRepresentation' in item_string:
            # Display the color bar attributes.
            attributes = [
                'Title',
                'ComponentTitle',
                'WindowLocation',
                'Orientation',
                'HorizontalTitle',
                'ScalarBarLength',
                'ScalarBarThickness',
                'Position',
                'UseCustomLabels',
                'AddRangeLabels',
                'DrawAnnotations',
                'AddRangeAnnotations',
                'AutomaticAnnotations',
                'DrawNanAnnotation',
                'DrawTickLabels'
                ]
            print('')
            _print_attibutes(arg, attributes, name)

        else:
            raise TypeError('Got unexpected type')


def set_print_view_state_color_bar(color_bar, name='color_bar'):
    """
    Set a scalar bar that will be included in print_view_state.
    """

    extra_args = get_function_arguments('extra_args',
        function_name='print_view_state')
    if extra_args is None:
        extra_args = []
        set_function_arguments('print_view_state', 'extra_args', extra_args)

    extra_args.append([name, color_bar])


def reset_print_view_state_color_bar():
    """
    Reset optional arguments for print_view_state.
    """

    delete_function_arguments('print_view_state', 'extra_args')


def delete_function_arguments(function_name, variable_name):
    """
    Remove an argument for a function that will be called from within
    ParaView.
    """

    if hasattr(paraview, 'pvutils_args'):
        if function_name in paraview.pvutils_args.keys():
            variable_dictionary = paraview.pvutils_args[function_name]
            if variable_name in variable_dictionary.keys():
                del variable_dictionary[variable_name]


def set_function_arguments(function_name, variable_name, variable_value,
        overwrite=False):
    """
    Set arguments for a function that will be called from within ParaView.

    Args
    ----
    function_name: str
        Name of the function that will use this argument.
    variable_name: str
        Name of the variable.
    variable_value:
        Name of the value for the variable.
    overwrite: bool
        If the variable should be overwritten if it is already set. Otherwise
        an error will be thrown it an existing variable is overwritten.
    """

    if not hasattr(paraview, 'pvutils_args'):
        paraview.pvutils_args = {}

    if function_name not in paraview.pvutils_args.keys():
        paraview.pvutils_args[function_name] = {}

    variable_dictionary = paraview.pvutils_args[function_name]
    if variable_name in variable_dictionary.keys() and not overwrite:
        raise KeyError(('The variable {} in the function {} is ' +
            'already set!').format(variable_name, function_name))

    variable_dictionary[variable_name] = variable_value


def reset_function_arguments():
    """
    Reset all function arguments for functions called from within ParaView.
    """
    if hasattr(paraview, 'pvutils_args'):
        paraview.pvutils_args = {}


def get_function_arguments(variable_name, default_value=None,
        function_name=None):
    """
    Get an argument that has previously been set. If it has not been set, the
    default value is returned.

    Args
    ----
    variable_name: str
        Name of the variable.
    default_value:
        If the variable has not been set, this variable is returned.
    function_name: str
        Name of the function that uses this variable. If this argument is not
        given, the name of the function will be determined from the stack.
    """

    if function_name is None:
        import inspect
        stack = inspect.stack()
        function_name = stack[1][3]

    if hasattr(paraview, 'pvutils_args'):
        if function_name in paraview.pvutils_args.keys():
            if variable_name in paraview.pvutils_args[function_name].keys():
                return paraview.pvutils_args[function_name][variable_name]

    return default_value


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


def update_scene():
    """
    This has the same effect, as if 'Apply' is clicked in ParaView, i.e. the
    time steps are loaded correctly in the GUI.
    """

    # Sometimes these commands are used, but it seems it it not needed. Keep
    # this comment for reference.
    # view = pvutils.get_view()
    # view.Update()
    # view.ResetCamera()
    # pa.UpdatePipeline()
    # pa.Render()

    animation_scene = pa.GetAnimationScene()
    animation_scene.UpdateAnimationUsingDataTimeSteps()
    return animation_scene


def set_color_range(field, val_min, val_max):
    """
    Set the max and min values for the contour plot color range.

    Args
    ----
    field: str
        Name of field quantity.
    val_min, val_max: float (or int)
        Minimum and maximum values for the color range
    """

    color_transfer_function = pa.GetColorTransferFunction(field)
    color_transfer_function.RescaleTransferFunction(float(val_min),
        float(val_max))
    opacity_transfer_function = pa.GetOpacityTransferFunction(field)
    opacity_transfer_function.RescaleTransferFunction(float(val_min),
        float(val_max))


def von_mises_stress(source, field_name='nodal_cauchy_stresses_xyz',
        field_type='POINTS'):
    """
    Add a field with the vonMises stress.

    Args
    ----
    source: pa.Item
        ParaView item that the filter should be applied to.
    field_name: str
        Name of the field. Should be a second-order tensor field.
    field_type: str:
        Type of the input field.
    """

    # Check if arrays are available.
    check_data(source, field_name, data_type=field_type)

    function = (
        'sqrt(' +
        '{0}_XX^2 + {0}_YY^2 + {0}_ZZ^2' +
        '- {0}_XX * {0}_YY - {0}_XX * {0}_ZZ - {0}_ZZ * {0}_YY' +
        '+ 6*({0}_XY^2 + {0}_YZ^2 + {0}_XZ^2))')

    calc = pa.Calculator(Input=source)
    calc.ResultArrayName = 'von_mises_stress_{}'.format(field_type)
    calc.Function = function.format(field_name)
    if field_type == 'POINTS':
        calc.AttributeType = 'Point Data'
    else:
        calc.AttributeType = 'Cell Data'
    return calc


def add_coordinate_axes(origin=None, basis=None, scale=1.0, resolution=20,
        show=True):
    """
    Add arrow representations for coordinate axes.

    Args
    ----
    origin: [float]
        Point where the axes will be drawn.
    basis: [[float]]
        List of basis vectors.
    scale: float
        Scale factor for displaying basis vectors.
    resolution: int
        ParaView internal resolution for arrow representation.
    show: bool
        If the basis vectors should be displayed.

    Return
    ----
    {
        'axis_source': ParaView data
            Source data for the bais vectors,
        'base_glyphs': ParaView data
            Axes visualization data.
    }
    """

    # Set default values for origin and basis vectors.
    if origin is None:
        origin = [0, 0, 0]
    if basis is None:
        basis = [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
            ]
    sorce = programmable_source('axes', origin=origin, basis=basis)
    base_glyphs = []
    for i, base in enumerate(basis):
        base_glyph = pvutils.glyph(sorce)
        base_glyph.GlyphType = 'Arrow'
        base_glyph.OrientationArray = ['POINTS', 'base_{}'.format(i + 1)]
        base_glyph.ScaleArray = ['POINTS', 'base_{}'.format(i + 1)]
        base_glyph.ScaleFactor = scale
        base_glyph.GlyphMode = 'All Points'
        base_glyph.GlyphType.TipResolution = resolution
        base_glyph.GlyphType.ShaftResolution = resolution

        # The arrow representation is scaled, so the arrow does not look
        # "thick" when getting longer.
        length = np.linalg.norm(base)
        base_glyph.GlyphType.TipRadius = 0.1 / length
        base_glyph.GlyphType.TipLength = 0.35 / length
        base_glyph.GlyphType.ShaftRadius = 0.03 / length

        if show:
            pa.Show(base_glyph)
        base_glyphs.append(base_glyph)

    return {
        'axis_source': sorce,
        'base_glyphs': base_glyphs
        }


def get_bounding_box(source):
    """
    Get the bounding box dimensions (using the outline filter).

    Return
    ----
    [[min_x, max_x], [min_y, max_y], [min_z, max_z]]
    """

    outline = pa.Outline(Input=source)
    vtk_data = get_vtk_data_as_numpy(outline, coordinates=True)
    position = vtk_data['coordinates']
    max_min_coordinates = [
        [np.min(position[:, i]), np.max(position[:, i])]
        for i in range(3)]
    pa.Delete(outline)
    return max_min_coordinates


def get_vtk_data_as_numpy(source, coordinates=False, point_data=False,
        cell_data=False, cell_connectivity=False):
    """
    Return all vtk data arrays.

    Args
    ----
    source:
        Paraview source object.
    coordinates:
        If point coordinates should be returned.
    point_data:
        If point data should be returned.
    cell_data:
        If cell data should be returned.
    cell_connectivity:
        If cell connectivity and cell types should be returned.
    Return
    ------
    {
    'coordinates': coordinates,
    'point_data': point_data,
    'cell_data': cell_data,
    'cell_types': cell_types,
    'cell_connectivity': cell_connectivity
    }
    """

    data = pa.servermanager.Fetch(source)

    if coordinates:
        coordinates = VN.vtk_to_numpy(data.GetPoints().GetData())
    else:
        coordinates = None

    def get_data_array(input_data):
        """
        Extract data from an array.
        """
        data_dir = {}
        for i in range(input_data.GetNumberOfArrays()):
            data_dir[input_data.GetArrayName(i)] = VN.vtk_to_numpy(
                input_data.GetArray(i))
        return data_dir

    if point_data:
        point_data = get_data_array(data.GetPointData())
    else:
        point_data = None
    if cell_data:
        cell_data = get_data_array(data.GetCellData())
    else:
        cell_data = None

    def vtk_id_to_list(vtk_id_list):
        """
        Convert an vtk id list to a python list.
        """
        return [int(vtk_id_list.GetId(i_id)) for i_id in
            range(vtk_id_list.GetNumberOfIds())]

    if cell_connectivity:
        n_cells = data.GetNumberOfCells()
        cell_types = [None] * n_cells
        cell_connectivity = [None] * n_cells
        for i in range(n_cells):
            cell_types[i] = data.GetCellType(i)
            cell_connectivity[i] = vtk_id_to_list(data.GetCell(i).GetPointIds())
    else:
        cell_types = None
        cell_connectivity = None

    return {
        'coordinates': coordinates,
        'point_data': point_data,
        'cell_data': cell_data,
        'cell_types': cell_types,
        'cell_connectivity': cell_connectivity
        }


def list_to_mathematica_string(data, name=None, string_format='{:.15f}'):
    """
    Convert a list to a Mathematica string.

    Args
    ----
    data: array-like
        The data that should be converted to a string.
    name: str
        Name of the Mathematica variable. If this argument is given the
        output is in the form: "name = {data};"
    string_format: str
        Format for numeric values.
    """

    def data_to_string(data):
        """Recursive function to convert array-like structures to a string."""
        if isinstance(data, list) or isinstance(data, np.ndarray):
            data_string = '{'
            for var in data:
                data_string += data_to_string(var)
            data_string = data_string[:-2] + '}, '
            return data_string
        else:
            return string_format.format(data) + ', '

    data_string = data_to_string(data)[:-2]
    if name is None:
        return data_string
    else:
        return '{} = {};'.format(name, data_string)


def export_to_tikz(name, view=None, dpi=300, color_transfer_functions=None,
        number_format='{$\\pgfmathprintnumber[sci,precision=1,sci generic={mantissa sep=,exponent={\\mathrm{e}{##1}}}]{\\tick}$}'):
    """
    Export a screenshot and wrap the color bars inside a TikZ axis.

    Args
    ----
    name: str
        Name of the output files. Can include directory paths.
    view: pa.View
        Active ParaView view.
    dpi: int
        Desired resolution of image in dots per inch.
    color_transfer_functions: [color_transfer_function]
        A list of the color transfer functions that should be wrapped by TikZ.
    number_format: str
        Number format to be used for the TeX ticks.
    """

    def rel_to_tikz(val, direction):
        """
        Convert a horizontal positioning to TikZ (cm).
        """
        return view.ViewSize[direction] * val / dpi * 2.54

    def dots_to_tikz(val):
        """
        Convert a distance in pt to cm.
        """
        return float(val) / dpi * 2.54

    def get_min_max_values(color_transfer_function):
        """
        Return the min and maximum value for a given color transfer function.
        """

        points = color_transfer_function.RGBPoints
        return points[0], points[-4]

    if view is None:
        view = get_view()

    # Name of image and TikZ file.
    image_name_full = name + '.png'
    image_name = os.path.split(image_name_full)[-1]
    tikz_name_full = name + '.tex'

    # Set options for colorbars.
    draw_tick_marks_old = []
    draw_tick_labels_old = []
    add_range_labels_old = []
    title_old = []
    if color_transfer_functions is not None:
        for color_transfer_function in color_transfer_functions:
            color_bar = pa.GetScalarBar(color_transfer_function, view)

            draw_tick_marks_old.append(color_bar.DrawTickMarks)
            draw_tick_labels_old.append(color_bar.DrawTickLabels)
            add_range_labels_old.append(color_bar.AddRangeLabels)
            title_old.append(color_bar.Title)
            color_bar.DrawTickMarks = 0
            color_bar.DrawTickLabels = 0
            color_bar.AddRangeLabels = 0
            color_bar.Title = ''

    # Save the screenshot.
    pa.SaveScreenshot(
        image_name_full,
        view,
        ImageResolution=view.ViewSize,
        OverrideColorPalette='WhiteBackground',
        TransparentBackground=0,
        FontScaling='Do not scale fonts'
        )

    # Create the TikZ code
    tikz_code = '''%% This file was created with pvutils
%% Use the following includes in the LaTeX header:
%\\usepackage{{tikz}}
%\\usepackage{{pgfplots}}
%% Optional:
%\\pgfplotsset{{compat=1.16}}
\\begin{{tikzpicture}}
\\node[anchor=south west,inner sep=0] (image) at (0,0) {{\\includegraphics[scale={scale}]{{{image_name}}}}};\n'''.format(
        scale=72.0 / dpi,
        image_name=image_name)

    if color_transfer_functions is not None:
        for i, color_transfer_function in enumerate(color_transfer_functions):
            color_bar = pa.GetScalarBar(color_transfer_function, view)

            rel_pos = color_bar.Position
            min_max = get_min_max_values(color_transfer_function)

            # Get the ticks.
            tick = []
            if add_range_labels_old[i] == 1:
                tick.extend(min_max)
            tick.extend(color_bar.CustomLabels)
            tick.sort()
            tick_str = ','.join(map(str, tick))

            # Add the code that is valid for all types of labels.
            tikz_code += '''\\begin{{axis}}[
scale only axis,
at={{({pos[0]}cm,{pos[1]}cm)}},
tick label style={{font=\\footnotesize}},
title={title},
xticklabel={number_format},
yticklabel={number_format},
ymin={min_max[0]},
ymax={min_max[1]},
xmin={min_max[0]},
xmax={min_max[1]},\n'''.format(
                pos=[rel_to_tikz(rel_pos[j], j) for j in range(2)],
                min_max=min_max,
                title=title_old[i],
                number_format=number_format
                )

            if color_bar.Orientation == 'Horizontal':
                tikz_code += '''ytick=\empty,
height={height}cm,
width={width}cm,
xtick={{{tick}}},
xtick pos=right,
xtick align=outside,
title style={{yshift=10pt,}},\n'''.format(
                    height=dots_to_tikz(color_bar.ScalarBarThickness),
                    width=rel_to_tikz(color_bar.ScalarBarLength, 0),
                    tick=tick_str
                    )

            else:
                tikz_code += '''xtick=\empty,
height={height}cm,
width={width}cm,
ytick={{{tick}}},
ytick pos=right,
ytick align=outside,\n'''.format(
                    width=dots_to_tikz(color_bar.ScalarBarThickness),
                    height=rel_to_tikz(color_bar.ScalarBarLength, 1),
                    tick=tick_str
                    )

            tikz_code += ']\n\end{axis}\n'

            # Reset the initial values for the color bar.
            color_bar.DrawTickMarks = draw_tick_marks_old[i]
            color_bar.DrawTickLabels = draw_tick_labels_old[i]
            color_bar.AddRangeLabels = add_range_labels_old[i]
            color_bar.Title = title_old[i]

    tikz_code += '\end{tikzpicture}%'

    # Write TikZ code to file.
    with open(tikz_name_full, 'w') as text_file:
        text_file.write(tikz_code)
