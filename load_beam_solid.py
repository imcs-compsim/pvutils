# -*- coding: utf-8 -*-
"""
Load a solid and a beam into paraview and apply some filters.
"""

# Import the simple module from the paraview.
from paraview.simple import GetActiveViewOrCreate, PVDReader, RenameSource, \
    WarpByVector, Show, _DisableFirstRenderCameraReset, ColorBy, Glyph, \
    ExtractSurface, Tube, Calculator, CellDatatoPointData, Hide, \
    XMLUnstructuredGridReader, ProgrammableFilter


# Python imports.
import os


def create_glyph(base, glyph_type='Arrow', scalars='None', vectors='None',
        scale_mode='off', scale_factor=1, glyph_mode='All Points',
        glyph_name=None):
    """Create a glyph in paraview and return the created object."""

    # Show and hide the base, otherewise there is a bug.
    Show(base, view)
    Hide(base, view)

    # Create the glyph.
    glyph = Glyph(Input=base, GlyphType=glyph_type)
    glyph.OrientationArray = ['POINTS', vectors]
    if scale_mode == 'scalar':
        glyph.ScaleArray = ['POINTS', scalars]
    else:
        glyph.ScaleArray = ['POINTS', vectors]
    glyph.ScaleFactor = scale_factor
    glyph.GlyphMode = glyph_mode

    # Rename if name is given.
    if glyph_name is not None:
        RenameSource(glyph_name, glyph)

    # Show the glyph.
    glyph_display = Show(glyph, view)

    return glyph, glyph_display


def color_by(data, field_name, is_baci):
    """Color the data by a certain field name."""
    if is_baci:
        ColorBy(data, ('POINTS', field_name))
    else:
        ColorBy(data, None)


def load_solid(solid_file):
    """Load the solid file in paraview."""

    # Check the extension of the file to see if it is from baci or meshpy.
    extension = os.path.splitext(solid_file)[1].lower()
    if extension == '.pvd':
        # Create the PVD reader.
        solid_reader = PVDReader(FileName=solid_file)
        is_baci = True
    elif extension == '.vtu':
        # Create the VTU reader.
        solid_reader = XMLUnstructuredGridReader(FileName=solid_file)
        is_baci = False
    else:
        raise TypeError('Extension has to be "pvd" or "vtu", '
            + 'got {}!'.format(extension))

    # Create the PVD reader.
    pvd_structure = solid_reader

    # Set the name in paraview to the filename.
    RenameSource(os.path.basename(solid_file), pvd_structure)

    # Apply displacements.
    if is_baci:
        warp_solid = WarpByVector(Input=pvd_structure)
        warp_solid.Vectors = ['POINTS', 'displacement']
        warp_solid.ScaleFactor = 1.0
        solid_display = Show(warp_solid, view)
    else:
        solid_display = Show(pvd_structure, view)

    # Set display style of solid.
    solid_display.Representation = 'Surface With Edges'

    # Set opacity.
    solid_display.Opacity = 0.5

    # Set nonlinear subdivision for quadratic elements.
    solid_display.NonlinearSubdivisionLevel = 4


def load_beam(beam_file):
    """Load the beam file in paraview."""

    # Check the extension of the file to see if it is from baci or meshpy.
    extension = os.path.splitext(beam_file)[1].lower()
    if extension == '.pvd':
        # Create the PVD reader.
        beam_reader = PVDReader(FileName=beam_file)
        is_baci = True
    elif extension == '.vtu':
        # Create the VTU reader.
        beam_reader = XMLUnstructuredGridReader(FileName=beam_file)
        is_baci = False
    else:
        raise TypeError('Extension has to be "pvd" or "vtu", '
            + 'got {}!'.format(extension))

    # Set the name in paraview to the filename.
    RenameSource(os.path.basename(beam_file), beam_reader)

    # Add cell data to point data filter.
    beam_cell_to_point = CellDatatoPointData(Input=beam_reader)

    # Extract the surface from the beam line.
    beam_extract_surface = ExtractSurface(Input=beam_cell_to_point)

    # Represent the beam curve as a tube.
    beam_tube = Tube(Input=beam_extract_surface)
    beam_tube.Scalars = ['POINTS', 'cross_section_radius']
    beam_tube.VaryRadius = 'By Absolute Scalar'
    beam_tube.NumberofSides = 10
    RenameSource('beam', beam_tube)

    # Show the tube.
    beam_tube_display = Show(beam_tube, view)
    color_by(beam_tube_display, 'displacement', is_baci)

    # Create glyphs for triads.
    create_glyph(beam_extract_surface, vectors='base_vector_1', scale_mode='scalar',
        scalars='cross_section_radius', scale_factor=6, glyph_name='base1')
    create_glyph(beam_extract_surface, vectors='base_vector_2', scale_mode='scalar',
        scalars='cross_section_radius', scale_factor=6, glyph_name='base2')
    create_glyph(beam_extract_surface, vectors='base_vector_3', scale_mode='scalar',
        scalars='cross_section_radius', scale_factor=6, glyph_name='base3')

    # Add programmable filter to get nodal points.
    node_filter = ProgrammableFilter(Input=beam_extract_surface)
    node_filter.Script = 'execfile("/home/ivo/dev/paraview-scripts/pf_extract_beam_nodes.py")'
    RenameSource('filter_nodes', node_filter)
    glyph, glyph_display = create_glyph(node_filter, glyph_type='Sphere',
        scale_mode='scalar', scalars='cross_section_radius', scale_factor=3,
        glyph_name='nodes')
    color_by(glyph_display, 'displacement', is_baci)


def load_beam_to_solid(beam_to_solid_file):
    """Load the beam file in paraview."""

    # Factor for glyph representation.
    glyph_factor = 0.01

    # Create the PVD reader.
    pvd_bts = PVDReader(FileName=beam_to_solid_file)

    # Show and hide so that the data for the glyps will be kept (bug in
    # paraview)
    Show(pvd_bts, view)
#        Hide(pvd_bts, view)

    # Set the name in paraview to the filename.
    RenameSource(os.path.basename(beam_to_solid_file), pvd_bts)

    # Add spheres for gauss points.
    create_glyph(pvd_bts, scale_mode='scalar',
        scale_factor=glyph_factor, scalars='gauss_point_value',
        glyph_name='gauss_points', glyph_type='Sphere')

    # Add spheres for segmentation points.
    create_glyph(pvd_bts, scale_mode='scalar',
        scale_factor=glyph_factor, scalars='segmentation_point_value',
        glyph_name='segmentation_points', glyph_type='Sphere')

    # Create fields for the glyph representation.
    search_value = Calculator(Input=pvd_bts)
    search_value.Function = 'ceil(0.5 * search_point_value)'.format(
        glyph_factor)
    search_value.ResultArrayName = 'search_point_glyph_value'
    Show(search_value, view)
    Hide(search_value, view)

    # Add spheres for search points.
    glyph, glyph_display = create_glyph(search_value, scale_mode='scalar',
        scale_factor=glyph_factor, scalars='search_point_glyph_value',
        glyph_name='search_points', glyph_type='Sphere')

    # Set the coloring of the search points.
    ColorBy(glyph_display, ('POINTS', 'search_point_value'))


if __name__ == '__main__':
    """Execute part of script."""

    # Disable reset on show.
    _DisableFirstRenderCameraReset()

    # Render view of application.
    view = GetActiveViewOrCreate('RenderView')

    # Get the filenames. Check if the files were given via environment
    # parameters.
    # Set the environment parameters with:
    # export PV_STRUCTURE="structure.pdv"
    parameters = [
        ['PV_STRUCTURE', 'xxx-structure.pvd'],
        ['PV_BEAM', 'xxx-structure-beams.pvd'],
        ['PV_BEAM_TO_SOLID', 'xxx-beam-to-solid-mesh-tying.pvd']
        ]
    file_names = []
    for key, default in parameters:
        if key in os.environ.keys():
            file_names.append(os.environ[key])
        else:
            file_names.append(default)

    # Check if absolute paths were given or relative paths.
    for i in range(len(file_names)):
        if not os.path.isabs(file_names[i]):
            file_names[i] = os.path.join(os.environ['PWD'], file_names[i])
    solid_file, beam_file, beam_to_solid_file = file_names

    # Load the files to paraview
    if os.path.isfile(solid_file):
        load_solid(solid_file)
    if os.path.isfile(beam_file):
        load_beam(beam_file)
    if os.path.isfile(beam_to_solid_file):
        load_beam_to_solid(beam_to_solid_file)

    # Reset view to loaded geometry.
    view.Update()
    view.ResetCamera()
