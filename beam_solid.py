# -*- coding: utf-8 -*-
"""
Load a solid and a beam into paraview and apply some filters.
"""

# Import the simple module from the paraview.
from paraview.simple import GetActiveViewOrCreate, PVDReader, RenameSource, \
    WarpByVector, Show, _DisableFirstRenderCameraReset, ColorBy, Glyph, \
    ExtractSurface, Tube, Calculator, CellDatatoPointData, Hide


# Python imports.
import os


def create_glyph(base, glyph_type='Arrow', scalars='None', vectors='None',
        scale_mode='off', scale_factor=1, glyph_mode='All Points',
        glyph_name=None):
    """Create a glyph in paraview and return the created object."""

    # Create the glyph.
    glyph = Glyph(Input=base, GlyphType=glyph_type)
    glyph.Scalars = ['POINTS', scalars]
    glyph.Vectors = ['POINTS', vectors]
    glyph.ScaleMode = scale_mode
    glyph.ScaleFactor = scale_factor
    glyph.GlyphMode = glyph_mode

    # Rename if name is given.
    if glyph_name is not None:
        RenameSource(glyph_name, glyph)

    # Show the glyph.
    glyph_display = Show(glyph, view)

    return glyph, glyph_display


def load_solid(solid_file):
    """Load the solid file in paraview."""

    if solid_file is not None:

        # Create the PVD reader.
        pvd_structure = PVDReader(FileName=solid_file)

        # Set the name in paraview to the filename.
        RenameSource(os.path.basename(solid_file), pvd_structure)

        # Apply displacements.
        warp_solid = WarpByVector(Input=pvd_structure)
        warp_solid.Vectors = ['POINTS', 'displacement']
        warp_solid.ScaleFactor = 1.0

        # Set display style of solid.
        solid_display = Show(warp_solid, view)
        solid_display.Representation = 'Surface With Edges'

        # Set opacity.
        solid_display.Opacity = 0.5


def load_beam(beam_file):
    """Load the beam file in paraview."""

    if beam_file is not None:

        # Create the PVD reader.
        pvd_beam = PVDReader(FileName=beam_file)

        # Set the name in paraview to the filename.
        RenameSource(os.path.basename(beam_file), pvd_beam)

        # Add cell data to point data filter.
        beam_cell_to_point = CellDatatoPointData(Input=pvd_beam)

        # Extract the surface from the beam line.
        beam_extract_surface = ExtractSurface(Input=beam_cell_to_point)

        # Create fields for the glyph representation.
        calculator = Calculator(Input=beam_extract_surface)
        calculator.Function = '2.0 * floor(node_value) * cross_section_radius'
        calculator.ResultArrayName = 'nodal_radius'

        # Represent the beam curve as a tube.
        beam_tube = Tube(Input=calculator)
        beam_tube.Scalars = ['POINTS', 'cross_section_radius']
        beam_tube.VaryRadius = 'By Absolute Scalar'
        beam_tube.NumberofSides = 10
        RenameSource('beam', beam_tube)

        # Show the tube.
        beam_tube_display = Show(beam_tube, view)
        ColorBy(beam_tube_display, ('POINTS', 'displacement'))

        # Create glyphs for triads.
        create_glyph(calculator, vectors='base_vector_1', scale_mode='scalar',
            scalars='cross_section_radius', scale_factor=6, glyph_name='base1')
        create_glyph(calculator, vectors='base_vector_2', scale_mode='scalar',
            scalars='cross_section_radius', scale_factor=6, glyph_name='base2')
        create_glyph(calculator, vectors='base_vector_3', scale_mode='scalar',
            scalars='cross_section_radius', scale_factor=6, glyph_name='base3')

        # Create glyphs for nodal points.
        glyph, glyph_display = create_glyph(calculator, glyph_type='Sphere',
            scale_mode='scalar', scalars='nodal_radius', scale_factor=2,
            glyph_name='nodes')
        ColorBy(glyph_display, ('POINTS', 'node_value'))


def load_beam_to_solid(beam_to_solid_file):
    """Load the beam file in paraview."""

    if beam_file is not None:

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
        search_value.Function = 'ceil(0.5*search_point_value)'.format(
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

    # Load data in paraview.
    solid_file = '/home/ivo/temp/xxx-structure.pvd'
    beam_file = '/home/ivo/temp/xxx-structure-beams.pvd'
    beam_to_solid_file = '/home/ivo/temp/xxx-beam-to-solid-mesh-tying.pvd'
    load_solid(solid_file)
    load_beam(beam_file)
    load_beam_to_solid(beam_to_solid_file)

    # Reset view to loaded geometry.
    view.Update()
    view.ResetCamera()
