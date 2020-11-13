# -*- coding: utf-8 -*-
"""
Load a beam into ParaView.
"""

# Python imports.
import os

# Local imports.
import utility_functions
import filter_wrapper

# ParaView imports.
import paraview.simple as pa


class BeamDisplay(object):
    """
    Class to handle the display of a beam in ParaView.
    """

    def __init__(self, beam_file, segments=8, merge_poly_lines=False,
            triads=True, factor_triads=3.0, show_triads=False,
            nodes=True, factor_nodes=3.0, show_nodes=False
            ):
        """
        Load a beam into ParaView.

        Args
        ----
        beam_file: str
            Path to the beam file that should be displayed.
        segments: int
            Number of segments to be used for the beam tube and the node
            spheres.
        merge_poly_lines: bool
            If the clean to grid and merge poly lines filter should be applied.
        triads: bool
            If the basis vectors should be displayed.
        factor_triads: float
            Factor between triad length and beam radius.
        show_triads: bool
            If triads are created, should they be visible or not.
        nodes: bool
            If the nodes of the beam elements (start and end node) should be
            displayed.
        factor_nodes: float
            Factor between node radius and beam radius.
        show_nodes: bool
            If nodes are created, should they be visible or not.
        """

        # Display items of this class.
        self.beam = utility_functions.load_file(beam_file)
        self.beam_cell_to_point = pa.CellDatatoPointData(Input=self.beam)
        if merge_poly_lines:
            self.beam_clean_to_grid = pa.CleantoGrid(
                Input=self.beam_cell_to_point)
            self.beam_merge_poly_line = utility_functions.programmable_filter(
                self.beam_clean_to_grid, 'merge_polylines')
            next_input = self.beam_merge_poly_line
        else:
            self.beam_clean_to_grid = None
            self.beam_merge_poly_line = None
            next_input = self.beam_cell_to_point
        self.beam_extract_surface = pa.ExtractSurface(Input=next_input)
        self.beam_tube = filter_wrapper.tube(self.beam_extract_surface)
        self.endpoints = None
        self.nodes = None
        self.base_vectors = []

        # Set the options for the tube filter.
        pa.UpdatePipeline()
        field_names = utility_functions.get_field_names(
            self.beam_cell_to_point)
        if ('cross_section_radius', 1) not in field_names['POINTS']:
            raise ValueError('Could not find cross_section_radius for tube!')
        self.beam_tube.Scalars = ['POINTS', 'cross_section_radius']
        self.beam_tube.VaryRadius = 'By Absolute Scalar'
        self.beam_tube.Radius = 1.0
        self.beam_tube.RadiusFactor = 1.0
        self.beam_tube.NumberofSides = segments
        utility_functions.display(self.beam_tube)

        # Display the nodes as spheres
        if nodes:
            self.endpoints = utility_functions.programmable_filter(
                self.beam_cell_to_point, 'get_polyline_endpoints')

            self.nodes = filter_wrapper.glyph(self.endpoints)
            self.nodes.GlyphType = 'Sphere'
            self.nodes.ScaleArray = ['POINTS', 'cross_section_radius']
            self.nodes.ScaleFactor = factor_nodes
            self.nodes.GlyphMode = 'All Points'
            self.nodes.GlyphType.ThetaResolution = segments
            self.nodes.GlyphType.PhiResolution = segments
            pa.RenameSource('nodes', self.nodes)
            if show_nodes:
                utility_functions.display(self.nodes)

        # Display the basis vector of the triads.
        if triads:
            for i in range(3):
                name = 'base_vector_{}'.format(i + 1)
                if (name, 3) in field_names['POINTS']:
                    base_vector = filter_wrapper.glyph(
                        self.beam_cell_to_point)
                    base_vector.GlyphType = 'Arrow'
                    base_vector.OrientationArray = ['POINTS', name]
                    base_vector.ScaleArray = ['POINTS',
                        'cross_section_radius']
                    base_vector.ScaleFactor = factor_triads
                    base_vector.GlyphMode = 'All Points'
                    pa.RenameSource(name, base_vector)
                    if show_triads:
                        utility_functions.display(base_vector)
                    self.base_vectors.append(base_vector)

        # Set the beam tube as active element.
        pa.SetActiveSource(self.beam_tube)
