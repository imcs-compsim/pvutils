# -*- coding: utf-8 -*-
"""
Load a beam into ParaView.
"""


# Python imports.
import os

# ParaView imports.
import pvutils
import paraview.simple as pa


class BeamDisplay(object):
    """
    Class to handle the display of a beam in ParaView.
    """

    def __init__(self, beam_file, triads=True, nodes=True, segments=8,
            factor_nodes=3.0, factor_triads=3.0):
        """
        Load a beam into ParaView.

        Args
        ----
        beam_file: str
            Path to the beam file that should be displayed.
        triads: bool
            If the basis vectors should be displayed.
        nodes: bool
            If the nodes of the beam elements (start and end node) should be
            displayed.
        segments: int
            Number of segments to be used for the beam tube and the node
            spheres.
        """

        # Display items of this class.
        self.beam = pvutils.load_file(beam_file)
        self.beam_cell_to_point = pa.CellDatatoPointData(Input=self.beam)
        self.beam_extract_surface = pa.ExtractSurface(
            Input=self.beam_cell_to_point)
        self.beam_tube = pvutils.tube(self.beam_extract_surface)
        self.endpoints = None
        self.nodes = None
        self.base_vectors = []

        # Set the options for the tube filter.
        pa.UpdatePipeline()
        field_names = pvutils.get_field_names(self.beam_cell_to_point)
        if ('cross_section_radius', 1) not in field_names['POINTS']:
            raise ValueError('Could not find cross_section_radius for tube!')
        self.beam_tube.Scalars = ['POINTS', 'cross_section_radius']
        self.beam_tube.VaryRadius = 'By Absolute Scalar'
        self.beam_tube.Radius = 1.0
        self.beam_tube.RadiusFactor = 1.0
        self.beam_tube.NumberofSides = segments
        pvutils.display(self.beam_tube)

        # Display the nodes as spheres
        if nodes:
            self.endpoints = pvutils.programmable_filter(
                self.beam_cell_to_point, 'get_polyline_endpoints')

            self.nodes = pvutils.glyph(self.endpoints)
            self.nodes.GlyphType = 'Sphere'
            self.nodes.ScaleArray = ['POINTS', 'cross_section_radius']
            self.nodes.ScaleFactor = factor_nodes
            self.nodes.GlyphMode = 'All Points'
            self.nodes.GlyphType.ThetaResolution = segments
            self.nodes.GlyphType.PhiResolution = segments
            pa.RenameSource('nodes', self.nodes)
            pvutils.display(self.nodes)

        # Display the basis vector of the triads.
        if triads:
            for i in range(3):
                name = 'base_vector_{}'.format(i + 1)
                if (name, 3) in field_names['POINTS']:
                    base_vector = pvutils.glyph(self.beam_cell_to_point)
                    base_vector.GlyphType = 'Arrow'
                    base_vector.OrientationArray = ['POINTS', name]
                    base_vector.ScaleArray = ['POINTS',
                        'cross_section_radius']
                    base_vector.ScaleFactor = factor_triads
                    base_vector.GlyphMode = 'All Points'
                    pa.RenameSource(name, base_vector)
                    pvutils.display(base_vector)
                    self.base_vectors.append(base_vector)
