# -*- coding: utf-8 -*-
"""
Load a beam into ParaView.
"""

# Python imports.
import os
import numpy as np

# Local imports.
from . import utility_functions
from . import filter_wrapper

# ParaView imports.
import paraview.simple as pa


class BeamDisplay(object):
    """
    Class to handle the display of a beam in ParaView.
    """

    def __init__(self, beam_file, segments=8, merge_poly_lines=False,
            merge_polylines_max_angle=None,
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
                self.beam_clean_to_grid, pvutils_filter='merge_polylines',
                max_angle=merge_polylines_max_angle)
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
                self.beam_cell_to_point,
                pvutils_filter='get_polyline_endpoints')

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


class Polynomial3DCurve(object):
    """
    This class helps get the polynomial representation of a third- order curve
    in 3D space.
    """

    def __init__(self, n_segments=5, order=3):
        """
        Setup the class and the interpolation variables.

        Args
        ----
        n_segments: int
            Number of segments used to visualize the curve. It is assumed that
            the segmentation points are equidistant in the parameter space.
        order: int
            Polynomial order of the curve.
        """

        self.order = order

        # Parameter positions of the points that will be given.
        self.xi = np.array(
            [-1, -1 + 2.0 / n_segments, 1 - 2.0 / n_segments, 1])

        # Interpolation matrix.
        A = np.zeros([self.order + 1, self.order + 1])
        for i_row in range(self.order + 1):
            for i_col in range(self.order + 1):
                A[i_row, i_col] = self.xi[i_row] ** i_col

        # Inverted interpolation matrix.
        self.A_inv = np.linalg.inv(A)

        # Coefficient matrix.
        self.coefficiens = np.zeros([self.order + 1, 3])

    def eval_r(self, cell_positions, xi):
        """
        Evaluate the position on the curve at xi.
        """

        for i_dir in range(3):
            self.coefficiens[:, i_dir] = np.dot(
                self.A_inv, cell_positions[:, i_dir])

        temp = np.zeros(3)
        for i_coeff in range(self.order + 1):
            temp += self.coefficiens[i_coeff, :] * xi ** i_coeff
        return temp

    def eval_rp(self, cell_positions, xi):
        """
        Calculate the derivative of the position at xi.
        """

        for i_dir in range(3):
            self.coefficiens[:, i_dir] = np.dot(
                self.A_inv, cell_positions[:, i_dir])

        temp = np.zeros(3)
        for i_coeff in range(1, self.order + 1):
            temp += (self.coefficiens[i_coeff, :] * xi ** (i_coeff - 1) *
                i_coeff)
        return temp
