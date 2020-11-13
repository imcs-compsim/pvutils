# -*- coding: utf-8 -*-
"""
Create the reference files for the merge polylines filter.
"""

# Python modules.
import os

# Meshpy modules.
from meshpy import Beam3rHerm2Line3, MaterialReissner, Mesh
from meshpy.mesh_creation_functions import create_beam_mesh_line


def create_reference_file_merge_polyline():
    """
    Create the half circle solid and beam.
    """

    base_dir = os.path.dirname(os.path.realpath(__file__))
    ref_dir = os.path.join(base_dir, 'reference-files')

    mesh = Mesh()
    mat = MaterialReissner(radius=0.1)

    n_el = 2

    # Create a single beam.
    z = -2
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 0, z],
        [1, 0, z],
        n_el=1)

    # Create two beams.
    z = -1
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 0, z],
        [1, 0, z],
        n_el=n_el)

    # Create a closed path.
    z = 0
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 0, z],
        [1, 0, z],
        n_el=n_el)
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [1, 0, z],
        [1, 1, z],
        n_el=n_el)
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [1, 1, z],
        [0, 1, z],
        n_el=n_el)
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 1, z],
        [0, 0, z],
        n_el=n_el)

    # Create check all possible connection points.
    z = 1
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 0, z],
        [0, 1, z],
        n_el=n_el)
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 1, z],
        [0.5, 2, z],
        n_el=n_el)

    z = 2
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 0, z],
        [0, 1, z],
        n_el=n_el)
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0.5, 2, z],
        [0, 1, z],
        n_el=n_el)

    z = 3
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 1, z],
        [0, 0, z],
        n_el=n_el)
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 1, z],
        [0.5, 2, z],
        n_el=n_el)

    z = 4
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 1, z],
        [0, 0, z],
        n_el=n_el)
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0.5, 2, z],
        [0, 1, z],
        n_el=n_el)

    # Create a simple bifurcation.
    z = 5
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 0, z],
        [0, 1, z],
        n_el=n_el)
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 1, z],
        [0.5, 2, z],
        n_el=n_el)
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 1, z],
        [-0.5, 2, z],
        n_el=n_el)

    # Create a bifurcation with a closed circle.
    z = 6
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 0, z],
        [0, 1, z],
        n_el=n_el)
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 1, z],
        [0.5, 2, z],
        n_el=n_el)
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 1, z],
        [-0.5, 2, z],
        n_el=n_el)
    create_beam_mesh_line(mesh, Beam3rHerm2Line3, mat,
        [0, 0, z],
        [-0.5, 2, z],
        n_el=n_el)

    mesh.write_vtk('merge_polylines_raw', ref_dir, ascii=True)


if __name__ == '__main__':
    """Execution part of script."""

    create_reference_file_merge_polyline()
