# -*- coding: utf-8 -*-
"""
Load a beam into ParaView.
"""

# Python imports.
import os

# ParaView imports.
import pvutils
import paraview.simple as pa


def load_beam_to_solid(solid_path, beam_path):
    """
    Load a beam and a solid together, and display them.
    """

    # Load the solid.
    solid = pvutils.load_file(solid_path)
    solid = pvutils.warp(solid)
    pvutils.display(solid, line_color=[1.0, 1.0, 1.0],
        representation='Surface With Edges')
    pvutils.contour(solid)

    # Load the beams.
    beam = pvutils.BeamDisplay(beam_path)
    pvutils.contour(beam.beam_tube)

    # Set to first time step.
    pvutils.set_timestep(pvutils.get_available_timesteps()[0])

    # Reset the view to fit the loaded data.
    view = pvutils.get_view()
    view.Update()
    view.ResetCamera()


if __name__ == '__main__':
    """Execute part of script."""

    # Get the working directory.
    work_dir = os.environ['PWD']

    # Get all pvd files in the working directory.
    pvd_files = [file for file in os.listdir(work_dir) if file.endswith('pvd')]

    # Check if there is ONE structure and beam file.
    structure_files = [file for file in pvd_files if '-structure.pvd' in file]
    beam_files = [file for file in pvd_files if '-structure-beams.pvd' in file]
    if not len(structure_files) == 1 and not len(beam_files) == 1:
        raise NameError('Could not find a single solid and beam file. '
            'Found {} solid and {} beam files.'.format(
                len(structure_files), len(beam_files)))

    # Call the script.
    load_beam_to_solid(
        os.path.join(work_dir, structure_files[0]),
        os.path.join(work_dir, beam_files[0])
        )
