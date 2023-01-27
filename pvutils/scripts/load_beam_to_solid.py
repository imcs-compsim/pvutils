# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# PVUTILS: Utility tools for the ParaView python interface
#
# Copyright 2023 Ivo Steinbrecher & Matthias Mayr
#                 Institute for Mathematics and Computer-based Simulation
#                 University of the Bundeswehr Munich
#                 https://www.unibw.de/imcs-en
#
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# -----------------------------------------------------------------------------
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
    pvutils.display(
        solid, line_color=[1.0, 1.0, 1.0], representation="Surface With Edges"
    )
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

    pvutils.update_scene()

    # Return the created objects.
    return solid, beam


def load_beam_to_solid_in_dir(path):
    """
    Load a beam and solid file in a given path.
    The beam file is expected to contain "-structure-beams.pvd", while the
    solid file should contain "-structure.pvd". There has to be exactly ONE
    beam and solid file in the working directory.
    """

    # Get all pvd files in the working directory.
    pvd_files = [file for file in os.listdir(path) if file.endswith("pvd")]

    # Check if there is ONE structure and beam file.
    structure_files = [file for file in pvd_files if "-structure.pvd" in file]
    beam_files = [file for file in pvd_files if "-structure-beams.pvd" in file]
    if not len(structure_files) == 1 and not len(beam_files) == 1:
        raise NameError(
            "Could not find a single solid and beam file. "
            "Found {} solid and {} beam files.".format(
                len(structure_files), len(beam_files)
            )
        )

    # Call the script.
    return load_beam_to_solid(
        os.path.join(path, structure_files[0]), os.path.join(path, beam_files[0])
    )


if __name__ == "__main__":
    """Execute part of script."""

    load_beam_to_solid_in_dir(os.environ["PWD"])
