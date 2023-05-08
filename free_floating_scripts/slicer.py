# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# PVUTILS: Utility tools for the ParaView python interface
#
# Copyright 2023 Ivo Steinbrecher & Matthias Mayr
#                Institute for Mathematics and Computer-based Simulation
#                University of the Bundeswehr Munich
#                https://www.unibw.de/imcs-en
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
Slices a 3d geometry with a series of planes.

This script is intended (or at least tested) to be called from the python shell
of the paraview GUI. The currently active item, e.g., highlighted in the
Pipeline Browser of the GUI, will be sliced.
"""

# Python imports.
import os
import numpy as np
from pathlib import Path

# ParaView imports.
import paraview.simple as pa

# Define the slicing with the following user inputs.
# Slicing planes are user-defined by the `slice_normal` and `number_of_slices`
# positions along this direction between `start_point` and `end_point`. Bounds
# may be perturbed to ensure that slices are within the volume.
slice_normal = [0, 0, 1]
number_of_slices = 5
start_point = 0.0
end_point = 1.0
perturb_bounds = True

# Results to be displayed on the slieced geometry.
coloring_specified = False
coloring = ("POINTS", "displacement", "Magnitude")

# Where to save the resulting slices?
output_path = str(Path.home())
output_name = "slice"


z_s = np.linspace(start_point, end_point, number_of_slices)

if perturb_bounds:
    # Perturb start and end point to ensure that slice will be inside volume.
    delta = (z_s[-1] - z_s[0]) / number_of_slices
    z_s[0] = z_s[0] + 1.0e-3 * delta
    z_s[-1] = z_s[-1] - 1.0e-3 * delta

# At this point we are looking at the volume data we want to slice.
source = GetActiveSource()

# Now, we hide it as it will block our view otherwise.
renderView1 = GetActiveViewOrCreate("RenderView")
Hide(source, renderView1)

# Apply slice.
slice = pa.Slice(source)
slice.SliceType = "Plane"
slice.SliceType.Normal = slice_normal

# Customize rendering view.
slice1 = GetActiveSource()
renderView1 = GetActiveViewOrCreate("RenderView")
slice1Display = Show(slice1, renderView1, "GeometryRepresentation")
Hide3DWidgets(proxy=slice1.SliceType)
if coloring_specified:
    ColorBy(slice1Display, coloring)

for i, z in enumerate(z_s):

    print("Slicing at z = " + str(z) + ".")
    slice.SliceType.Origin = [x * z for x in slice_normal]

    renderView1.Update()

    SaveScreenshot(
        os.path.join(output_path, output_name + str(i) + ".png"), renderView1
    )
