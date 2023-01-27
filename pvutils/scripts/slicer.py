#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Slices a 3d geometry with a series of planes.
"""

# Python imports.
import numpy as np
from pathlib import Path

# ParaView imports.
import paraview.simple as pa

# Define slicing
start_point = 0.0
end_point = 1.0
number_of_slices = 5
slice_normal = [0, 0, 1]
perturb_bounds = True
coloring_specified = False
coloring = ("POINTS", "displacement", "Magnitude")
output_path = str(Path.home()) + "/"
output_name = "slice"

z_s = np.linspace(z_0, z_N, number_of_slices)

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
s = pa.Slice(source)
s.SliceType = "Plane"
s.SliceType.Normal = slice_normal

# Customize rendering view.
slice1 = GetActiveSource()
renderView1 = GetActiveViewOrCreate("RenderView")
slice1Display = Show(slice1, renderView1, "GeometryRepresentation")
# Turn of plane rendering.
Hide3DWidgets(proxy=slice1.SliceType)
# set scalar coloring
if coloring_specified:
    ColorBy(slice1Display, coloring)

for i, z in enumerate(z_s):

    print("Slicing at z = " + str(z) + ".")
    s.SliceType.Origin = [x * z for x in slice_normal]

    renderView1.Update()

    SaveScreenshot(output_path + output_name + str(i) + ".png", renderView1)
