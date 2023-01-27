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
This script generates a point with vector data. Can be used to visualize
coordinate axes.

Use with programmable source filter in paraview and execute with:
>> exec(open("<path to file>").read())
"""

from paraview import vtk


def get_keyword_argument(name, kwargs_id):
    """
    Return the wanted keyword argument.
    """
    import paraview

    if hasattr(paraview, "programmable_source_kwargs"):
        if name in paraview.programmable_source_kwargs[kwargs_id].keys():
            return paraview.programmable_source_kwargs[kwargs_id][name]
    raise ValueError('Keyword argument "{}" not found'.format(name))


# Get keyword arguments.
origin = get_keyword_argument("origin", kwargs_id)
basis = get_keyword_argument("basis", kwargs_id)

# Add the origin.
output = self.GetOutput()
newPts = vtk.vtkPoints()
newPts.InsertPoint(0, origin[0], origin[1], origin[2])
output.SetPoints(newPts)

# Add the basis vectors.
for i, base in enumerate(basis):
    base_data = vtk.vtkFloatArray()
    base_data.SetNumberOfComponents(3)
    base_data.SetName("base_{}".format(i + 1))
    base_data.InsertNextTuple3(base[0], base[1], base[2])
    output.GetPointData().AddArray(base_data)
