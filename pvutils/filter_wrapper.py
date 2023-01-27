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
Function wrappers for ParaView functions.
"""


# Import ParaView module.
import paraview.simple as pa

# Import pvutils functions.
from . import utility_functions


def glyph(data):
    """
    Apply a glyph filter to data.

    Args
    ----
    data: ParaView data object
        ParaView item the glyph filter will be applied to.
    """

    # This show and hide is needed, because otherwise the scalar data settings
    # in the glyph will not be applied correctly (ParaView bug).
    pa.Show(data)
    pa.Hide(data)

    # Create the filter and set the desired options.
    glyph = pa.Glyph(Input=data)
    return glyph


def warp(data, field="displacement", scale_factor=1.0):
    """
    Warp the data by a vector field.

    Args
    ----
    data: ParaView data object
        ParaView item the warp filter will be applied to.
    field: String
        Name of vector field to be used as "displacement" of the warp.
    scale_factor: Float
        Scaling factor.

    Return
    ------
    warp: ParaView data object
        The resulting ParaView item after applying the filter.
    """

    utility_functions.check_data(data, field, dimension=3)
    warp = pa.WarpByVector(Input=data)
    warp.Vectors = ["POINTS", field]
    warp.ScaleFactor = scale_factor
    return warp


def transform(data, translate=None, rotate=None, scale=None):
    """
    Apply a 'Transform' filter to a ParaView item.

    Args
    ----
    data: ParaView data object
        ParaView item the transpose filter will be applied to.
    transform: 3D-array
        Translation increment vector.
    rotate: 3D array
        Rotation increment vector.
    scale:
        Scaling in each Cartesian direction.

    Return
    ------
    transform: ParaView data object
        The resulting ParaView item after applying the filter.
    """

    if translate is None:
        translate = [0.0, 0.0, 0.0]
    if rotate is None:
        rotate = [0.0, 0.0, 0.0]
    if scale is None:
        scale = [1.0, 1.0, 1.0]

    transform = pa.Transform(Input=data)
    transform.Transform.Translate = translate
    transform.Transform.Rotate = rotate
    transform.Transform.Scale = scale
    return transform


def threshold(data, field="displacement", data_type="POINTS", threshold_range=None):
    """
    Apply a 'Threshold' filter to a ParaView item.

    Args
    ----
    data: ParaView data object
        ParaView item the threshold filter will be applied to.
    field: String
        Name of field whose values will be subject to thresholding.
    data_type: String
        Type/topology of data, encoded using ParaView-style namings.
    threshold_range: Array with 2 entries
        Array with min and max values of valid range of values.

    Return
    ------
    threshold: ParaView data object
        The resulting ParaView item after applying the filter.
    """

    # This show and hide is needed, because otherwise the scalar data settings
    # in the threshold will not be applied correctly (ParaView bug).
    # Todo: check if this fixes it
    # pa.UpdatePipeline()
    pa.Show(data)
    pa.Hide(data)

    threshold = pa.Threshold(Input=data)
    threshold.Scalars = [data_type, field]

    if threshold_range is None:
        threshold.LowerThreshold = -1.0e12
        threshold.UpperThreshold = 1.0e12
    else:
        threshold.LowerThreshold = threshold_range[0]
        threshold.UpperThreshold = threshold_range[1]

    return threshold


def clip(data, clip_type="Plane", origin=None, normal=None, invert=True):
    """
    Apply a 'clip' filter to a ParaView item.

    Args
    ----
    data: ParaView data object
        ParaView item the threshold filter will be applied to.
    clip_type: str
        Type of clipping operation. Not all types implemented, yet.
    origin: list
        Coordinates of origin of the clipping geometry
    normal: list
        Normal vector to define the clipping geometry
    invert: bool
        Flag to indicate "invert" of the clipping output
    """

    if clip_type != "Plane":
        raise Exception(
            "Clipt types other than 'plane' are not available, yet. Please implement them."
        )
    if origin is None:
        origin = [0.0, 0.0, 0.0]
    if normal is None:
        normal = [1.0, 0.0, 0.0]

    clip = pa.Clip(Input=data)
    clip.ClipType = "Plane"
    clip.ClipType.Origin = origin
    clip.ClipType.Normal = normal
    clip.Invert = int(invert)
    return clip


def tube(data, slices=8):
    """
    Apply the tube filter to a beam.
    """

    utility_functions.check_data(data, "cross_section_radius", dimension=1)
    tube = pa.Tube(Input=data)
    tube.Scalars = [None, "cross_section_radius"]
    tube.VaryRadius = "By Absolute Scalar"
    tube.NumberofSides = slices
    return tube


def temporal_interpolator(data):
    """
    Apply the temporal data interpolator to an item. It is necessary to show
    and hide the resulting item, as otherwise ParaView will raise an error as
    soon as field data is requested.
    """

    temp_interpolator = pa.TemporalInterpolator(Input=data)
    pa.Show(temp_interpolator)
    pa.Hide(temp_interpolator)
    return temp_interpolator
