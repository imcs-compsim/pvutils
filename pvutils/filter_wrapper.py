# -*- coding: utf-8 -*-
"""
Function  wrappers for ParaView functions.
"""


# Import ParaView module.
import paraview.simple as pa

# Import pvutils functions.
from .utils import check_data


def warp(data, field='displacement', scale_factor=1.0):
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

    check_data(data, field, dimension=3)
    warp = pa.WarpByVector(Input=data)
    warp.Vectors = ['POINTS', field]
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


def threshold(data, field='displacement', data_type='POINTS',
        threshold_range=None):
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

    if threshold_range is None:
        threshold_range = [-1.0e+12, 1.0e+12]

    threshold = pa.Threshold(Input=data)
    threshold.Scalars = [data_type, field]
    threshold.ThresholdRange = threshold_range
    return threshold


def tube(data, slices=8):
    """
    Apply the tube filter to a beam.
    """

    check_data(data, 'cross_section_radius', dimension=1)
    tube = pa.Tube(Input=data)
    tube.Scalars = [None, 'cross_section_radius']
    tube.VaryRadius = 'By Absolute Scalar'
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
