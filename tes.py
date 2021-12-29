# -*- coding: utf-8 -*-
"""
This script is used to test the functionality of pvutils.
"""

# Python imports.
import unittest
import os
import shutil
import hashlib
import numpy as np
import sys

# ParaView imports.
import pvutils
import paraview.simple as pa
import vtk
from vtk.util import numpy_support as vtk_numpy

def _empty_temp_testing_directory():
    """Delete all files in the testing directory, if it exists."""
    if os.path.isdir(testing_temp):
        for the_file in os.listdir(testing_temp):
            file_path = os.path.join(testing_temp, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(e)

def test_full_script_beam_to_solid_volume_meshtying():
    """
    Load a patch test of a solid structure with beams, apply some filters
    to it and check the created screenshot. This will test if the complete
    workflow works.
    """

    # Load the solid and beam mesh.
    solid = pvutils.load_file(os.path.join(testing_reference,
        'solid_cube_case', 'solid.case'))
    beam = pvutils.load_file(os.path.join(testing_reference,
        'beam_helix_pvd', 'beam.pvd'))
    scene = pvutils.update_scene()
    scene.GoToLast()


    # Apply the tube filter to the beam and display the curvature.
    scale_factor = 50.0
    beam = pvutils.warp(beam, scale_factor=scale_factor)
    beam = pa.CellDatatoPointData(Input=beam)
    beam = pa.ExtractSurface(Input=beam)
    beam = pvutils.tube(beam, slices=20)
    pvutils.display(beam)
    pvutils.contour(beam, 'curvature_2_GPs', vector_type='Y')

    # The solid has some unwanted lines in it. Use the add_id filter to get
    # the cell ids and apply a threshold to the ids.
    solid_filter = pvutils.programmable_filter(solid, 'add_id')
    solid_filter = pa.Threshold(Input=solid_filter)
    pa.UpdatePipeline()
    solid_filter.Scalars = ['CELLS', 'cell_id']
    solid_filter.ThresholdRange = [0, 111]
    pvutils.display(solid_filter, solid_color=[0.0, 0.0, 0.0],
        line_width=2, representation='Outline')

    # Cut half of the solid, so the beam inside it will be visible.
    solid_cut = pa.Clip(Input=solid_filter)
    solid_cut.ClipType = 'Plane'
    solid_cut.ClipType.Origin = [0.00001, -0.00001, 0.0]
    solid_cut.ClipType.Normal = [1.0, -1.0, 0.0]
    solid_cut = pvutils.warp(solid_cut, scale_factor=scale_factor)
    pvutils.contour(solid_cut, 'nodal_2PK_stresses_xyz', vector_type='ZZ')
    pvutils.display(solid_cut, representation='Surface With Edges',
        line_color=[1, 1, 1])

    # Set the view options.
    view = pa.GetActiveViewOrCreate('RenderView')
    view.CameraPosition = [2.76685, -3.24653, 2.75741]
    view.CameraFocalPoint = [0.304706, 0.19705, 0.957474]
    view.CameraViewUp = [-0.230935, 0.320594, 0.918634]
    view.CameraViewAngle = 30
    view.CameraParallelScale = 1.23587
    view.OrientationAxesVisibility = 0
    view.CameraParallelProjection = 0

    # Set the size of the desired output picture.
    dpi = 300
    size = [6, 6]  # This is in cm.
    size_pixel = pvutils.get_size_pixel(size, dpi)
    font_size = [10, 8]
    view.ViewSize = size_pixel

    # Set and place the color map for the solid.
    display = pa.GetDisplayProperties(solid_cut, view=view)
    display.SetScalarBarVisibility(view, True)
    color_function_sigma = pa.GetColorTransferFunction(
        'nodal_2PK_stresses_xyz')
    color_bar_stress = pa.GetScalarBar(color_function_sigma, view)
    color_bar_stress.Title = '$S_{ZZ}$'
    color_bar_stress.ComponentTitle = ''
    color_bar_stress.WindowLocation = 'AnyLocation'
    color_bar_stress.Position = [0.69, 0.6]
    color_bar_stress.ScalarBarLength = 0.2
    pvutils.set_colorbar_font(color_bar_stress, font_size, dpi, font='TeX')

    # Set and place the color map for the beam.
    display = pa.GetDisplayProperties(beam, view=view)
    display.SetScalarBarVisibility(view, True)
    color_function_kappa = pa.GetColorTransferFunction('curvature_2_GPs')
    color_bar_kappa = pa.GetScalarBar(color_function_kappa, view)
    color_bar_kappa.Title = '$\kappa$'
    color_bar_kappa.ComponentTitle = ''
    color_bar_kappa.WindowLocation = 'AnyLocation'
    color_bar_kappa.Position = [0.69, 0.2]
    color_bar_kappa.ScalarBarLength = 0.2
    pvutils.set_colorbar_font(color_bar_kappa, font_size, dpi, font='TeX')



    # For the TikZ images, set one color bar to be horizontal, so both
    # cases are tested. Change the placement of the color bar, so
    # ParaView updates the position (this is not done by just changing
    # the orientation.
    color_bar_stress.Orientation = 'Horizontal'
    color_bar_stress.Position = [0.7, 0.6]
    color_bar_stress.Position = [0.69, 0.6]

    # Create the TikZ wrapped image.
    base_path = os.path.join(testing_temp, 'test')
    pvutils.export_to_tikz(
        os.path.join(testing_temp, 'test'),
        dpi=dpi,
        color_transfer_functions=[
            color_function_sigma,
            color_function_kappa
            ],
        number_format='$\\pgfmathprintnumber{\\tick}$')

    # Compare the with the reference image.

    # Compare the created TikZ code.
    with open(base_path + '.tex', 'r') as tikz_file:
        tikz_code = tikz_file.read()


    color_bar_stress_proxy = color_bar_stress.SMProxy
    repr_proxy = color_bar_stress_proxy.GetRepresentationProxy()
    for i in dir(repr_proxy):
        print(i)


    tmp = repr_proxy.GetProperty('ScalarBarActor')

    print('')

    for i in dir(tmp):
        print(i)

    actor = tmp.GetProxy(0).GetClientSideObject()


    print('')

    for i in dir(actor):
        print(i)
    print(actor)


    print(actor.GetDrawTickLabels())
    print(actor.GetDrawTickMarks())


    #print(actor.GetCustomLabels())



if __name__ == '__main__':
    # Execution part of script.

    pa._DisableFirstRenderCameraReset()

    # Define the testing paths.
    testing_path = os.path.abspath(os.getcwd()) + '/tests'
    testing_reference = os.path.join(testing_path, 'reference-files')
    testing_temp = os.path.join(testing_path, 'testing-tmp')

    # Check and clean the temp directory.
    if not os.path.isdir(testing_temp):
        os.makedirs(testing_temp)
    _empty_temp_testing_directory()

    # Perform the tests.
    test_full_script_beam_to_solid_volume_meshtying()
