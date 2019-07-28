# -*- coding: utf-8 -*-
"""
This script is used to test the functionality of pvutils.
"""

# Python imports.
import unittest
import os
import shutil
import hashlib

# ParaView imports.
import pvutils
import paraview.simple as pa
pa._DisableFirstRenderCameraReset()


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


def _get_file_hash(path):
    """Return the md5 hash of a file."""
    hasher = hashlib.md5()
    with open(path, 'rb') as f:
        buf = f.read()
        hasher.update(buf)
    return hasher.hexdigest()


def _is_gitlab():
    """Check if the environment variable for gitlab testing is set."""
    if os.getenv('GITLAB_TESTING', '0') == '1':
        return True
    else:
        return False


class TestPvutils(unittest.TestCase):
    """Test various stuff from the pvutils module."""

    def setUp(self):
        """
        This method is called before each test and resets the ParaView state.
        """

        # Set default values for global parameters.
        pvutils.reset_paraview()

    def _get_test_name(self):
        """
        Return the name of the current test method, without the leading test_.
        """

        name = unittest.TestCase.id(self)
        split_name = name.split('.')[-1].split('_', 1)
        if not split_name[0] == 'test':
            raise ValueError(('The test name {} does not match the expected '
                + 'format').format(name))
        return split_name[1]

    def test_full_script_beam_to_solid_volume_meshtying(self):
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
        scene = pa.GetAnimationScene()
        scene.UpdateAnimationUsingDataTimeSteps()
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
        color_function = pa.GetColorTransferFunction('nodal_2PK_stresses_xyz')
        color_bar = pa.GetScalarBar(color_function, view)
        color_bar.Title = '$S_{ZZ}$'
        color_bar.ComponentTitle = ''
        color_bar.WindowLocation = 'AnyLocation'
        color_bar.Position = [0.69, 0.6]
        color_bar.ScalarBarLength = 0.2
        pvutils.set_colorbar_font(color_bar, font_size, dpi, font='TeX')

        # Set and place the color map for the beam.
        display = pa.GetDisplayProperties(beam, view=view)
        display.SetScalarBarVisibility(view, True)
        color_function = pa.GetColorTransferFunction('curvature_2_GPs')
        color_bar = pa.GetScalarBar(color_function, view)
        color_bar.Title = '$\kappa$'
        color_bar.ComponentTitle = ''
        color_bar.WindowLocation = 'AnyLocation'
        color_bar.Position = [0.69, 0.2]
        color_bar.ScalarBarLength = 0.2
        pvutils.set_colorbar_font(color_bar, font_size, dpi, font='TeX')

        # This part will only be executed on local jobs, since GitLab can not
        # open a X window.
        if not _is_gitlab():
            # Export screenshot.
            screenshot_path = os.path.join(testing_temp, '{}_temp.png'.format(
                self._get_test_name()))
            pa.SaveScreenshot(screenshot_path,
                view,
                ImageResolution=size_pixel,
                OverrideColorPalette='WhiteBackground',
                TransparentBackground=0,
                FontScaling='Do not scale fonts'
                )

            # Compare the created image with the reference image by comparing
            # their hashes. The image created with ParaView is slightly
            # different from the one created with pvpython. The eyeball test
            # passes, therefore two reference images will be compared.
            if pvutils.is_pvpython():
                interpreter = 'pvpython'
            else:
                interpreter = 'paraview'
            self.assertEqual(
                _get_file_hash(screenshot_path),
                _get_file_hash(os.path.join(testing_reference,
                    '{}_ref_{}.png'.format(
                        self._get_test_name(), interpreter))))


if __name__ == '__main__':
    # Execution part of script.

    # Define the testing paths.
    testing_path = os.path.abspath(os.getcwd())
    testing_reference = os.path.join(testing_path, 'reference-files')
    testing_temp = os.path.join(testing_path, 'testing-tmp')

    # Check and clean the temp directory.
    if not os.path.isdir(testing_temp):
        os.makedirs(testing_temp)
    _empty_temp_testing_directory()

    # Perform the tests.
    unittest.TextTestRunner().run(
        unittest.TestLoader().loadTestsFromTestCase(TestPvutils))
