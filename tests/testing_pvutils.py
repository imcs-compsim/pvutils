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


def compare_data(data1, data2, raise_error=False, tol_float=None):
    """
    Compare the ParaView data1 and data2, by compairing the stored data.

    Args
    ----
    raise_error: bool
        If true, then an error will be raised in case the files do not match.
        Otherwise False will be returned.
    tol_float: None / float
        If given, numbers will be considered equal if the difference between
        them is smaller than tol_float.
    """

    # Default value for the numerical tolerance.
    if tol_float is None:
        tol_float = 1e-16

    def compare_arrays(array1, array2, name=None):
        """
        Compare two vtk arrays.
        """

        diff = vtk_numpy.vtk_to_numpy(array1) - vtk_numpy.vtk_to_numpy(array2)
        if not np.max(np.abs(diff)) < tol_float:
            error_string = 'VTK array comparison failed!'
            if name is not None:
                error_string += ' Name of the array: {}'.format(name)
            raise ValueError(error_string)

    def compare_data_sets(data1, data2):
        """
        Compare data sets obtained from vtk objects.
        """

        # Both data sets need to have the same number of arrays.
        if not data1.GetNumberOfArrays() == data2.GetNumberOfArrays():
            raise ValueError('Length of vtk data objects do not match!')

        # Compare each array.
        for i in range(data1.GetNumberOfArrays()):

            # Get the arrays with the same name.
            name = data1.GetArrayName(i)
            array1 = data1.GetArray(name)
            array2 = data2.GetArray(name)
            compare_arrays(array1, array2, name=name)

    # Perform all checks, catch errors.
    try:

        # Compare the point positions.
        compare_arrays(
            data1.GetPoints().GetData(),
            data2.GetPoints().GetData(),
            name='point_positions'
            )

        # Compare the cell and point data of the array.
        compare_data_sets(data1.GetCellData(), data2.GetCellData())
        compare_data_sets(data1.GetPointData(), data2.GetPointData())

        # Compare the cell connectivity.
        compare_arrays(
            data1.GetCells().GetData(),
            data2.GetCells().GetData(),
            name='cell_connectivity')

        # Compare the cell types.
        compare_arrays(
            data1.GetCellTypesArray(),
            data2.GetCellTypesArray(),
            name='cell_type')

    except Exception as error:
        if raise_error:
            raise error
        return False

    return True


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

    def _save_screenshot_and_compare(self, view, **kwargs):
        """
        Save the view to a screenshot and compare it with the reference image.
        The reference image is stored in the 'reference=files' directory. When
        the testing is performed via GitLab, this part is skipped, as currently
        we can not open a X window and therefore create screenshots with the
        GitLab runner.

        Args
        ----
        view: ParaView view object
            View that will be written to an image.
        kwargs:
            Will be passed to pa.SaveScreenshot.
        """

        if not _is_gitlab():

            from matplotlib.image import imread

            # Export screenshot.
            screenshot_path = os.path.join(testing_temp,
                '{}_temp.png'.format(self._get_test_name()))
            pa.SaveScreenshot(screenshot_path, view, **kwargs)

            # Compare the created image with the reference image.
            reference_path = os.path.join(testing_reference,
                 '{}_ref.png'.format(self._get_test_name()))
            test_image = imread(screenshot_path)
            ref_image = imread(reference_path)

            # Get the average difference.
            # https://www.pyimagesearch.com/2014/09/15/python-compare-two-images/
            err = np.sum((test_image.astype("float")
                - ref_image.astype("float")) ** 2)
            err /= float(test_image.shape[0] * test_image.shape[1])

            self.assertTrue(err < 5e-4)

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

        # Test that the correct name can be read from the source.
        self.assertEqual('solid', pvutils.get_source_name(solid))

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

        # Compare the current view with the reference image.
        self._save_screenshot_and_compare(view,
                ImageResolution=size_pixel,
                OverrideColorPalette='WhiteBackground',
                TransparentBackground=0,
                FontScaling='Do not scale fonts'
                )

    def test_time_step(self):
        """
        Test the set time step functions. Two meshes are loaded with different
        time values. Both meshes contain a solid block with the same (constant
        per time step) displacement.
        """

        three_steps = pvutils.load_file(os.path.join(testing_reference,
            'solid_cantilever_pvd', 'solid_cantilever_3-steps.pvd'))
        four_steps = pvutils.load_file(os.path.join(testing_reference,
            'solid_cantilever_pvd', 'solid_cantilever_4-steps.pvd'))

        # Check if all time steps are found.
        time_step_error = np.linalg.norm(
            np.array(pvutils.get_available_timesteps())
            - [0.3, 0.4, 0.6, 0.8, 0.9, 1.2]
            )
        self.assertTrue(time_step_error < 1e-10)

        # The rest will only be tested locally since it requires an X window
        # (for the temporal time filter).
        if not _is_gitlab():
            # Apply filters to the meshes and display them. The temporal
            # interpolator is applied, so both meshes have the same
            # displacement (interpolated).
            three_steps = pvutils.temporal_interpolator(three_steps)
            three_steps = pvutils.transform(three_steps, translate=[0, 0, 2])
            three_steps = pvutils.warp(three_steps)
            pvutils.display(three_steps)
            four_steps = pvutils.warp(four_steps)
            pvutils.display(four_steps)

            # Set the time. The three step mesh will be interpolated, as it
            # does not have a discrete time step here.
            pvutils.set_timestep(0.9, fail_on_not_available_time=True)

            # Set the view.
            view = pa.GetActiveViewOrCreate('RenderView')
            view.CameraPosition = [11.429, 2.4, 1]
            view.CameraFocalPoint = [0, 2.4, 1]
            view.CameraViewUp = [0, 0, 1]
            view.CameraViewAngle = 30
            view.CameraParallelScale = 2.95804
            view.OrientationAxesVisibility = 0
            view.CameraParallelProjection = 0
            view.InteractionMode = '2D'
            view.ViewSize = [400, 400]

            # Compare the current view with the reference image.
            self._save_screenshot_and_compare(view,
                    OverrideColorPalette='WhiteBackground',
                    TransparentBackground=0
                    )

    def test_beam_display(self):
        """
        Test the BeamDisplay object.
        """

        # Load the beam.
        beam = pvutils.BeamDisplay(os.path.join(testing_reference,
            'beam_cantilever_pvd', 'cantilever.pvd'), segments=30,
            factor_nodes=4, factor_triads=10)

        # Set the view.
        view = pvutils.get_view()
        view.CameraPosition = [-0.907414, 1.60167, 1.27251]
        view.CameraFocalPoint = [0.161918, -0.0745694, 0.503964]
        view.CameraViewUp = [0.208326, -0.294696, 0.932606]
        view.CameraViewAngle = 30
        view.CameraParallelScale = 1.73
        view.OrientationAxesVisibility = 1
        view.CameraParallelProjection = 0
        view.ViewSize = [600, 800]
        view.InteractionMode = '3D'

        # Color the beam parts.
        pvutils.contour(beam.beam_tube)
        pvutils.contour(beam.nodes)
        for triad in beam.base_vectors:
            pvutils.contour(triad)

        # Compare the current view with the reference image.
        self._save_screenshot_and_compare(view,
                OverrideColorPalette='WhiteBackground',
                TransparentBackground=0
                )

    def test_get_parents(self):
        """
        Test the get_parents function.
        """

        # Load a beam representation and get the parents of some entries.
        beam = pvutils.BeamDisplay(os.path.join(testing_reference,
            'beam_cantilever_pvd', 'cantilever.pvd'), segments=30,
            factor_nodes=4, factor_triads=10)
        node_parents = pvutils.get_parents(beam.nodes)
        tube_parents = pvutils.get_parents(beam.beam_tube)

        # Compare with the expected parents.
        self.assertEqual(node_parents,
            [beam.nodes, beam.endpoints, beam.beam_cell_to_point, beam.beam])
        self.assertEqual(tube_parents,
            [
                beam.beam_tube, beam.beam_extract_surface,
                beam.beam_cell_to_point, beam.beam
            ])

    def test_script_load_beam_to_solid(self):
        """
        Load a patch test of a solid structure with beams, apply some filters
        to it and check the created screenshot. This will test if the complete
        workflow works.
        """

        # Call the script.
        _solid, beam = pvutils.scripts.load_beam_to_solid_in_dir(
            os.path.join(testing_reference, 'beam_and_solid_cantilever')
            )

        # Display some features on the beam.
        pvutils.contour(beam.nodes)
        pa.Show(beam.base_vectors[1])

        # Setup the view.
        view = pvutils.get_view()
        view.CameraPosition = [37.0083, 21.8822, -37.7023]
        view.CameraFocalPoint = [4.54484, 3.0232, -2.89383]
        view.CameraViewUp = [-0.0915005, 0.90881, 0.407052]
        view.CameraViewAngle = 30
        view.CameraParallelScale = 16.0335
        view.OrientationAxesVisibility = 0
        view.CameraParallelProjection = 0
        view.ViewSize = [696, 654]
        view.InteractionMode = '3D'

        # Compare the current view with the reference image.
        self._save_screenshot_and_compare(view,
                OverrideColorPalette='WhiteBackground',
                TransparentBackground=0
                )


    def test_merge_polylines_filter(self):
        """
        Test the merge polylines programmable filter.
        """

        raw_file = os.path.join(testing_reference,
            'beam_merge_polylines_raw.vtu')
        ref_file = os.path.join(testing_reference,
            'beam_merge_polylines_reference.vtu')
        test_file = os.path.join(testing_temp, 'beam_merge_polylines.vtu')

        # Load the beam (with the merge polylines filter)
        beam = pvutils.BeamDisplay(raw_file, merge_poly_lines=True)

        # Load the reference file.
        beam_ref = pvutils.load_file(ref_file)

        # Compare the vtk file with the reference file.
        self.assertTrue(compare_data(
            pa.servermanager.Fetch(beam_ref),
            pa.servermanager.Fetch(beam.beam_merge_poly_line),
            raise_error=True))


if __name__ == '__main__':
    # Execution part of script.

    pa._DisableFirstRenderCameraReset()

    # Define the testing paths.
    testing_path = os.path.abspath(os.getcwd())
    testing_reference = os.path.join(testing_path, 'reference-files')
    testing_temp = os.path.join(testing_path, 'testing-tmp')

    # Check and clean the temp directory.
    if not os.path.isdir(testing_temp):
        os.makedirs(testing_temp)
    _empty_temp_testing_directory()

    # Perform the tests.
    run = unittest.TextTestRunner().run(
        unittest.TestLoader().loadTestsFromTestCase(TestPvutils))
    sys.exit(not (run.wasSuccessful() and len(run.skipped) == 0))
