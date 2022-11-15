# -*- coding: utf-8 -*-
"""
Functions to convert meshes to a formats that can be imported in blender.
"""

# Import python modules.
import os
import numpy as np

# Load the paraview utility functions.
import pvutils

# Import paraview module.
import paraview.simple as pa


class ThirdOrderCurve(object):
    """
    This class helps to evaluate the tangent on a third order curve that is
    only given by discrete points (4 are needed).
    """

    def __init__(self, n_segments=5):
        """
        Setup the class and the interpolation variables.
        """

        # Parameter positions of the points that will be given.
        self.xi = np.array([-1, -1 + 2.0 / n_segments, 1 - 2.0 / n_segments, 1])

        # Interpolation matrix.
        A = np.zeros([4, 4])
        for i_row in range(4):
            for i_col in range(4):
                A[i_row, i_col] = self.xi[i_row] ** i_col

        # Inverted interpolation matrix.
        self.A_inv = np.linalg.inv(A)

        # Coefficient matrix.
        self.coefficiens = np.zeros([4, 3])

    def get_tangent(self, cell_positions, xi):
        """
        Calculate the tangent.
        """

        for i_dir in range(3):
            self.coefficiens[:, i_dir] = np.dot(self.A_inv, cell_positions[:, i_dir])

        temp = np.zeros(3)
        for i_coeff in range(1, 4):
            temp += self.coefficiens[i_coeff, :] * xi ** (i_coeff - 1) * i_coeff
        return temp


def surface_to_blender(
    item,
    name,
    base_dir,
    time_steps=None,
    frames=None,
    frame_interval=1,
    set_time_step_function=None,
):
    """
    Convert a surface mesh to blender. The topology of the mesh has to stay
    constant for all time steps. It is assumed, that the source has already an
    applied warp filter, i.e. the point coordinates are already at the current
    position.

    Note: To convert a volume to its boundary surfaces the following ParaView
        filters have to be applied:
            - CleantoGrid -> (optional for BACI output) Removes double nodes,
                so the volume cells are connected to each other.
            - ExtractSurface -> Extracts the boundaries of the volume, i.e.,
                the surfaces.
    """

    # If no time steps are given, use all available.
    if time_steps is None:
        time_steps = pvutils.get_available_timesteps()

    # Save the mesh to disk. The displacement at this stage does not matter,
    # as it will be overwritten in blender.
    pa.SaveData(
        os.path.join(base_dir, name + "_mesh.ply"), proxy=item, FileType="Binary"
    )

    # Get group data.
    vtk_data = pvutils.get_vtk_data_as_numpy(
        item, cell_data=True, cell_connectivity=True
    )
    cell_data = vtk_data["cell_data"]
    cell_connectivity = vtk_data["cell_connectivity"]
    blender_cell_group = {}
    if "blender_cell_group" in cell_data.keys():

        # Get all vertices for all grouped cells.
        for i_cell, cell_group in enumerate(cell_data["blender_cell_group"]):
            cell_group_int = int(cell_group)
            if cell_group_int not in blender_cell_group.keys():
                blender_cell_group[cell_group_int] = []
            blender_cell_group[cell_group_int].extend(cell_connectivity[i_cell])

        # Extract the unique vertices.
        for group_id in blender_cell_group.keys():
            blender_cell_group[group_id] = list(set(blender_cell_group[group_id]))

    for i_time, time in enumerate(time_steps):
        if set_time_step_function is None:
            pvutils.set_timestep(time)
        else:
            set_time_step_function(time)
        pvutils.get_view().Update()
        vtk_data = pvutils.get_vtk_data_as_numpy(item, coordinates=True)

        if i_time == 0:
            position = np.zeros([len(time_steps), len(vtk_data["coordinates"]), 3])

        position[i_time, :] = vtk_data["coordinates"]
        print("Got surface coordinates for {} at time {}".format(name, time))

    if frames is not None:
        if not len(frames) == len(time_steps):
            raise ValueError("Frames must have the same length as time steps")
    else:
        frames = [1 + i * frame_interval for i in range(len(time_steps))]

    # Save data.
    np.save(
        os.path.join(base_dir, name + "_data"),
        {
            "position": position,
            "frames": frames,
            "blender_cell_group": blender_cell_group,
        },
    )


def fibers_to_blender(
    item,
    file_name,
    time_steps=None,
    n_segments=5,
    frames=None,
    frame_interval=1,
    set_time_step_function=None,
):
    """
    Convert a fiber representation in ParaView (with polylines) to blender.

    Args
    ----
    item: ParaView source
        Item that sould be converted to blender. The source can only contain
        polyline cells.
    file_name: str
        File to be used for output files.
    n_segments: int
        Number of visualization segments that are used for each beam finite
        element.
    frames: [int]
        List of frames. If this is not given, then the frame interval will be
        used.
    frame_interval: int
        Number of frames between two successive time steps.
    set_time_step_function: Function
        Function to set the time step. This can be used to manually change
        certain frames.
    """

    # If no time steps are given, use all available.
    if time_steps is None:
        time_steps = pvutils.get_available_timesteps()

    # Get the vtk data.
    vtk_data = pvutils.get_vtk_data_as_numpy(
        item, coordinates=True, point_data=True, cell_connectivity=True
    )

    if not np.max(np.abs(np.array(vtk_data["cell_types"]) - 4)) == 0:
        raise ValueError('Only polylines allowed in "fibers_to_blender"')

    # Set the connectivity of the vtk data.
    vtk_connectivity = vtk_data["cell_connectivity"]

    # Get the output data.
    n_fibers = len(vtk_connectivity)

    def get_n_elements_per_fiber(i):
        n_connectivity = len(vtk_connectivity[i])
        if (n_connectivity - 1) % n_segments == 0:
            return len(vtk_connectivity[i]) // n_segments
        else:
            raise ValueError("Number of points and number of segments does not match!")

    n_elements_per_fiber = [get_n_elements_per_fiber(i) for i in range(n_fibers)]
    n_elements_total = np.sum(n_elements_per_fiber)
    n_points_total = n_elements_total + n_fibers
    coordinates = np.zeros([len(time_steps), n_points_total, 3])
    tangents = np.zeros([len(time_steps), n_points_total, 3])

    # Get the topology of the output.
    topology = []
    last_index = 0
    cell_index = 0
    cell_relevant_nodes = np.zeros([n_elements_total, 4], dtype=int)
    cell_to_point = np.zeros([n_elements_total, 2], dtype=int)
    fiber_radii = np.zeros(n_fibers)
    for i_fiber, local_connectivity in enumerate(vtk_connectivity):
        fiber_radii[i_fiber] = vtk_data["point_data"]["cross_section_radius"][
            local_connectivity[0]
        ]

        cell_topology = [
            last_index + i for i in range(n_elements_per_fiber[i_fiber] + 1)
        ]
        topology.append(cell_topology)
        last_index = cell_topology[-1] + 1

        for i_cell in range(n_elements_per_fiber[i_fiber]):
            cell_relevant_nodes[cell_index] = [
                local_connectivity[n_segments * i_cell],
                local_connectivity[n_segments * i_cell + 1],
                local_connectivity[n_segments * (i_cell + 1) - 1],
                local_connectivity[n_segments * (i_cell + 1)],
            ]
            cell_to_point[cell_index, 0] = cell_topology[i_cell]
            if i_cell == n_elements_per_fiber[i_fiber] - 1:
                cell_to_point[cell_index, 1] = cell_topology[-1]
            else:
                cell_to_point[cell_index, 1] = -1
            cell_index += 1

    curve = ThirdOrderCurve()
    for i_time, time in enumerate(time_steps):
        if set_time_step_function is None:
            pvutils.set_timestep(time)
        else:
            set_time_step_function(time)
        pvutils.get_view().Update()
        vtk_data = pvutils.get_vtk_data_as_numpy(item, coordinates=True)

        for i_cell in range(n_elements_total):

            # Get the relevant positions for this cell.
            cell_positions = vtk_data["coordinates"][cell_relevant_nodes[i_cell]]

            # Set the position and tangent function.
            for i_start_end, i_local, xi_local in [[0, 0, -1.0], [1, -1, 1.0]]:
                temp = cell_to_point[i_cell, i_start_end]
                if temp >= 0:
                    coordinates[i_time, temp] = cell_positions[i_local]
                    # The scaling is to match the tangent definition in
                    # blender.
                    tangents[i_time, temp] = (
                        curve.get_tangent(cell_positions, xi_local) * 2.0 / 3.0
                    )

        print("Got fiber coordinates at time {}".format(time))

    if frames is not None:
        if not len(frames) == len(time_steps):
            raise ValueError("Frames must have the same length as time steps")
    else:
        frames = [1 + i * frame_interval for i in range(len(time_steps))]

    np.save(
        file_name,
        {
            "connectivity": topology,
            "coordinates": coordinates,
            "tangents": tangents,
            "fiber_radii": fiber_radii,
            "frames": frames,
        },
    )
