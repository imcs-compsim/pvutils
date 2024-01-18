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
Load stored data in blender.
"""

import bpy
import numpy as np
import os


def set_keyframe_sequence(shape_key, frames, i_frame):
    """
    Set the values for a shape key value at the relevant frames.
    """

    attribute = "value"

    if i_frame == 0:
        return
    else:
        if i_frame == len(frames) - 1:
            values = [0.0, 1.0]
            values_frames = [frames[i_frame + i_inner] for i_inner in [-1, 0]]
        else:
            values = [0.0, 1.0, 0.0]
            values_frames = [frames[i_frame + i_inner] for i_inner in [-1, 0, 1]]

    for v, f in zip(values, values_frames):
        setattr(shape_key, attribute, v)
        shape_key.keyframe_insert(data_path=attribute, frame=f)


def set_fiber_position(shape_key, nodes, coordinates, tangents):
    """
    Set the position for a single fiber.
    """

    # Get all points and tangents we need for this fiber.
    relevant_coordinates = coordinates[nodes]
    relevant_tangents = tangents[nodes]

    # Set the node positions and tangents.
    data = shape_key.data
    for i_node in range(len(nodes)):
        pos = relevant_coordinates[i_node]
        tan = relevant_tangents[i_node]
        data[i_node].co = pos
        if i_node > 0:
            data[i_node].handle_left = pos - tan
        if i_node < len(nodes) - 1:
            data[i_node].handle_right = pos + tan


def create_fiber(nodes, pv_data, i_fiber):
    """
    Create a single fiber.
    """

    # Add cross-section.
    bpy.ops.curve.primitive_bezier_circle_add(
        radius=pv_data["fiber_radii"][i_fiber],
        enter_editmode=False,
        location=(0, 0, 0),
    )
    cs = bpy.context.active_object
    cs.name = "cross_section_{}".format(i_fiber + 1)

    # Create curve.
    bpy.ops.curve.primitive_bezier_curve_add(enter_editmode=False)
    curve = bpy.context.active_object
    curve.name = "beam_{}".format(i_fiber + 1)

    # Get the spline and the point collection.
    spline = curve.data.splines[0]
    points = spline.bezier_points

    # Get number of curve Bezier points.
    n_points = len(nodes)

    # Set the number of nodes and free the first and last one.
    points.add(n_points - 2)
    points[0].handle_left_type = "FREE"
    points[-1].handle_right_type = "FREE"

    # Add shape keys for each time step.
    obj = bpy.context.active_object
    frames = pv_data["frames"]
    for i_frame in range(len(frames)):
        if i_frame == 0:
            name = "Basis"
        else:
            name = f"frame_{i_frame}"
        shape_key = obj.shape_key_add(name=name, from_mix=False)

        # Set the position for this shape key.
        set_fiber_position(
            shape_key,
            nodes,
            pv_data["coordinates"][i_frame],
            pv_data["tangents"][i_frame],
        )

        # Set the key frame sequence.
        set_keyframe_sequence(shape_key, frames, i_frame)

    # Set the tube.
    bpy.context.object.data.bevel_mode = "OBJECT"
    bpy.context.object.data.bevel_object = bpy.data.objects[cs.name]


def load_fibers(file_name):
    """
    Load fibers stored on disk.
    """

    # We have to set the cursor to the origin, so the fiber placement is
    # correct. Save the current location and reset it at the end of this
    # function.
    current_cursor_location = bpy.context.scene.cursor.location
    bpy.context.scene.cursor.location = (0.0, 0.0, 0.0)

    # Load the data written from ParaView.
    pv_data = np.load(file_name, allow_pickle=True, encoding="bytes")[()]

    # Loop over beams.
    connectivity = pv_data["connectivity"]
    for i_fiber, nodes in enumerate(connectivity):
        create_fiber(nodes, pv_data, i_fiber)
        print("Finished fiber {}/{}".format(i_fiber + 1, len(connectivity)))

    # Reset the cursor position.
    bpy.context.scene.cursor.location = current_cursor_location


def load_surface(name, base_dir):
    """
    Load an animated surface in blender.
    """

    # Load the the data.
    bpy.ops.import_mesh.ply(filepath=os.path.join(base_dir, name + "_mesh.ply"))
    obj = bpy.context.active_object
    obj.name = name
    pv_data = np.load(
        os.path.join(base_dir, name + "_data.npy"), allow_pickle=True, encoding="bytes"
    )[()]
    frames = pv_data["frames"]
    position = pv_data["position"]

    # Add shape keys for each time step.
    for i_frame in range(len(frames)):
        if i_frame == 0:
            frame_name = "Basis"
        else:
            frame_name = f"frame_{i_frame}"
        shape_key = obj.shape_key_add(name=frame_name, from_mix=False)

        # Set position for this shape key.
        for i_vert, vert in enumerate(shape_key.data):
            vert.co = position[i_frame, i_vert]

        # Set the keyframe values.
        set_keyframe_sequence(shape_key, frames, i_frame)

        print(
            'Finished surface "{}" for frame {}/{}'.format(
                name, i_frame + 1, len(frames)
            )
        )

    # Set the groups.
    blender_cell_group = pv_data["blender_cell_group"]
    for key in blender_cell_group.keys():
        new_vertex_group = obj.vertex_groups.new(name="cell_group_" + str(key))
        vertex_group_data = blender_cell_group[key]
        new_vertex_group.add(vertex_group_data, 1.0, "ADD")


if __name__ == "__main__":
    """Execution part of script"""

    load_fibers("path_to_the_fiber_file.npy")
    load_surface("name_of_the_surface_in_blender", "path_to_the_surf_dir")
