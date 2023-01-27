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
Define the functions and wrappers that should be available outside of this
module.
"""

# Load utility functions.
from .utility_functions import (
    add_coordinate_axes,
    check_data,
    contour,
    display,
    export_to_tikz,
    get_available_timesteps,
    get_base,
    get_bounding_box,
    get_display,
    get_field_names,
    get_function_arguments,
    get_parents,
    get_size_pixel,
    get_source_name,
    get_view,
    get_vtk_data_as_numpy,
    is_pvpython,
    list_to_mathematica_string,
    load_file,
    print_view_state,
    programmable_filter,
    programmable_source,
    reset_layout,
    reset_paraview,
    reset_print_view_state_color_bar,
    set_categorized_colorbar,
    set_color_range,
    set_colorbar_font,
    set_function_arguments,
    set_print_view_state_color_bar,
    set_timestep,
    setup_view,
    update_scene,
    von_mises_stress,
)

# Load the filter wrappers.
from .filter_wrapper import (
    clip,
    glyph,
    temporal_interpolator,
    threshold,
    transform,
    tube,
    warp,
)

# Load utility classes.
from .utility_classes import BeamDisplay

# Load scripts.
from . import scripts
