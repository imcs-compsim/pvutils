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
Macros for the ParaView application.
"""


# Import paraview module.
import paraview.simple as pa
from pvutils import print_view_state, get_view


def macro_print_view_state():
    """Print the state of the view in the ParaView console."""
    print_view_state(get_view())


def macro_reload_sources():
    """Reload all sources that depend on a file."""

    # Get all sources in the current session.
    sources = pa.GetSources()
    for _key, value in sources.items():
        if (
            "FileName" in value.ListProperties()
            or "FileNames" in value.ListProperties()
        ):
            # If a source is linked to a file, reload it.
            return_value = pa.ReloadFiles(value)
            if not return_value:
                raise ValueError("Error in reload file script!")
