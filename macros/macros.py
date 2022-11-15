# -*- coding: utf-8 -*-
"""
Print the state of the current view in ParaView.
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
