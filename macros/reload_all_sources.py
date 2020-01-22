# -*- coding: utf-8 -*-
"""
Reload all sources that depend on a file.
"""

# Import paraview module.
import paraview.simple as pa

# Get all sources in the current session.
sources = pa.GetSources()
for _key, value in sources.items():
    if ('FileName' in value.ListProperties()
            or 'FileNames' in value.ListProperties()):
        # If a source is linked to a file, reload it.
        return_value = pa.ReloadFiles(value)
        if not return_value:
            raise ValueError('Error in reload file script!')
