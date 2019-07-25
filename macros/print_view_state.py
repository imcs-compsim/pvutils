# -*- coding: utf-8 -*-
"""
Print the state of the current view in ParaView.
"""


# Import paraview module.
import paraview.simple as pa
from pvutils.utils import print_view_state

view = view = pa.GetActiveViewOrCreate('RenderView')
print_view_state(view)
