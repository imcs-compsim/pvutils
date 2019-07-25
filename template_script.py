"""
Use this script as a template for paraview scripts.
ParaView 5.6 was used to develope this package.
"""


# Load the paraview utility functions. The path to this module has to be
# added to the paraview python system paths. Load the paraview utility
# functions.
import sys
sys.path.append('/home/ivo/dev/python/pvutils')
import pvutils.utils as pvutils

# Import paraview module.
import paraview.simple as pa
pa._DisableFirstRenderCameraReset()
