"""
Use this script as a template for ParaView scripts.
ParaView 5.6 was used to develop this package.
"""


# Load the ParaView utility functions. The path to this module has to be
# added to the ParaView python system paths. Alternatively this can be done by
# setting the PYTHONPATH environment variable.
# export PYTHONPATH="${PYTHONPATH}:/home/ivo/dev/python/pvutils"
import sys
print(sys.path)
sys.path.append('/home/ivo/dev/python/pvutils')

# Load the paraview utility functions.
import pvutils.utils as pvutils

# Import paraview module.
import paraview.simple as pa
pa._DisableFirstRenderCameraReset()
