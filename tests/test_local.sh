#!/bin/bash

# For this script to work, the variable PARAVIEW_PATH has to be set and python has to find the pvutils package.

# Set paths and variables.
ORIGINAL_PYTHONPATH=$PYTHONPATH
ORIGINAL_LD_LIBRARY_PATH=$LD_LIBRARY_PATH

# Test with pvpython.
echo "Test pvutils with pvpython"
export PYTHONPATH="${ORIGINAL_PYTHONPATH}"
export LD_LIBRARY_PATH="${ORIGINAL_LD_LIBRARY_PATH}"
export IS_PVPYTHON=1
${PARAVIEW_PATH}/bin/pvpython testing_pvutils.py

# Test with system python.
echo ""
echo "Test pvutils with system python"
export PYTHONPATH="${ORIGINAL_PYTHONPATH}:${PARAVIEW_PATH}/lib/python2.7/site-packages"
export LD_LIBRARY_PATH="${ORIGINAL_LD_LIBRARY_PATH}:${PARAVIEW_PATH}/lib"
export IS_PVPYTHON=1
python2 testing_pvutils.py

# Test with paraview. This will open a GUI that has to be closed by the user.
echo ""
echo "Test pvutils with paraview (will open GUI)"
export PYTHONPATH="${ORIGINAL_PYTHONPATH}"
export LD_LIBRARY_PATH="${ORIGINAL_LD_LIBRARY_PATH}"
export IS_PVPYTHON=0
${PARAVIEW_PATH}/bin/paraview --script=testing_pvutils.py
