#!/bin/bash

# For this script to work, the variable PARAVIEW_PATH has to be set and python has to find the pvutils package.

# Set paths and variables.
ORIGINAL_PYTHONPATH=$PYTHONPATH
ORIGINAL_LD_LIBRARY_PATH=$LD_LIBRARY_PATH

# Test with pvpython.
echo "Test pvutils with pvpython"
export PYTHONPATH="${ORIGINAL_PYTHONPATH}"
export LD_LIBRARY_PATH="${ORIGINAL_LD_LIBRARY_PATH}"
${PARAVIEW_PATH}/bin/pvpython testing_pvutils.py

# Test with system python.
echo ""
if [ -z ${PYTHON_EXE+x} ];
then
    echo "Python executable is not set, tests will be skipped"
else
    echo "Test pvutils with system python"
    export PYTHONPATH="${ORIGINAL_PYTHONPATH}:${PARAVIEW_PATH}/lib/python3.9/site-packages"
    export LD_LIBRARY_PATH="${ORIGINAL_LD_LIBRARY_PATH}:${PARAVIEW_PATH}/lib"
    ${PYTHON_EXE} testing_pvutils.py
fi

# Test with paraview. This will open a GUI that has to be closed by the user.
echo ""
echo "Test pvutils with paraview (will open GUI)"
export PYTHONPATH="${ORIGINAL_PYTHONPATH}"
export LD_LIBRARY_PATH="${ORIGINAL_LD_LIBRARY_PATH}"
${PARAVIEW_PATH}/bin/paraview --script=testing_pvutils.py
