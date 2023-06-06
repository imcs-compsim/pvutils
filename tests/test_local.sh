#!/bin/bash

# For this script to work, the variable PARAVIEW_PATH has to be set and python has to find the pvutils package.

# Set paths and variables.
ORIGINAL_PYTHONPATH=$PYTHONPATH
ORIGINAL_LD_LIBRARY_PATH=$LD_LIBRARY_PATH
PVUTILS_DIR="$(dirname $(pwd))"

# Test with pvpython.
echo "Test pvutils with pvpython"
export PYTHONPATH="${ORIGINAL_PYTHONPATH}:${PVUTILS_DIR}"
export LD_LIBRARY_PATH="${ORIGINAL_LD_LIBRARY_PATH}"
${PARAVIEW_PATH}/bin/pvpython testing_pvutils.py

# Test with system python.
echo ""
if [ -z ${PYTHON_EXE+x} ];
then
    echo "Python executable is not set, tests will be skipped"
else
    echo "Test pvutils with system python"
    # Don't set the pvutils path here, since we assume pvutils is installed in the provied python path
    export PYTHONPATH="${ORIGINAL_PYTHONPATH}:${PARAVIEW_PATH}/lib/python3.9/site-packages"
    export LD_LIBRARY_PATH="${ORIGINAL_LD_LIBRARY_PATH}:${PARAVIEW_PATH}/lib"
    ${PYTHON_EXE} testing_pvutils.py
fi

# Test with paraview. This will open a GUI that has to be closed by the user.
echo ""
echo "Test pvutils with paraview (will open GUI)"
export PYTHONPATH="${ORIGINAL_PYTHONPATH}:${PVUTILS_DIR}"
export LD_LIBRARY_PATH="${ORIGINAL_LD_LIBRARY_PATH}"
${PARAVIEW_PATH}/bin/paraview --script=testing_pvutils.py
