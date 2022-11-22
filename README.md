# pvutils

A set of functions that should simplify the usability of `python` scripts in ParaView.

- The directory `pvuitls` contains the files for this module, utility functions and ParaView filter wrappers. 
- Some useful programmable filters can be found in the directory `filters`.
To use one of the programmable filters in a script, use the function `programmable_filter` from `pvutils`.


## ParaView version

Unfortunately the ParaView scripting interface changes with minor version.
The current version of this module is developed with [ParaView 5.10.1](https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.10&type=binary&os=Linux&downloadFile=ParaView-5.10.1-MPI-Linux-Python3.9-x86_64.tar.gz).


## Execute scripts with `pvutils`

Three two ways exist to execute a python ParaView script

- Execute with the ParaView python interpreter (make sure that the path to the root `pvutils` directory is in `PYTHONPATH`):
  - `<path to ParaView>/bin/pvpython <path to script>`: execute the script with the ParaView internal python interpreter.
  Per default no graphical window is opened.
  - `<path to ParaView>/bin/paraview --script=<path to script>`: execute the script with the default ParaView application.
  A GUI will be opened.

- Execute with the system python interpreter `python3 <path to script>`.
This the recommended way to add the ParaView script to an existing python workflow.
For this to work, certain things have to be considered:
  - Add the path to the ParaView python interface
    ```
    export PYTHONPATH="${PYTHONPATH}:${<path to ParaView>}/lib/python3.9/site-packages"
    ```
  - Add the path to the ParaView python libraries
    ```
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${<path to ParaView>}/lib"
    ```
  - If your external python workflow depends on binary packages, e.g., `cython` or `numpy`, make sure the libraries provided by ParaView are compatible with your binaries.


## Testing

A change in the repository also starts a pipeline in GitLab.

The unit tests are defined in `tests/testing.pvutils.py`.

To setup local testing, modify your `~/.bashrc` to
- point the environment variable `PARAVIEW_PATH` pointing to the root directory of the ParaView installation (the one containing `bin` and `lib`)
- include th repository `pvutils` into your `PYTHONPATH`.

Then, all tests can be executed by
- navigating to `tests/`
- executing `test_local.sh` 
