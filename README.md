# pvutils

A set of functions that should simplify the usability of `python` scripts in ParaView.

- The directory `pvuitls` contains the files for this module, utility functions and ParaView filter wrappers. 
- Some useful programmable filters can be found in the directory `filters`. To use one of the programmable filters in a script, use the function `programmable_filter` from `pvutils`.


## ParaView version

Unfortunately the ParaView scripting interface changes with minor version.
The current version of this module is developed with [ParaView 5.10.1](https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.10&type=binary&os=Linux&downloadFile=ParaView-5.10.1-MPI-Linux-Python3.9-x86_64.tar.gz).


## Execute scripts with `pvutils`

Three different ways exist to execute a python ParaView script

- `<path to ParaView>/bin/pvpython <path to script>`: execute the script with the ParaView internal python interpreter. Per default no grapical window is opened.
- `<path to ParaView>/bin/paraview --script=<path to script>`: execute the script with the default ParaView application. A GUI will be opened.
- `python2 <path to script>`: execute the script with the system python interpreter. For this to work, certain path have to be in the environment variables. For details regarding the paths, have a look in the files `.gitlab.yml` and `tests/test_local.sh`.

If the script uses `pvutils` the path to the root `pvutils` directory has to be in `PYTHONPATH`.


## Testing

A change in the repository also starts a pipeline in GitLab.
Tests that require a graphical output (e.g. everything that writes an image) can not be performed via GitLab, therefore a local run of `tests/test_local.sh` should be performed before EVERY commit.

The unit tests are defined in `tests/testing.pvutils.py`.

To setup local testing, modify your `~/.bashrc` to
- point the environment variable `PARAVIEW_PATH` pointing to the root directory of the ParaView installation (the one containing `bin` and `lib`)
- include th repository `pvutils` into your `PYTHONPATH`.

Then, all tests can be executed by
- navigating to `tests/`
- executing `test_local.sh` 
