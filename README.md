# pvutils

A set of functions that should simplify the usability of `python` scripts in ParaView.

## Setup

Three different ways exist to execute a python ParaView script

- `<path to ParaView>/bin/pvpython <path to script>`: execute the script with the ParaView internal python interpreter. Per default no grapical window is opened.
- `<path to ParaView>/bin/paraview --script=<path to script>`: execute the script with the default ParaView application. A GUI will be opened.
- `python2 <path to script>`: execute the script with the system python interpreter. For this to work, certain path have to be in the environment variables. For details regarding the paths, have a look in the files `.gitlab.yml` and `tests/test_local.sh`.

If the script uses `pvutils` the path to the root `pvutils` directory has to be in `PYTHONPATH`.


## Testing

The unit tests are defined in `tests/testing.pvutils.py`.
All tests can be t with the script `tests/test_local.sh`. The environment variable `PARAVIEW_PATH` - pointing to the root directory of the ParaView installation (the one containing `bin` and `lib`), has to be set.
The path to `pvutils` has to be in `PYTHONPATH`.
A change in the repository also starts a pipeline in GitLab.
Tests that require a graphical output (e.g. everything that writes an image) can not be performed via GitLab, therefore a local run of `tests/test_local.sh` should be performed before EVERY commit.
