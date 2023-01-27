# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# PVUTILS: Utility tools for the ParaView python interface
#
# Copyright 2023 Ivo Steinbrecher & Matthias Mayr
#                 Institute for Mathematics and Computer-based Simulation
#                 University of the Bundeswehr Munich
#                 https://www.unibw.de/imcs-en
#
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# -----------------------------------------------------------------------------
"""
Create the input file for the beam and solid test case.
"""

# Import mesh utilities.
from meshpy import (
    InputFile,
    mpy,
    MaterialReissner,
    Beam3rHerm2Line3,
    Function,
    BoundaryCondition,
)
from meshpy.mesh_creation_functions import create_beam_mesh_line
from meshpy.header_functions import set_header_static, set_runtime_output
from cubitpy import cupy, CubitPy
from cubitpy.mesh_creation_functions import create_brick

# Problem parameters.
load = 0.01

# Create the cubit cantilever.
cubit = CubitPy()
solid_cantilever = create_brick(
    cubit, 1, 1, 30, mesh_interval=[2, 2, 5], element_type=cupy.element_type.hex27
)
cubit.add_node_set(
    solid_cantilever.surfaces()[1],
    name="solid_fix",
    bc_type=cupy.bc_type.dirichlet,
    bc_description="NUMDOF 3 ONOFF 1 1 1 VAL 0 0 0 FUNCT 0 0 0",
)
cubit.add_node_set(
    solid_cantilever.surfaces()[0],
    name="solid_load",
    bc_type=cupy.bc_type.neumann,
    bc_description="NUMDOF 3 ONOFF 1 1 0 VAL {0} {0} 0 FUNCT 1 1 0".format(load),
)

# Set the head string.
cubit.head = """
    --------------------------------------------------------------MATERIALS
    MAT 1 MAT_Struct_StVenantKirchhoff YOUNG 100 NUE 0.3 DENS 0.0
    """

# Add the beam.
mesh = InputFile(cubit=cubit)

# Create material and time function.
mat = MaterialReissner(radius=0.6, youngs_modulus=100, nu=0.3)
fun = Function("COMPONENT 0 FUNCTION t")
mesh.add(mat, fun)

# Add beam.
beam_set = create_beam_mesh_line(
    mesh, Beam3rHerm2Line3, mat, [4, 0, -15], [4, 0, 15], n_el=5
)
mesh.add(
    BoundaryCondition(
        beam_set["start"],
        (
            "NUMDOF 9 ONOFF 1 1 1 1 1 1 0 0 0 VAL 0 0 0 0 0 0 0 0 0 "
            + "FUNCT 0 0 0 0 0 0 0 0 0"
        ),
        bc_type=mpy.bc.dirichlet,
    )
)
mesh.add(
    BoundaryCondition(
        beam_set["end"],
        (
            "NUMDOF 9 ONOFF 1 1 0 0 0 0 0 0 0 VAL {0} {0} 0 0 0 0 0 0 0 "
            + "FUNCT 1 1 0 0 0 0 0 0 0"
        ),
        format_replacement=[load],
        bc_type=mpy.bc.neumann,
    )
)

# Set headers.
set_header_static(
    mesh, time_step=1.0, n_steps=1, tol_residuum=1e-10, tol_increment=1e-10
)
set_runtime_output(mesh)

# Create the input file.
mesh.write_input_file("/home/ivo/temp/beam_and_solid_cantilever.dat")
