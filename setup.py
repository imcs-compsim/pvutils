# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# PVUTILS: Utility tools for the ParaView python interface
#
# Copyright 2023 Ivo Steinbrecher & Matthias Mayr
#                Institute for Mathematics and Computer-based Simulation
#                University of the Bundeswehr Munich
#                https://www.unibw.de/imcs-en
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

from setuptools import setup

setup(
    name="pvutils",
    version="0.1",
    author="Matthias Mayr, Ivo Steinbrecher",
    author_email="matthias.mayr@unibw.de, ivo.steinbrecher@unibw.de",
    description="Utils for ParaView",
    packages=["blender", "pvutils"],
    install_requires=[
        "numpy==1.21.1",  # This is the version provided with ParaView 5.10.1
    ],
)
