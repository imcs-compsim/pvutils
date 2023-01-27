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
"""
This script is used to test that all headers in the repository are correct.

This file is adapted from LaTeX2AI (https://github.com/stoani89/LaTeX2AI).
"""

# Import python modules.
import os
import subprocess
import unittest


def get_repository_dir():
    """
    Get the root directory of this repository.
    """

    script_path = os.path.realpath(__file__)
    root_dir = os.path.dirname(os.path.dirname(script_path))
    return root_dir


def get_license_text():
    """
    Return the license text as a string.
    """

    return """PVUTILS: Utility tools for the ParaView python interface

Copyright 2023 Ivo Steinbrecher & Matthias Mayr
               Institute for Mathematics and Computer-based Simulation
               University of the Bundeswehr Munich
               https://www.unibw.de/imcs-en


Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License."""


def get_all_source_files():
    """
    Get all source files that should be checked for license headers.
    """

    # Get the files in the git repository.
    repo_dir = get_repository_dir()
    process = subprocess.Popen(
        ["git", "ls-files"], stdout=subprocess.PIPE, cwd=repo_dir
    )
    out, _err = process.communicate()
    files = out.decode("UTF-8").strip().split("\n")

    source_line_endings = [".py"]
    source_ending_types = {".py": "py"}
    source_files = {"py": []}
    for file in files:
        extension = os.path.splitext(file)[1]
        if extension not in source_line_endings:
            pass
        else:
            source_files[source_ending_types[extension]].append(
                os.path.join(repo_dir, file)
            )
    return source_files


def license_to_source(license_text, source_type):
    """
    Convert the license text to a text that can be written to source code.
    """

    header = None
    start_line = "-" * 77
    if source_type == "py":
        header = "# -*- coding: utf-8 -*-"
        comment = "#"
    else:
        raise ValueError("Wrong extension!")

    source = []
    if header is not None:
        source.append(header)
    source.append(comment + " " + start_line)
    for line in license_text.split("\n"):
        if len(line) > 0:
            source.append(comment + " " + line)
        else:
            source.append(comment + line)
    source.append(comment + " " + start_line)
    return "\n".join(source)


def check_license():
    """
    Check the license for all source files.
    """

    license_text = get_license_text()
    source_files = get_all_source_files()

    skip_list = []
    wrong_headers = []

    for key in source_files:
        header = license_to_source(license_text, key)
        for file in source_files[key]:
            for skip in skip_list:
                if file.endswith(skip):
                    break
            else:
                with open(file) as source_file:
                    source_text = source_file.read()
                    if not source_text.startswith(header):
                        wrong_headers.append(file)

    return wrong_headers


class TestHeaders(unittest.TestCase):
    """This class tests the headers in the repository."""

    def test_headers(self):
        """
        Check if all headers are correct.
        """

        wrong_headers = check_license()
        wrong_headers_string = "Wrong headers in: " + ", ".join(wrong_headers)
        self.assertTrue(len(wrong_headers) == 0, wrong_headers_string)


if __name__ == "__main__":
    # Execution part of script.
    unittest.main()
