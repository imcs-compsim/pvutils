# -*- coding: utf-8 -*-
"""
This script is used to test the functionality of pvutils.
"""

# Python imports.
import unittest
import os


# ParaView imports.
import pvutils.utils as pvutils


class TestPvutils(unittest.TestCase):
    """Test various stuff from the pvutils module."""

    def setUp(self):
        """
        This method is called before each test and resets the ParaView state.
        """

        # Set default values for global parameters.
        pvutils.reset_paraview()

    def test_pvpython(self):
        """
        Check if the is is_pvpython function works correctly.
        """

        # Get the environment variable IS_PVPYTHON. If it is not defined,
        # assume False.
        is_pvpython_environment_string = os.environ.get('IS_PVPYTHON', '0')
        if is_pvpython_environment_string == '0':
            is_pvpython_environment = False
        elif is_pvpython_environment_string == '1':
            is_pvpython_environment = True
        else:
            raise ValueError('Got wrong environment variable IS_PVPYTHON!')
        self.assertEqual(is_pvpython_environment, pvutils.is_pvpython())

    def test_load_document(self):
        """
        Check that the implemented mesh types can be loaded.
        """

        pass


if __name__ == '__main__':
    # Execution part of script.

    unittest.TextTestRunner().run(
        unittest.TestLoader().loadTestsFromTestCase(TestPvutils))
