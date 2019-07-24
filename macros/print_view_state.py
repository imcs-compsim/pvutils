# -*- coding: utf-8 -*-
"""
Functions to simplify the use of scripts in paraview.
"""


# Import paraview module.
import paraview
import paraview.simple as pa


def _print_attibutes(obj, attributes, variable_name):
    """
    Print attributes of an object, given by a list of strings.
    """

    for att in attributes:
        attribute = getattr(obj, att)
        if isinstance(attribute, str):
            att_str = '\'{}\''.format(attribute)
            att_str = att_str.replace('\\', '\\\\')
        elif isinstance(attribute, float):
            att_str = '{0:.6g}'.format(attribute)
        elif type(attribute) == paraview.servermanager.VectorProperty:
            att_str_list = []
            for val in attribute:
                att_str = att_str_list.append('{0:.6g}'.format(val))
            att_str = '[{}]'.format(', '.join(att_str_list))
        else:
            att_str = str(attribute)
        print('{}.{} = {}'.format(variable_name, att, att_str))


def print_view_state(view, *args, **kwargs):
    """
    Allow the user to setup and return the relevant values.
    """

    # Display the view attributes.
    attributes = [
        'CameraPosition',
        'CameraFocalPoint',
        'CameraViewUp',
        'CameraViewAngle',
        'CameraParallelScale',
        'OrientationAxesVisibility',
        'CameraParallelProjection',
        'ViewSize'
        ]
    _print_attibutes(view, attributes, 'view')
    print('')

    # If additional items are given to this function, print their properties.
    for arg in args:

        item_string = str(arg)
        if 'ScalarBarWidgetRepresentation' in item_string:
            # Display the color bar attributes.
            attributes = [
                'Title',
                'ComponentTitle',
                'WindowLocation',
                'Orientation',
                'ScalarBarLength',
                'ScalarBarThickness',
                'Position'
                ]
            print('')
            _print_attibutes(arg, attributes, 'color_bar')


view = view = pa.GetActiveViewOrCreate('RenderView')
print_view_state(view)
