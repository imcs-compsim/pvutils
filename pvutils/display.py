


import os

from paraview.simple import *


def load_file(path):
    
    _dummy, extension = os.path.splitext(path)
    extension = extension.split('.')[-1].lower()
    
    if extension == 'pvd':
        return PVDReader(FileName=path)



def display_data(data):
    
    pass
