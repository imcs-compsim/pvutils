

import sys
sys.path
sys.path.append('/home/kubuntu/Desktop/pv')



import numpy as np

from pvutils.display import *
paraview.simple._DisableFirstRenderCameraReset()

data = load_file('/home/kubuntu/Desktop/xxx/xxx-structure-beams.pvd')
data = CellDatatoPointData(Input=data)
data = ExtractSurface(Input=data)

data = Tube(Input=data)
data.Scalars = [None, 'cross_section_radius']
data.VaryRadius = 'By Absolute Scalar'
data.NumberofSides = 20


data = load_file('/home/kubuntu/Desktop/xxx/xxx-structure.pvd')
data = WarpByVector(Input=data)
data.Vectors = ['POINTS', 'displacement']
data.ScaleFactor = 2.


def pp(a):
    for item in dir(a):
        print(item)
        
pp(data)
pp(data.GetFieldDataInformation())
print(data.GetPointDataInformation().keys())
print(data.GetPointDataInformation().values())
print(data.GetPointDataInformation()['displacement'])
pp(data.GetPointDataInformation()['displacement'])
print(data.GetPointDataInformation()['displacement'].GetRange(2))
print(data.GetPointDataInformation()['displacement'].FieldData.GetNumberOfArrays())
print('')
pp(data.PointData.GetArray('displacement'))
print(data.GetPointDataInformation()['displacement'].Proxy)
#print(data.GetFieldDataInformation().GetFieldData())

#view = GetActiveViewOrCreate('RenderView')
#display = Show(data, view)
#view.Update()
