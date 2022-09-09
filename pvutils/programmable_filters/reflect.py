"""
Reflect a dataset in ParaView. Sometimes ParaView has troubles reflecting
polylines.

Use with programmable filter in paraview and execute with:
>> exec(open('path to file').read())
"""

# Import paraview modules.
import paraview
import numpy as np

# Get "keyword arguments".
origin = np.array([0.0, 0.0, 0.0])
normal_vector = np.array([1.0, 0.0, 0.0])
if hasattr(paraview, 'programmable_filter_kwargs'):
    pv_dir = paraview.programmable_filter_kwargs[kwargs_id]
    if 'origin' in pv_dir.keys():
        origin = np.array(pv_dir['origin'])
    if 'normal_vector' in pv_dir.keys():
        normal_vector = np.array(pv_dir['normal_vector'])

# Get input and output objects.
pdi = self.GetInput()
pdo = self.GetOutput()

# Add old fields.
for i in range(pdi.GetCellData().GetNumberOfArrays()):
    pdo.GetCellData().AddArray(pdi.GetCellData().GetArray(i))
for i in range(pdi.GetPointData().GetNumberOfArrays()):
    pdo.GetPointData().AddArray(pdi.GetPointData().GetArray(i))

# Reflect origin and normal
origin = np.array([0.0, 0.0, 0.0])
normal_vector = np.array([1.0, 0.0, 0.0])
normal_vector = np.array(normal_vector / np.linalg.norm(normal_vector))

# Get the reflection matrix A.
A = np.eye(3) - 2. * np.dot(
    np.transpose(np.asmatrix(normal_vector)),
    np.asmatrix(normal_vector)
    )

# Reflect the points.
n_points = pdo.GetNumberOfPoints()
for i in range(n_points):
    pos_old = np.array(pdo.GetPoints().GetPoint(i))
    pos_old -= origin
    pos_new = pos_old * A
    pos_new += origin
    pdo.GetPoints().SetPoint(i, (pos_new[0, 0], pos_new[0, 1], pos_new[0, 2]))
