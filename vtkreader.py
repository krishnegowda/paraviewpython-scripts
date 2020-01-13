#!/usr/bin/pvpython

import vtk
import numpy as np
import foamIO
import paraview.simple as paraview 
from vtk.numpy_interface import dataset_adapter as dsa
import vtkio 
import matplotlib.pyplot as plt

folder='/scratch/krivi/nozzle-simulations/0_05_0_1_million/'
reader=foamIO.OpenFOAMReader(folder+'0.foam')  

reader.readLast()
data=vtkio.getBlockByName(reader.getDataSet(), 'internalMesh') 
#Reads data set and name of block

# Generate a line
lineSource = vtk.vtkLineSource()
lineSource.SetPoint1([0, 0, 0])
lineSource.SetPoint2([0.00115, 0, 0])
lineSource.SetResolution(100)
lineSource.Update()

line = lineSource.GetOutput()

# Probe the dataset to the line
probeFilter = vtk.vtkProbeFilter()
probeFilter.SetSourceData(data)    # Dataset to be probed from
probeFilter.SetInputData(line)     # Dataset to be probed to
probeFilter.Update()

probedData = probeFilter.GetOutput()

# Wrap data in a dataset_adapter for ease of access
ds = dsa.WrapDataObject(probedData)

points = ds.Points
U = ds.PointData['U']
p = ds.PointData['p']

plt.plot(points[:, 0], U[:, 0])
plt.show()


X, Y = np.meshgrid(np.arange(0, 2 * np.pi, .2), np.arange(0, 2 * np.pi, .2))
U = np.cos(X)
V = np.sin(Y)

plt.figure()
plt.title('Arrows scale with plot width, not view')
Q = plt.quiver(X, Y, U, V, units='width')
qk = plt.quiverkey(Q, 0.9, 0.9, 2, r'$2 \frac{m}{s}$', labelpos='E',
                   coordinates='figure')

plt.figure()
plt.title("pivot='mid'; every third arrow; units='inches'")
Q = plt.quiver(X[::3, ::3], Y[::3, ::3], U[::3, ::3], V[::3, ::3],
               pivot='mid', units='inches')
qk = plt.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E',
                   coordinates='figure')
plt.scatter(X[::3, ::3], Y[::3, ::3], color='r', s=5)

plt.figure()
plt.title("pivot='tip'; scales with x view")
M = np.hypot(U, V)
Q = plt.quiver(X, Y, U, V, M, units='x', pivot='tip', width=0.022,
               scale=1 / 0.15)
qk = plt.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E',
                   coordinates='figure')
plt.scatter(X, Y, color='k', s=5)

plt.show()
