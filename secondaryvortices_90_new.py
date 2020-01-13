# !/usr/bin/pvpython
import sys
import os
sys.path.insert(0, os.path.expanduser('~/.local/lib/python2.7/site-packages'))
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import vtk
#from vtktools import foamIO
from vtktools import foamIO
from vtktools import vtkio
import paraview.simple as paraview 
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
#import vtkio 
from vtk.numpy_interface import algorithms as algs
from paraview.servermanager import *

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
#paraview.simple._DisableFirstRenderCameraReset()


#create a new 'OpenFOAMReader' reading a datafile
a0foam = OpenFOAMReader('/Users/krishna/test_python_installation/data/0.foam')
if a0foam:
    print("Success")
else:
   print("Failed")

a0foam.MeshRegions = ['internalMesh']
a0foam.CellArrays = ['U', 'alpha.water', 'nu', 'p', 'p_rgh']

#view.ViewTime
#for finding last timestep in a series of data 
a0foam.TimestepValues
tsteps = a0foam.TimestepValues
a0foam.UpdatePipeline(max(tsteps))
#RenderView1 = GetRenderView()
#RenderView1.ViewTime = max(tsteps)
#  or
animationScene1 = GetAnimationScene()
animationScene1.GoToLast()
#animationScene1.UpdateAnimationUsingDataTimeSteps=tsteps[-1]

#Take slices
x_list = np.arange(0.006001,0.016,0.0005)
FileName='/scratch/krivi/geometry_simulations/geomet_post/Velocity_gradients_Orientation/150-deg/'
#value1=1
#value2=1
#value3=0

#figsize=(10,10)

#fig = plt.figure()
#fig,axlist = plt.subplots(4,5,sharex=True, sharey=True,figsize=figsize )
#print(axlist.flat)
#fig.subplots_adjust(hspace=5,wspace=5)
min_el_array=np.zeros(shape=20)           
iter =0
max_el_array=np.zeros(shape=20)
for i in range(len(x_list)): #taking slices
      slice1 = Slice(Input=a0foam)
      slice1.SliceType = 'Plane'
      slice1.SliceOffsetValues = [0.0]
      slice1.SliceType.Origin = [x_list[i],0.0055,0.0005]
      slice1.SliceType.Normal = [1,0,0]
      slice1.Triangulatetheslice = 0
      print(iter)
      #Use vtk wrapper to convert to numpy array
      slice_data = vtkio.getBlockByName(servermanager.Fetch(slice1), 'internalMesh')
      sliceds = dsa.WrapDataObject(slice_data)
      points = np.copy(sliceds.Points)
      inds = np.argsort(points[:, 1] + 1e4*points[:, 2]) ##sorting using weighted functin
      points = points[inds, :]   
      Ny=np.where(points[:,1] == points[0,1])[0][1] ## Selecting the 1st index value of y-cordinate
      Nz=points.shape[0]//Ny  
      U = np.copy(sliceds.PointData['U'][inds, :])
      #U=np.array(U)
      alpha = np.copy(sliceds.PointData['alpha.water'])
      Ux=U[:,0]*1000
      Uy=U[:,1]*1000
      Uz=U[:,2]*1000
      #print(points[:,1])
      #print(points[:,2])
      Umag = np.sqrt(Ux * Ux + Uy * Uy+ Uz*Uz)
      #print(np.max(Umag))
      Uy=(U[:,1]*1000)/np.max(Umag)
      Uz=(U[:,2]*1000)/np.max(Umag)
      Umag=np.reshape(Umag,(Nz,Ny))
      Uy=np.reshape(U[:,1]*1000,(Nz,Ny))
      Uz=np.reshape(U[:,2]*1000,(Nz,Ny))
      Y=np.reshape(points[:,1]*1000-5.5,(Nz,Ny))
      Z=np.reshape(points[:,2]*1000-0.5,(Nz,Ny))
      #Y,Z= np.meshgrid(np.linspace(0,51,0.1), np.linspace(1,51,0.1))
      
      #print(value1)
      #ax=plt.subplot(3,6,value2)

      maxmag= np.max(Umag)
      minmag= np.min(Umag)


      max_el_array[iter] = maxmag
      min_el_array[iter] = minmag

      max_value = np.amax(max_el_array)
      
      min_value = np.amin(min_el_array)
      iter=iter+1
      

print("MaxValue: " + str(np.amax(max_el_array)))
print("MaxValue: " + str(np.amin(min_el_array)))

#x_list = np.arange(0.006001,0.016,0.0005)
#FileName='/scratch/krivi/geometry_simulations/geomet_post/Velocity_gradients_Orientation/150-deg/'
value1=1
value2=1
value3=0

figsize=(10,10)

#fig = plt.figure()
fig,axlist = plt.subplots(4,5,sharex=True, sharey=True,figsize=figsize )
#print(axlist.flat)
#fig.subplots_adjust(hspace=5,wspace=5)
#min_el_array=np.zeros(shape=20)           
#iter =0
#max_el_array=np.zeros(shape=20)
for ax, i in zip(axlist.flat,x_list.flat):
      slice1 = Slice(Input=a0foam)
      slice1.SliceType = 'Plane'
      slice1.SliceOffsetValues = [0.0]
      slice1.SliceType.Origin = [i,0.0055,0.0005]
      slice1.SliceType.Normal = [1,0,0]
      slice1.Triangulatetheslice = 0
      print(iter)
      #Use vtk wrapper to convert to numpy array
      slice_data = vtkio.getBlockByName(servermanager.Fetch(slice1), 'internalMesh')
      sliceds = dsa.WrapDataObject(slice_data)
      points = np.copy(sliceds.Points)
      inds = np.argsort(points[:, 1] + 1e4*points[:, 2]) ##sorting using weighted functin
      points = points[inds, :]   
      Ny=np.where(points[:,1] == points[0,1])[0][1] ## Selecting the 1st index value of y-cordinate
      Nz=points.shape[0]//Ny  
      U = np.copy(sliceds.PointData['U'][inds, :])
      #U=np.array(U)
      alpha = np.copy(sliceds.PointData['alpha.water'])
      Ux=U[:,0]*1000
      Uy=U[:,1]*1000
      Uz=U[:,2]*1000
      #print(points[:,1])
      #print(points[:,2])
      Umag = np.sqrt(Ux * Ux + Uy * Uy+ Uz*Uz)
      #print(np.max(Umag))
      Uy=(U[:,1]*1000)/np.max(Umag)
      Uz=(U[:,2]*1000)/np.max(Umag)
      Umag=np.reshape(Umag,(Nz,Ny))
      Uy=np.reshape(U[:,1]*1000,(Nz,Ny))
      Uz=np.reshape(U[:,2]*1000,(Nz,Ny))
      Y=np.reshape(points[:,1]*1000-5.5,(Nz,Ny))
      Z=np.reshape(points[:,2]*1000-0.5,(Nz,Ny))
      #Y,Z= np.meshgrid(np.linspace(0,51,0.1), np.linspace(1,51,0.1))
      
      #print(value1)
      #ax=plt.subplot(3,6,value2)

      #maxmag= np.max(Umag)
      #minmag= np.min(Umag)


      #max_el_array[iter] = maxmag
      #min_el_array[iter] = minmag

      #max_value = np.amax(max_el_array)
      
      #min_value = np.amin(min_el_array)
      #iter=iter+1
      #cont = ax.contourf(Y, Z, Umag, np.linspace(0,25,30), cmap='jet')
      
      cont = ax.contourf(Y, Z, Umag, np.linspace(min_value,max_value,20), cmap='jet')
      skip=3
      quiv = ax.quiver(Y[::skip,::skip], Z[::skip,::skip], Uy[::skip, ::skip], Uz[::skip, ::skip], width=0.003, headwidth=3, headlength=5,)
      #quiv = ax.quiver(Y[::3,::3], Z[::3,::3], Uy[::3, ::3], Uz[::3, ::3], width=0.001)
      #quiv = plt.quiver(Y, Z, Uy, Uz, width=0.0001, scale=10,headwidth=3., headlength=4., color='black')
      #ax= plt.gca()
      #qk = plt.quiverkey(quiv, 0.9, 09, 1,'',labelpos='E', coordinates='figure')
      #plt.colorbar()
      #plt.gca().set_aspect('equal')
      #fig.colorbar(cont, orientation='vertical', ticks=limit)
      ax.set_title('$x/h={}$'.format(value1))
      value1=value1+0.5
      #ax.set_xlim(-0.5,0.5)
      #ax.set_ylim(-0.5,0.5)
      ax.set_xticks(np.arange(-0.5,0.6,0.25))
      ax.set_yticks(np.arange(-0.5,0.6,0.25))
      #ax.set_rasterized(True)
      #ax.set_aspect('equal')
      
      #value2=value2+1
    
#for index_row, row in enumerate(axlist.flat):
#    for index_column, column in enumerate(row):
#        if index_row == 2:
#            ax.set_xlabel('$y/h$')
      #break

#print("MaxValue: " + str(np.amax(max_el_array)))
#print(enumerate(row))
plt.tight_layout()
limit=np.linspace(min_value,max_value,20) #np.arange(0,30,5)
cbaxes= fig.add_axes([0.99,0.3,0.02,0.4])
cb = plt.colorbar(cont,cax=cbaxes, orientation='vertical', ticks=limit)
#fig.colorbar(cont, orientation='vertical', ticks=limit, location='right',shrink=0.6)
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = "10"
fig.text(0.52,0.0,'$y/h$', ha='center', va='center')
fig.text(0.01,0.52,'$z/h$', ha='center', va='center', rotation='vertical')
fig.text(0.52,1.00,'$90$', ha='center', va='center')
fig.set_size_inches(11.69,8.27)
plt.savefig('velocity_new_90.eps', format='eps',bbox_inches='tight', bbox_padding=0.5, dpi=300)
plt.show(block=False) 
plt.close(fig)           
            #ax.set_xlabel("$y/h$")
        #if index_column == 0:
            #ax.set_ylabel('$z/h$')
        #ax.set_aspect('equal')            
            
