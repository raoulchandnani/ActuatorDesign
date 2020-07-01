from abaqus import *
from abaqusConstants import *
import __main__
import section
import odbSection
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import connectorBehavior
import displayGroupOdbToolset as dgo
from math import atan, sin, cos, tan
from Post_P_Script import getResults,getResults2,find_strain,find_stress
import numpy as np
from math import atan2,cos,sin,pi
with open("Aero_Data.txt", "r") as a:
    Data= [[float(num) for num in line.split(' ')] for line in a]
Data=np.array(Data)    

## Slice morphing data ##
count=0
Datalow=[[0.000000 for x in range(7)] for y in range(len(Data))]
Datalow=np.array(Datalow)
for i in range(len(Data)):
    if Data[i,6]!=0:
        if abs(Data[i,1])<0.15:
            Datalow[count]=Data[i]
            count=count+1
Datalow=Datalow[~np.all(Datalow == 0, axis=1)]

s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
	type=DEFORMABLE_BODY)
d2 = p.datums
loftsections=[]
for xloc in range(10):
	if xloc>0:
		DP=p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0315*xloc)
		DPP=DP.id
		p = mdb.models['Model-1'].parts['Part-1']
		DP=p.DatumAxisByTwoPoint(point1=(Datalow[-1,1],Datalow[0,2]-0.0315*xloc,0.0315*xloc), point2=(Datalow[-1,1],Datalow[0,2]+0.0315*xloc,0.0315*xloc))
		DPE=DP.id
		t = p.MakeSketchTransform(sketchPlane=d2[DPP], sketchUpEdge=d2[DPE], 
			sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0, 
			0, 0.0315*xloc))
		s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
			sheetSize=198.85, gridSpacing=4.97, transform=t)
		g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
		s.setPrimaryObject(option=SUPERIMPOSE)
		p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
	XLOC=11.85+0.0315*xloc
	count=0
	Dataxloc=[[0.000000 for x in range(7)] for y in range(len(Datalow))]
	Dataxloc=np.array(Dataxloc)
	for i in range(len(Datalow)):
		if abs(Datalow[i,0]-XLOC)<1E-06:
			Dataxloc[count]=Datalow[i]
			count=count+1
	Dataxloc=Dataxloc[~np.all(Dataxloc == 0, axis=1)]	
	for i in range(len(Dataxloc)-2):
		s.Spline(points=((Dataxloc[i,1], Dataxloc[i,2]-Datalow[0,2]), (Dataxloc[i+1,1], Dataxloc[i+1,2]-Datalow[0,2])))
	s.Line(point1=(Dataxloc[0,1], Dataxloc[0,2]-Datalow[0,2]), point2=(Dataxloc[0,1], Dataxloc[0,2]-Datalow[0,2]+0.1))
	s.Line(point1=(Dataxloc[0,1], Dataxloc[0,2]-Datalow[0,2]+0.1), point2=(Dataxloc[-2,1], Dataxloc[-2,2]-Datalow[0,2]+0.1))
	s.Line(point1=(Dataxloc[-2,1], Dataxloc[-2,2]-Datalow[0,2]), point2=(Dataxloc[-2,1], Dataxloc[-2,2]-Datalow[0,2]+0.1))	
	p = mdb.models['Model-1'].parts['Part-1']
	if xloc==0:
		p.BaseWire(sketch=s)
	else:
		p.Wire(sketchPlane=d2[DPP], sketchUpEdge=d2[DPE], sketchPlaneSide=SIDE1, 
			sketchOrientation=RIGHT, sketch=s)
	s.unsetPrimaryObject()
	p = mdb.models['Model-1'].parts['Part-1']
	del mdb.models['Model-1'].sketches['__profile__']
	p = mdb.models['Model-1'].parts['Part-1']
	p = mdb.models['Model-1'].parts['Part-1']
	e = p.edges
	edges = e.getByBoundingBox(-1000000,-1000000,0.0315*xloc-0.0001,1000000,1000000,0.0315*xloc+0.0001)
	p.Set(edges=edges, name='E%i'%xloc)
	loftsections.append(p.sets['E%i'%xloc].edges[:])
loftsections=tuple(loftsections)
p.SolidLoft(loftsections=loftsections, startCondition=NONE, 
    endCondition=NONE)
stahp

