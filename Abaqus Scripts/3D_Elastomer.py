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
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
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
with open("Aero_Data_I.txt", "r") as a:
    Datalow= [[float(num) for num in line.split(' ')] for line in a]
Datalow=np.array(Datalow)
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
p = mdb.models['Model-1'].Part(name='Elastomer', dimensionality=THREE_D, 
	type=DEFORMABLE_BODY)
d2 = p.datums
loftsectionsA=[]
for xloc in range(100):
	if xloc>0:
		DP=p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0315*xloc)
		DPP=DP.id
		p = mdb.models['Model-1'].parts['Elastomer']
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
	Dataxloc=[[0.000000 for x in range(3)] for y in range(len(Datalow))]
	Dataxloc=np.array(Dataxloc)
	for i in range(len(Datalow)):
		if abs(Datalow[i,0]-XLOC)<1E-06:
			Dataxloc[count]=Datalow[i]
			count=count+1
	Dataxloc=Dataxloc[~np.all(Dataxloc == 0, axis=1)]
	Pad_thk=-(min(Datalow[:,2])-min(Dataxloc[:,2]))-0.04#-(3.15-0.0315*xloc)*tan(0.06841917290255273)+0.04
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP8].pointOn
	Pad_thk=DVct2[2]+Datalow[0,2]-0.18+(3.15-0.0315*xloc)*tan(0.06841917290255273)
	for i in range(len(Dataxloc)-1):
#		s.Spline(points=((Dataxloc[i,1], -(Dataxloc[i,2]-min(Datalow[:,2])-DVct2[2])), (Dataxloc[i+1,1], -(Dataxloc[i+1,2]-min(Datalow[:,2])-DVct2[2]))))
		s.Spline(points=((Dataxloc[i,1], -(Dataxloc[i,2])), (Dataxloc[i+1,1], -(Dataxloc[i+1,2]))))
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP8].pointOn
	s.Line(point1=(DVct[0]-Pad, DVct2[2]-Pad_thk), point2=(Dataxloc[0,1], -(Dataxloc[0,2])))
#	s.Line(point1=(DVct[0]-Pad, DVct2[2]-Pad_thk), point2=(DVct[0]-Pad, -(Dataxloc[0,2])))
	s.Line(point1=(DVct[0]-Pad, DVct2[2]-Pad_thk), point2=(DVct2[0], DVct2[2]-Pad_thk))
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP5].pointOn
	s.Line(point1=(DVct2[0], DVct2[2]-Pad_thk), point2=(DVct[0], DVct[2]-Pad_thk))
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP8].pointOn
	s.Line(point1=(DVct[0], DVct[2]-Pad_thk), point2=(DVct2[0], DVct2[2]-Pad_thk))
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP5].pointOn
	s.Line(point1=(DVct2[0], DVct2[2]-Pad_thk), point2=(DVct[0], DVct[2]-Pad_thk))
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP8].pointOn
	s.Line(point1=(DVct[0], DVct[2]-Pad_thk), point2=(DVct2[0], DVct2[2]-Pad_thk))
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP5].pointOn
	s.Line(point1=(DVct2[0], DVct2[2]-Pad_thk), point2=(DVct[0], DVct[2]-Pad_thk))
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP8].pointOn
	s.Line(point1=(DVct[0], DVct[2]-Pad_thk), point2=(DVct2[0], DVct2[2]-Pad_thk))
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP5].pointOn
	s.Line(point1=(DVct2[0], DVct2[2]-Pad_thk), point2=(DVct[0], DVct[2]-Pad_thk))
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP8].pointOn
	s.Line(point1=(DVct[0], DVct[2]-Pad_thk), point2=(DVct2[0], DVct2[2]-Pad_thk))
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP5].pointOn
	s.Line(point1=(DVct2[0], DVct2[2]-Pad_thk), point2=(DVct[0], DVct[2]-Pad_thk))
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP9].pointOn
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP8].pointOn
	DVct3=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
	s.Line(point1=(DVct[0], DVct[2]-Pad_thk), point2=(-DVct3[0]+Pad, DVct2[2]-Pad_thk))
#	s.Line(point1=(Dataxloc[-1,1], -(Dataxloc[-1,2])), point2=(-DVct3[0]+Pad, -(Dataxloc[-2,2])))
	s.Line(point1=(Dataxloc[-1,1], -(Dataxloc[-1,2])), point2=(-DVct3[0]+Pad, DVct2[2]-Pad_thk))

	p = mdb.models['Model-1'].parts['Elastomer']
	if xloc==0:
		p.BaseWire(sketch=s)
	else:
		p.Wire(sketchPlane=d2[DPP], sketchUpEdge=d2[DPE], sketchPlaneSide=SIDE1, 
			sketchOrientation=RIGHT, sketch=s)
	s.unsetPrimaryObject()
	p = mdb.models['Model-1'].parts['Elastomer']
	del mdb.models['Model-1'].sketches['__profile__']
	p = mdb.models['Model-1'].parts['Elastomer']
	p = mdb.models['Model-1'].parts['Elastomer']
	e = p.edges
	edges = e.getByBoundingBox(-1000000,-1000000,0.0315*xloc-0.0001,1000000,1000000,0.0315*xloc+0.0001)
	p.Set(edges=edges, name='E%i'%xloc)
xlist=[]
for xloc in range(100):
	xlist.append(xloc)
p.SolidLoft(loftsections=[p.sets['E%i'%xloc].edges[:] for xloc in xlist], startCondition=NORMAL, endCondition=NORMAL)
# p.SolidLoft(loftsections=[p.sets['E%i'%xloc].edges[:] for xloc in range(43,71)], startCondition=SPECIFIED, 
    # startTangent=90.0, startMagnitude=5.0, endCondition=SPECIFIED, 
    # endTangent=90.0, endMagnitude=5.0)
