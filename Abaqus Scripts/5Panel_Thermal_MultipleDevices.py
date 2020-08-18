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
from math import atan, sin, cos, tan,atan2,pi
import numpy as np
from Post_P_Script import getResults,getResults2,find_strain,find_stress
##########################

# Variables Assembly ##
ND = 4 # Number of devices

D=0.015/4 # Hinge distance from edge
Total = 0.28 # Total spanwise length
Pad=0.01 # Excess padding in spanwise dierection
LD=0.0315 # Distance between devices in keelwise direction

## Angles ##
phiL2A=13.25 
phiR2A=13.25 
phiLA=5.79
phiRA=5.79
phiMA=0;
phiL2B=12.07; 
phiR2B=12.07;
phiLB=3.07;
phiRB=3.07;
phiMB=0;
phiL2C=7.78; 
phiR2C=7.78;
phiLC=9.16;
phiRC=9.16;
phiMC=0;

phiL2D=12.07; 
phiR2D=12.07; 
phiLD=3.07;
phiRD=3.07;
phiMD=0;


## Variables for Conformal Surface ##
start=30 # start location for conformal surface from aero data
end=98 # end location for conformal surface from aero data
min_thk=0.01 # minimum thickness of conformal surface

## Load conformal surface aero data ##
with open("Aero_Data_I.txt", "r") as a:
    Datalow= [[float(num) for num in line.split(' ')] for line in a]
Datalow=np.array(Datalow)
 
## Variables for TT ##
T1A=1.2
T3A=1.2
T1B=1.32
T3B=1.32
T1C=1.078
T3C=1.078

T1D=1.32
T3D=1.32

ORTT=0.014
IRTT=0.007
TTL=0.1275

## Variables for Plates ##
VA=0.663; # Vertical plate lengths
VB=0.537*0.5
VC=0.631
VD=0.537*0.5

GA=9.2/1000
HNGA=-8.3/1000
H1A=0.021
H2A=0.065
H3A=Total-2*(H1A+GA)*cos(phiL2A*pi/180)-2*(H2A+GA)*cos(phiLA*pi/180)-GA
H4A=H2A
H5A=H1A

GB=7.7/1000
HNGB=-6.0/1000
H1B=0.039
H2B=0.060
H3B=Total-2*(H1B+GB)*cos(phiL2B*pi/180)-2*(H2B+GB)*cos(phiLB*pi/180)-GB
H4B=H2B
H5B=H1B

GC=0.01
HNGC=-4.4/1000
H1C=0.014
H2C=0.074
H3C=Total-2*(H1C+GC)*cos(phiL2C*pi/180)-2*(H2C+GC)*cos(phiLC*pi/180)-GC
H4C=H2C
H5C=H1C

GD=7.7/1000
HNGD=-6.0/1000
H1D=0.039
H2D=0.060
H3D=Total-2*(H1D+GD)*cos(phiL2D*pi/180)-2*(H2D+GD)*cos(phiLD*pi/180)-GD
H4D=H2D
H5D=H1D

## Stops code before submitting job if test=1 ##
test=1

rownamelist=[]
for i in range(ND):
	rownamelist.append(chr(65+i))

## Adjusts plate widths to match aero data slices ##
VA=ceil(VA/0.0315738727171977)*0.0315738727171977
VB=ceil(VB/0.0315738727171977)*0.0315738727171977
VC=ceil(VC/0.0315738727171977)*0.0315738727171977
VB=ceil(VD/0.0315738727171977)*0.0315738727171977
LD=ceil(LD/0.0315738727171977)*0.0315738727171977

##########################
### Write data file column headings
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

##########################
###Sketch Geometry and Create Parts###
print 'Sketching/Creating the part'
Mdb()

## Sketch geometry of TT ##

s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=30.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.Spot(point=(0.0, 0.0))
s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(ORTT, 0.0))
s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(IRTT, 0.0))
p = mdb.models['Model-1'].Part(name='Torque_Tube', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Torque_Tube']
p.BaseSolidExtrude(sketch=s, depth=0.5*(TTL))
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Torque_Tube']
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Torque_Tube']
f = p.faces
p.Mirror(mirrorPlane=f[3], keepOriginal=ON)

## Partition the geometry##
p = mdb.models['Model-1'].parts['Torque_Tube']
c = p.cells
pickedCells = c.findAt(((0.5*(IRTT+ORTT), 0, 0), ))
v1, e, d1 = p.vertices, p.edges, p.datums
p.PartitionCellByPlaneThreePoints(point2=v1.findAt(coordinates=(ORTT, 0.0, 
    (0.5*TTL))), point3=v1.findAt(coordinates=(ORTT, 0.0, -(0.5*TTL))), cells=pickedCells, 
    point1=p.InterestingPoint(edge=e.findAt(coordinates=(0,-ORTT, (0.5*TTL))), 
    rule=MIDDLE))
p = mdb.models['Model-1'].parts['Torque_Tube']
c = p.cells
pickedCells = c.findAt(((0, 0.5*(IRTT+ORTT), 0), ), ((0, 
     -0.5*(IRTT+ORTT), 0), ))
v2, e1, d2 = p.vertices, p.edges, p.datums
p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=p.InterestingPoint(
    edge=e1.findAt(coordinates=(0, ORTT, (0.5*TTL))), rule=MIDDLE), 
    point2=p.InterestingPoint(edge=e1.findAt(coordinates=(0, -ORTT, 
    (0.5*TTL))), rule=MIDDLE), point3=p.InterestingPoint(edge=e1.findAt(
    coordinates=(0, ORTT, -(0.5*TTL))), rule=MIDDLE))
p = mdb.models['Model-1'].parts['Torque_Tube']
DP=p.DatumPointByCoordinate(coords=(0, 0, (0.5*TTL)))
DP_TT1=DP.id
p = mdb.models['Model-1'].parts['Torque_Tube']
DP=p.DatumPointByCoordinate(coords=(0, 0, -(0.5*TTL)))
DP_TT2=DP.id

## Create sets##
p = mdb.models['Model-1'].parts['Torque_Tube']
e = p.edges
edges = e.findAt(((0.0, ORTT, 0), ), ((0.0, -ORTT, 0), ), ((-ORTT, 0.0, 
    0), ), ((ORTT, 0.0, 0), ))
p.Set(edges=edges, name='Longitudinal_edges')

p = mdb.models['Model-1'].parts['Torque_Tube']
e = p.edges
edges = e.findAt(((0.0, 0.5*(IRTT+ORTT), -(0.5*TTL)), ), ((0.0, 0.5*(IRTT+ORTT), (0.5*TTL)), ), ((0.0, 
    -0.5*(IRTT+ORTT), (0.5*TTL)), ), ((0.0, -0.5*(IRTT+ORTT), -(0.5*TTL)), ), ((0.5*(IRTT+ORTT), 0.0, (0.5*TTL)), ), ((
    -0.5*(IRTT+ORTT), 0.0, -(0.5*TTL)), ), ((0.5*(IRTT+ORTT), 0.0, -(0.5*TTL)), ), ((-0.5*(IRTT+ORTT), 0.0, (0.5*TTL)), ))
p.Set(edges=edges, name='Radial_edges')

p = mdb.models['Model-1'].parts['Torque_Tube']
e = p.edges
edges = e.findAt(((ORTT*cos(pi/4), ORTT*sin(pi/4), -(0.5*TTL)), ), ((-ORTT*cos(pi/4), ORTT*sin(pi/4), 
    (0.5*TTL)), ), ((ORTT*cos(pi/4), ORTT*sin(pi/4), (0.5*TTL)), ), ((-ORTT*cos(pi/4), -ORTT*sin(pi/4), -(0.5*TTL)), 
    ), ((ORTT*cos(pi/4), -ORTT*sin(pi/4), (0.5*TTL)), ), ((-ORTT*cos(pi/4), -ORTT*sin(pi/4), (0.5*TTL)), ), (
    (-ORTT*cos(pi/4), ORTT*sin(pi/4), -(0.5*TTL)), ), ((ORTT*cos(pi/4), -ORTT*sin(pi/4), -(0.5*TTL)), ))
p.Set(edges=edges, name='Circumferential_edges')

## Sketch geometry of Plates ##
for rowname in rownamelist:
	for plate_num in [1,2,3,4,5]:
		s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
			sheetSize=200.0)
		g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
		s.setPrimaryObject(option=STANDALONE)
		s.rectangle(point1=((-0.5*eval('H%i%s'%(plate_num,rowname))), (-0.5*eval('V%s'%rowname))), point2=(0.5*eval('H%i%s'%(plate_num,rowname)),0.5*eval('V%s'%rowname)))
		p = mdb.models['Model-1'].Part(name='Plate%i%s'%(plate_num,rowname), dimensionality=THREE_D, 
			type=DEFORMABLE_BODY)
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		p.BaseShell(sketch=s)
		s.unsetPrimaryObject()
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		del mdb.models['Model-1'].sketches['__profile__']

		#Defining the face partitions
		print 'Partitioning part'
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		f1, e2, d2 = p.faces, p.edges, p.datums
		t = p.MakeSketchTransform(sketchPlane=f1.findAt(coordinates=(0, 
			0, 0.0), normal=(0.0, 0.0, 1.0)), sketchUpEdge=e2.findAt(
			coordinates=(0.5*eval('H%i%s'%(plate_num,rowname)),0, 0.0)), sketchPlaneSide=SIDE1, origin=(0, 0,0.0))
		s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
			sheetSize=459.88, gridSpacing=11.49, transform=t)
		g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
		s.setPrimaryObject(option=SUPERIMPOSE)
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

		s.rectangle(point1=(-(0.5*eval('H%i%s'%(plate_num,rowname))-D), -(0.5*eval('V%s'%rowname)-D)), point2=((0.5*eval('H%i%s'%(plate_num,rowname))-D), (0.5*eval('V%s'%rowname)-D)))
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		f = p.faces
		pickedFaces = f.findAt(((0,0, 0.0), ))
		e, d1 = p.edges, p.datums
		p.PartitionFaceBySketch(sketchUpEdge=e.findAt(coordinates=(0.5*eval('H%i%s'%(plate_num,rowname)),0, 0.0)), 
			faces=pickedFaces, sketch=s)
		s.unsetPrimaryObject()
		del mdb.models['Model-1'].sketches['__profile__']
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		s = p.faces
		side1Faces = s.findAt(((0, 0, 0.0), ))
		# side1Faces = s.findAt(((0, 0, 0.0), ), (((0.5*H-0.5*D), 0.25*V, 0.0), ), (((0.5*H-0.5*D), -0.25*V, 0.0), ), ((
			# -(0.5*H-0.5*D), 0.25*V, 0.0), ), ((-(0.5*H-0.5*D), -0.25*V, 0.0), ))
		p.Surface(side2Faces=side1Faces, name='Surf-Bottom')
		###Datum Points
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		DP=p.DatumPointByCoordinate(coords=((0.5*eval('H%i%s'%(plate_num,rowname))-D), 0.5*eval('V%s'%rowname)-D, 0))
		DP2=DP.id
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		DP=p.DatumPointByCoordinate(coords=((0.5*eval('H%i%s'%(plate_num,rowname))-D), -0.5*eval('V%s'%rowname)+D, 0))
		DP1=DP.id
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		DP=p.DatumPointByCoordinate(coords=(-(0.5*eval('H%i%s'%(plate_num,rowname))-D), 0.5*eval('V%s'%rowname)-D, 0))
		DP3=DP.id
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		DP=p.DatumPointByCoordinate(coords=(-(0.5*eval('H%i%s'%(plate_num,rowname))-D), -0.5*eval('V%s'%rowname)+D, 0))
		DP4=DP.id
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		DP=p.DatumPointByCoordinate(coords=((0.5*eval('H%i%s'%(plate_num,rowname))), 0.5*eval('V%s'%rowname), 0))
		DP6=DP.id
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		DP=p.DatumPointByCoordinate(coords=((0.5*eval('H%i%s'%(plate_num,rowname))), -0.5*eval('V%s'%rowname), 0))
		DP5=DP.id
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		DP=p.DatumPointByCoordinate(coords=(-(0.5*eval('H%i%s'%(plate_num,rowname))), 0.5*eval('V%s'%rowname), 0))
		DP7=DP.id
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		DP=p.DatumPointByCoordinate(coords=(-(0.5*eval('H%i%s'%(plate_num,rowname))), -0.5*eval('V%s'%rowname), 0))
		DP8=DP.id
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		DP=p.DatumPointByCoordinate(coords=((0.5*(eval('H%i%s'%(plate_num,rowname))+eval('G%s'%rowname))), 0.5*eval('V%s'%rowname)-D, eval('HNG%s'%rowname)))
		DP10=DP.id
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		DP=p.DatumPointByCoordinate(coords=((0.5*(eval('H%i%s'%(plate_num,rowname))+eval('G%s'%rowname))), -0.5*eval('V%s'%rowname)+D, eval('HNG%s'%rowname)))
		DP9=DP.id
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		DP=p.DatumPointByCoordinate(coords=(-(0.5*(eval('H%i%s'%(plate_num,rowname))+eval('G%s'%rowname))), 0.5*eval('V%s'%rowname)-D, eval('HNG%s'%rowname)))
		DP11=DP.id
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		DP=p.DatumPointByCoordinate(coords=(-(0.5*(eval('H%i%s'%(plate_num,rowname))+eval('G%s'%rowname))), -0.5*eval('V%s'%rowname)+D, eval('HNG%s'%rowname)))
		DP12=DP.id
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		DP=p.DatumPointByCoordinate(coords=(0, 0, 0))
		DP13=DP.id

#Assemble Parts
print 'Placing Parts in Space'
#Create Instances here

for rowname in rownamelist:
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByDefault(CARTESIAN)
	p = mdb.models['Model-1'].parts['Plate1%s'%rowname]
	a.Instance(name='Left_Plate_2%s'%rowname, part=p, dependent=ON)
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByDefault(CARTESIAN)
	p = mdb.models['Model-1'].parts['Plate2%s'%rowname]
	a.Instance(name='Left_Plate%s'%rowname, part=p, dependent=ON)
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByDefault(CARTESIAN)
	p = mdb.models['Model-1'].parts['Plate3%s'%rowname]
	a.Instance(name='Middle_Plate%s'%rowname, part=p, dependent=ON)
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByDefault(CARTESIAN)
	p = mdb.models['Model-1'].parts['Plate4%s'%rowname]
	a.Instance(name='Right_Plate%s'%rowname, part=p, dependent=ON)
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByDefault(CARTESIAN)
	p = mdb.models['Model-1'].parts['Plate5%s'%rowname]
	a.Instance(name='Right_Plate_2%s'%rowname, part=p, dependent=ON)
	Vct=[0,0,0]	
	DVct=[-0.5*Total,-0.5*eval('V%s'%rowname)+D,0]
	for loc in range(3):
		Vct[loc]= DVct[loc] - mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP12].pointOn[loc]
	a = mdb.models['Model-1'].rootAssembly
	a.translate(instanceList=('Left_Plate_2%s'%rowname, ), vector=(Vct[0],Vct[1], Vct[2]))
	a = mdb.models['Model-1'].rootAssembly
	a.rotate(instanceList=('Left_Plate_2%s'%rowname, ), axisPoint=(DVct[0],DVct[1], DVct[2]), 
		axisDirection=(0.0, -(2*eval('V%s'%rowname)), 0.0), angle=(eval('phiL2%s'%rowname)))	
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP9].pointOn

	for loc in range(3):
		Vct[loc]= DVct[loc] - mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP12].pointOn[loc]
	a = mdb.models['Model-1'].rootAssembly
	a.translate(instanceList=('Left_Plate%s'%rowname, ), vector=(Vct[0],Vct[1], Vct[2]))
	a = mdb.models['Model-1'].rootAssembly
	a.rotate(instanceList=('Left_Plate%s'%rowname, ), axisPoint=(DVct[0],DVct[1], DVct[2]), 
		axisDirection=(0.0, -(2*eval('V%s'%rowname)), 0.0), angle=(eval('phiL%s'%rowname)))		
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP9].pointOn
	for loc in range(3):
		Vct[loc]= DVct[loc] - mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP12].pointOn[loc]
	a = mdb.models['Model-1'].rootAssembly
	a.translate(instanceList=('Middle_Plate%s'%rowname, ), vector=(Vct[0],Vct[1], Vct[2]))
	a = mdb.models['Model-1'].rootAssembly
	a.rotate(instanceList=('Middle_Plate%s'%rowname, ), axisPoint=(DVct[0],DVct[1], DVct[2]), 
		axisDirection=(0.0, -(2*eval('V%s'%rowname)), 0.0), angle=(eval('phiM%s'%rowname)))			

	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP9].pointOn
	for loc in range(3):
		Vct[loc]= DVct[loc] - mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP12].pointOn[loc]
	a = mdb.models['Model-1'].rootAssembly
	a.translate(instanceList=('Right_Plate%s'%rowname, ), vector=(Vct[0],Vct[1], Vct[2]))
	a = mdb.models['Model-1'].rootAssembly
	a.rotate(instanceList=('Right_Plate%s'%rowname, ), axisPoint=(DVct[0],DVct[1], DVct[2]), 
		axisDirection=(0.0, -(2*eval('V%s'%rowname)), 0.0), angle=(-eval('phiR%s'%rowname)))	

	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP9].pointOn
	for loc in range(3):
		Vct[loc]= DVct[loc] - mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP12].pointOn[loc]
	a = mdb.models['Model-1'].rootAssembly
	a.translate(instanceList=('Right_Plate_2%s'%rowname, ), vector=(Vct[0],Vct[1], Vct[2]))
	a = mdb.models['Model-1'].rootAssembly
	a.rotate(instanceList=('Right_Plate_2%s'%rowname, ), axisPoint=(DVct[0],DVct[1], DVct[2]), 
		axisDirection=(0.0, (2*eval('V%s'%rowname)), 0.0), angle=(eval('phiR2%s'%rowname)))

## Create rows ##
	if rowname=='A':
		Dist=0
	else:
		Dist=Dist+0.5*eval('V%s'%rowname)+LD
	a = mdb.models['Model-1'].rootAssembly
	a.translate(instanceList=('Middle_Plate%s'%rowname,'Right_Plate%s'%rowname,'Right_Plate_2%s'%rowname,'Left_Plate%s'%rowname,'Left_Plate_2%s'%rowname, ), vector=(0,Dist, 0))
	Dist=Dist+0.5*eval('V%s'%rowname)
	
	p = mdb.models['Model-1'].parts['Torque_Tube']
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByDefault(CARTESIAN)
	a.Instance(name='Torque_Tube_1%s'%rowname, part=p, dependent=ON)
	
#Create Reference Points
print 'Defining Reference Points'

## Reference points##
for rowname in rownamelist:
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP10].pointOn
	a = mdb.models['Model-1'].rootAssembly													# RPs for Middle to Right 1
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#8,9
	exec('RP1%s =RP.id'%rowname)
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP2%s =RP.id'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP9].pointOn
	a = mdb.models['Model-1'].rootAssembly													# RPs for Middle to Right 2
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))													#10,11
	exec('RP3%s =RP.id'%rowname)															
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP4%s =RP.id'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP11].pointOn
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP6%s =RP.id'%rowname)
	a = mdb.models['Model-1'].rootAssembly													# RPs for Middle to Left 1
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#12,13
	exec('RP5%s =RP.id'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP12].pointOn
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP7%s =RP.id'%rowname)
	a = mdb.models['Model-1'].rootAssembly														# RPs for Middle to Left 2
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#14,15
	exec('RP8%s =RP.id'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP10].pointOn
	a = mdb.models['Model-1'].rootAssembly															# RPs for Right to Right
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#16,17
	exec('RP9%s =RP.id'%rowname)
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP10%s =RP.id'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP9].pointOn
	a = mdb.models['Model-1'].rootAssembly															# RPs for Right to Right 2
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#18,19
	exec('RP11%s =RP.id'%rowname)
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP12%s =RP.id'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP11].pointOn
	a = mdb.models['Model-1'].rootAssembly														# RPs for Left to Left 
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#RP13,21
	exec('RP13%s =RP.id'%rowname)
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP14%s =RP.id'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP12].pointOn
	a = mdb.models['Model-1'].rootAssembly														# RPs for Left to Left  2
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#22,23
	exec('RP15%s =RP.id'%rowname)
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP16%s =RP.id'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP13].pointOn
	a = mdb.models['Model-1'].rootAssembly														# RPs for Rigid Body Constraints								#24,25,26
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP17%s =RP.id'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP13].pointOn
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP18%s =RP.id'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP13].pointOn
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP19%s =RP.id'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP10].pointOn
	a = mdb.models['Model-1'].rootAssembly																	# RPs for Hinge Support at endswith								#27,28,29,30
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP20%s =RP.id'%rowname)
	a = mdb.models['Model-1'].rootAssembly																	# RPs for Hinge Support at endswith								#27,28,29,30
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP26%s =RP.id'%rowname)
	r1 = a.referencePoints
	refPoints1=(r1[eval('RP20%s'%rowname)], )
	a.Set(referencePoints=refPoints1, name='RP_Z3+%s'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP9].pointOn
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1]+(eval('V%s'%rowname)-2*D-TTL), DVct[2]))
	exec('RP21%s =RP.id'%rowname)
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP27%s =RP.id'%rowname)
	r1 = a.referencePoints
	refPoints1=(r1[eval('RP21%s'%rowname)], )
	a.Set(referencePoints=refPoints1, name='RP_Z3-%s'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP11].pointOn
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1]-(eval('V%s'%rowname)-2*D-TTL), DVct[2]))
	exec('RP22%s =RP.id'%rowname)
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP28%s =RP.id'%rowname)
	r1 = a.referencePoints
	refPoints1=(r1[eval('RP22%s'%rowname)], )
	a.Set(referencePoints=refPoints1, name='RP_Z1-%s'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP12].pointOn
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP23%s =RP.id'%rowname)
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP29%s =RP.id'%rowname)
	r1 = a.referencePoints
	refPoints1=(r1[eval('RP23%s'%rowname)], )
	a.Set(referencePoints=refPoints1, name='RP_Z1+%s'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP13].pointOn
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP24%s =RP.id'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP13].pointOn
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	exec('RP25%s =RP.id'%rowname)

## Create Elastomer
Presloc=[0,0,0]
LMD=0
for rowname in rownamelist:
	LMD=LMD+eval('V%s'%rowname)
LMD=LMD+LD*(ND-1)
LEL=(0.0315*(end-1)-0.0315*start)/cos(6.841917290255273E-02)
for rowname in rownamelist:
	exec('NS%s=V%s*cos(6.841917290255273E-02)/0.0315'%(rowname,rowname))
LDS=LD*cos(6.841917290255273E-02)/0.0315

for rowname in rownamelist:
	if rowname=='A':
		startA=start+int(round(0.5*(LEL-LMD)*cos(6.841917290255273E-02)/0.0315,0))
	else:
		exec('start%s=Dist'%(rowname))
	exec('end%s=start%s+1+int(round(V%s*cos(6.841917290255273E-02)/0.0315,0))'%(rowname,rowname,rowname))
	Dist=eval('end%s'%rowname)-1+int(round(LD*cos(6.841917290255273E-02)/0.0315,0))

wire_count=0
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
	sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
p = mdb.models['Model-1'].Part(name='Elastomer', dimensionality=THREE_D, 
	type=DEFORMABLE_BODY)
d2 = p.datums
p.DatumCsysByThreePoints(name='Datum_Aero_Data', coordSysType=CARTESIAN, origin=(
    0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), point2=(1.0, 1.0, 0.0))

Curveloc=np.zeros((1,3))
Datamin=np.zeros(100)
for xloc in range(0,100):
	XLOC=11.85+0.0315*xloc
	count=0
	Dataxloc=[[0.000000 for x in range(7)] for y in range(len(Datalow))]
	Dataxloc=np.array(Dataxloc)
	for i in range(len(Datalow)):
		if abs(Datalow[i,0]-XLOC)<1E-06:
			Dataxloc[count]=Datalow[i]
			count=count+1
	Dataxloc=Dataxloc[~np.all(Dataxloc == 0, axis=1)]
	Datamin[xloc]=(min(-Dataxloc[:,2])-(0.0315*xloc)*tan(6.841917290255273E-02))
lowest=min(Datamin)

for rowname in rownamelist:
	Pathloc=(0,0,0,0,0,0)
	Setloc =np.zeros((5,3))	
	for xloc in range(eval('start%s'%rowname),eval('end%s'%rowname)):
		if xloc>0:
			DP=p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0315*xloc)
			DPP=DP.id
			exec('DPP%i=DP.id'%xloc)				
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
		Dataxloc=[[0.000000 for x in range(7)] for y in range(len(Datalow))]
		Dataxloc=np.array(Dataxloc)
		for i in range(len(Datalow)):
			if abs(Datalow[i,0]-XLOC)<1E-06:
				Dataxloc[count]=Datalow[i]
				count=count+1
		Dataxloc=Dataxloc[~np.all(Dataxloc == 0, axis=1)]
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2A'].datums[DP13].pointOn
		Pad_thk=DVct2[2]+Datalow[0,2]-0.18+(3.15-0.0315*xloc)*tan(6.841917290255273E-02)
		Pad_thk=(0.0315*xloc)*tan(6.841917290255273E-02)+lowest-DVct2[2]-min_thk
		j=0
		if xloc==(eval('start%s'%rowname)+eval('end%s'%rowname))/2:
			for i in range(len(Dataxloc)-1):
				Presloc=np.vstack((Presloc,[0.7*(Dataxloc[i,1])+0.3*(Dataxloc[i+1,1]),0.7*(-(Dataxloc[i,2]))-0.3*(Dataxloc[i+1,2]),0.0315*xloc]))
				Presloc=np.vstack((Presloc,[0.3*(Dataxloc[i,1])+0.7*(Dataxloc[i+1,1]),0.3*(-(Dataxloc[i,2]))-0.7*(Dataxloc[i+1,2]),0.0315*xloc]))
				s.Line(point1=(Dataxloc[i,1], -(Dataxloc[i,2])), point2=(Dataxloc[i+1,1], -(Dataxloc[i+1,2])))				
			Presloc=Presloc[~np.all(Presloc == 0, axis=1)]
		else:		
			for i in range(len(Dataxloc)-1):
				s.Spline(points=((Dataxloc[i,1], -(Dataxloc[i,2])), (Dataxloc[i+1,1], -(Dataxloc[i+1,2]))))
		for i in range(len(Dataxloc)):
			Curveloc=np.vstack((Curveloc,[Dataxloc[i,1],-(Dataxloc[i,2]),0.0315*xloc]))
		Curveloc=Curveloc[~np.all(Curveloc == 0, axis=1)]		
		DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP12].pointOn
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP8].pointOn
		s.Line(point1=(DVct[0]-Pad, DVct2[2]+Pad_thk), point2=(Dataxloc[0,1], -(Dataxloc[0,2])))
		s.Line(point1=(DVct[0]-Pad, DVct2[2]+Pad_thk), point2=(DVct2[0], DVct2[2]+Pad_thk))
		if xloc==eval('start%s'%rowname) or xloc==eval('end%s'%rowname)-1:
			Path=(11.85+xloc*0.0315,DVct[0]-Pad, DVct2[2]+Pad_thk)
		DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP5].pointOn
		s.Line(point1=(DVct2[0], DVct2[2]+Pad_thk), point2=(DVct[0], DVct[2]+Pad_thk))
		PlatePath=np.vstack(([0.0315*xloc,DVct2[0], DVct2[2]+Pad_thk],[0.0315*xloc,DVct[0], DVct[2]+Pad_thk]))	
		if xloc==(eval('start%s'%rowname)+eval('end%s'%rowname))/2:
			Setloc[0]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])+Pad_thk)
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP8].pointOn
		s.Line(point1=(DVct[0], DVct[2]+Pad_thk), point2=(DVct2[0], DVct2[2]+Pad_thk))
		DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP5].pointOn
		s.Line(point1=(DVct2[0], DVct2[2]+Pad_thk), point2=(DVct[0], DVct[2]+Pad_thk))
		PlatePath=np.vstack((PlatePath,[0.0315*xloc,DVct2[0], DVct2[2]+Pad_thk],[0.0315*xloc,DVct[0], DVct[2]+Pad_thk]))	
		if xloc==(eval('start%s'%rowname)+eval('end%s'%rowname))/2:
			Setloc[1]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])+Pad_thk)
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP8].pointOn
		s.Line(point1=(DVct[0], DVct[2]+Pad_thk), point2=(DVct2[0], DVct2[2]+Pad_thk))
		DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP5].pointOn
		s.Line(point1=(DVct2[0], DVct2[2]+Pad_thk), point2=(DVct[0], DVct[2]+Pad_thk))
		PlatePath=np.vstack((PlatePath,[0.0315*xloc,DVct2[0], DVct2[2]+Pad_thk],[0.0315*xloc,DVct[0], DVct[2]+Pad_thk]))					
		if xloc==(eval('start%s'%rowname)+eval('end%s'%rowname))/2:
			Setloc[2]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])+Pad_thk)			
		if xloc==eval('start%s'%rowname):
			Axis_start=(DVct2[0], DVct2[2]+Pad_thk,0.0315*xloc)
			Axis_end=(DVct[0], DVct[2]+Pad_thk,0.0315*xloc)
			Axis_1_avg=(0.5*(DVct[0]+DVct2[0]), 0.5*(DVct[2]+DVct2[2])+Pad_thk,0.0315*xloc)
			if rowname=='A':
				DP=p.DatumPointByCoordinate(coords=(0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])+Pad_thk , 0.0315*xloc))
				DPloc=DP.id
		if xloc==eval('end%s'%rowname)-1:
			Axis_2_avg=(0.5*(DVct[0]+DVct2[0]), 0.5*(DVct[2]+DVct2[2])+Pad_thk,0.0315*xloc)
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP8].pointOn
		s.Line(point1=(DVct[0], DVct[2]+Pad_thk), point2=(DVct2[0], DVct2[2]+Pad_thk))
		DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP5].pointOn
		s.Line(point1=(DVct2[0], DVct2[2]+Pad_thk), point2=(DVct[0], DVct[2]+Pad_thk))
		PlatePath=np.vstack((PlatePath,[0.0315*xloc,DVct2[0], DVct2[2]+Pad_thk],[0.0315*xloc,DVct[0], DVct[2]+Pad_thk]))		
		if xloc==(eval('start%s'%rowname)+eval('end%s'%rowname))/2:
			Setloc[3]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])+Pad_thk)
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP8].pointOn
		s.Line(point1=(DVct[0], DVct[2]+Pad_thk), point2=(DVct2[0], DVct2[2]+Pad_thk))
		DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP5].pointOn
		s.Line(point1=(DVct2[0], DVct2[2]+Pad_thk), point2=(DVct[0], DVct[2]+Pad_thk))
		PlatePath=np.vstack((PlatePath,[0.0315*xloc,DVct2[0], DVct2[2]+Pad_thk],[0.0315*xloc,DVct[0], DVct[2]+Pad_thk]))		
		if xloc==(eval('start%s'%rowname)+eval('end%s'%rowname))/2:
			Setloc[4]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])+Pad_thk)
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP9].pointOn
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP8].pointOn
		DVct3=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP12].pointOn
		s.Line(point1=(DVct[0], DVct[2]+Pad_thk), point2=(-DVct3[0]+Pad, DVct2[2]+Pad_thk))
		s.Line(point1=(Dataxloc[-1,1], -(Dataxloc[-1,2])), point2=(-DVct3[0]+Pad, DVct2[2]+Pad_thk))
		if xloc==(eval('start%s'%rowname)+eval('end%s'%rowname))/2:
			DVct1=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP13].pointOn
			DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP13].pointOn
			DVct3=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP13].pointOn
			DVct4=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP13].pointOn
			DVct5=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP13].pointOn
			A=[[DVct1[0],DVct1[2]+Pad_thk,0.0315*xloc],[DVct2[0],DVct2[2]+Pad_thk,0.0315*xloc],[DVct3[0],DVct3[2]+Pad_thk,0.0315*xloc],[DVct4[0],DVct4[2]+Pad_thk,0.0315*xloc],[DVct5[0],DVct5[2]+Pad_thk,0.0315*xloc]]
		if xloc==eval('start%s'%rowname) or xloc==eval('end%s'%rowname)-1:
			Path=np.hstack((Path,(11.85+xloc*0.0315,-DVct3[0]+Pad, DVct2[2]+Pad_thk)))
			Pathloc=np.vstack((Pathloc,Path))
		p = mdb.models['Model-1'].parts['Elastomer']
		if xloc==0:
			p.BaseWire(sketch=s)
			wire_count+=1
		else:
			p.Wire(sketchPlane=d2[DPP], sketchUpEdge=d2[DPE], sketchPlaneSide=SIDE1, 
				sketchOrientation=RIGHT, sketch=s)
			wire_count+=1
		s.unsetPrimaryObject()
		p = mdb.models['Model-1'].parts['Elastomer']
		del mdb.models['Model-1'].sketches['__profile__']
		p = mdb.models['Model-1'].parts['Elastomer']
		p = mdb.models['Model-1'].parts['Elastomer']
		e = p.edges
		edges = e.getByBoundingBox(-1000000,-1000000,0.0315*xloc-0.0001,1000000,1000000,0.0315*xloc+0.0001)
		p.Set(edges=edges, name='E%i'%xloc)
		exec('Dataxloc%i=Dataxloc'%xloc)
		exec('PlatePath%i=PlatePath'%xloc)
	exec('Setloc%s=Setloc'%rowname)
	for i in range(29):
		j=i+1
		DP=p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=eval('Dataxloc%i[j,1]'%eval('start%s'%rowname)))
		DPP=DP.id
		p = mdb.models['Model-1'].parts['Elastomer']
		d2 = p.datums
		DP=p.DatumAxisByTwoPoint(point1=(Datalow[-1,0],0,0), point2=(Datalow[-1,0],10,0))
		DPE=DP.id
		t = p.MakeSketchTransform(sketchPlane=d2[DPP], sketchUpEdge=d2[DPE], 
			sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0, 
			0, 0))
		s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
			sheetSize=198.85, gridSpacing=4.97, transform=t)
		g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
		s.setPrimaryObject(option=SUPERIMPOSE)
		p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
		s.Spline(points=[(-0.0315*(k), -eval('Dataxloc%i[j,2]'%k)) for k in range(eval('start%s'%rowname),eval('end%s'%rowname))])
		W=p.Wire(sketchPlane=d2[DPP], sketchUpEdge=d2[DPE], sketchPlaneSide=SIDE1, 
			sketchOrientation=RIGHT, sketch=s)
		wire_count+=1
		s.unsetPrimaryObject()
		p.Set(edges=e.getByBoundingBox(eval('Dataxloc%i[j,1]'%eval('start%s'%rowname)),0,0.0315*eval('start%s'%rowname),eval('Dataxloc%i[j,1]'%eval('start%s'%rowname)),100,0.0315*eval('end%s'%rowname)), name='PathWire%s%i'%(rowname,i))

	DP=p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=-0.15)
	DPP=DP.id
	p = mdb.models['Model-1'].parts['Elastomer']
	d2 = p.datums
	DP=p.DatumAxisByTwoPoint(point1=(Datalow[-1,0],0,0), point2=(Datalow[-1,0],10,0))
	DPE=DP.id
	t = p.MakeSketchTransform(sketchPlane=d2[DPP], sketchUpEdge=d2[DPE], 
		sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0, 
		0, 0))
	s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
		sheetSize=198.85, gridSpacing=4.97, transform=t)
	g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
	s.setPrimaryObject(option=SUPERIMPOSE)
	p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
	s.Line(point1=(-(Pathloc[1,0]-11.85), Pathloc[1,2]), point2=(-(Pathloc[2,0]-11.85), Pathloc[2,2]))
	W=p.Wire(sketchPlane=d2[DPP], sketchUpEdge=d2[DPE], sketchPlaneSide=SIDE1, 
		sketchOrientation=RIGHT, sketch=s)
	wire_count+=1
	s.unsetPrimaryObject()
	p.Set(edges=[e.findAt(((Pathloc[1,1], Pathloc[1,2]+0.0315*(loca-0.5)*tan(6.841917290255273E-02), Pathloc[1,0]-11.85+0.0315*(loca-0.5)), ),  ) for loca in range(1,eval('end%s'%rowname)-eval('start%s'%rowname))], name='PathWire%s%i'%(rowname,29))

	DP=p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.15)
	DPP=DP.id
	p = mdb.models['Model-1'].parts['Elastomer']
	p.featuresById[W.id]
	d2 = p.datums
	DP=p.DatumAxisByTwoPoint(point1=(Datalow[-1,0],0,0), point2=(Datalow[-1,0],10,0))
	DPE=DP.id
	t = p.MakeSketchTransform(sketchPlane=d2[DPP], sketchUpEdge=d2[DPE], 
		sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0, 
		0, 0))
	s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
		sheetSize=198.85, gridSpacing=4.97, transform=t)
	g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
	s.setPrimaryObject(option=SUPERIMPOSE)
	p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
	s.Line(point1=(-(Pathloc[1,3]-11.85), Pathloc[1,5]), point2=(-(Pathloc[2,3]-11.85), Pathloc[2,5]))
	W=p.Wire(sketchPlane=d2[DPP], sketchUpEdge=d2[DPE], sketchPlaneSide=SIDE1, 
		sketchOrientation=RIGHT, sketch=s)
	wire_count+=1
	s.unsetPrimaryObject()
	p.Set(edges=[e.findAt(((Pathloc[1,4], Pathloc[1,5]+0.0315*(loca-0.5)*tan(6.841917290255273E-02), Pathloc[1,3]-11.85+0.0315*(loca-0.5)), ),  ) for loca in range(1,eval('end%s'%rowname)-eval('start%s'%rowname))], name='PathWire%s%i'%(rowname,30))	
	p.SolidLoft(loftsections=[p.sets['E%i'%xloc].edges[:] for xloc in range(eval('start%s'%rowname),eval('end%s'%rowname))],paths=[p.sets['PathWire%s%i'%(rowname,locc)].edges[:] for locc in range(31)],globalSmoothing=ON,keepInternalBoundaries=ON)
for rowname in [rownamelist[0],rownamelist[-1]]:
	Pathloc=(0,0,0,0,0,0)
	if rowname==rownamelist[0]:
		begin=start
		finish=eval('start%s'%rowname)+1
	else:
		begin=eval('end%s'%rowname)-1
		finish=end
	for xloc in range(begin,finish):
		if xloc>0:
			DP=p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0315*xloc)
			DPP=DP.id
			exec('DPP%i=DP.id'%xloc)				
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
		Dataxloc=[[0.000000 for x in range(7)] for y in range(len(Datalow))]
		Dataxloc=np.array(Dataxloc)
		for i in range(len(Datalow)):
			if abs(Datalow[i,0]-XLOC)<1E-06:
				Dataxloc[count]=Datalow[i]
				count=count+1
		Dataxloc=Dataxloc[~np.all(Dataxloc == 0, axis=1)]
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2A'].datums[DP13].pointOn
		Pad_thk=DVct2[2]+Datalow[0,2]-0.18+(3.15-0.0315*xloc)*tan(6.841917290255273E-02)
		Pad_thk=(0.0315*xloc)*tan(6.841917290255273E-02)+lowest-DVct2[2]-min_thk
		j=0
		if xloc==(begin+finish)/2:
			for i in range(len(Dataxloc)-1):
				Presloc=np.vstack((Presloc,[0.7*(Dataxloc[i,1])+0.3*(Dataxloc[i+1,1]),0.7*(-(Dataxloc[i,2]))-0.3*(Dataxloc[i+1,2]),0.0315*xloc]))
				Presloc=np.vstack((Presloc,[0.3*(Dataxloc[i,1])+0.7*(Dataxloc[i+1,1]),0.3*(-(Dataxloc[i,2]))-0.7*(Dataxloc[i+1,2]),0.0315*xloc]))
				s.Line(point1=(Dataxloc[i,1], -(Dataxloc[i,2])), point2=(Dataxloc[i+1,1], -(Dataxloc[i+1,2])))				
			Presloc=Presloc[~np.all(Presloc == 0, axis=1)]
		else:		
			for i in range(len(Dataxloc)-1):
				s.Spline(points=((Dataxloc[i,1], -(Dataxloc[i,2])), (Dataxloc[i+1,1], -(Dataxloc[i+1,2]))))				
		DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP12].pointOn
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP8].pointOn
		s.Line(point1=(DVct[0]-Pad, DVct2[2]+Pad_thk), point2=(Dataxloc[0,1], -(Dataxloc[0,2])))
		s.Line(point1=(DVct[0]-Pad, DVct2[2]+Pad_thk), point2=(DVct2[0], DVct2[2]+Pad_thk))
		if xloc==begin or xloc==finish-1:
			Path=(11.85+xloc*0.0315,DVct[0]-Pad, DVct2[2]+Pad_thk)
		DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP5].pointOn
		s.Line(point1=(DVct2[0], DVct2[2]+Pad_thk), point2=(DVct[0], DVct[2]+Pad_thk))
		PlatePath=np.vstack(([0.0315*xloc,DVct2[0], DVct2[2]+Pad_thk],[0.0315*xloc,DVct[0], DVct[2]+Pad_thk]))
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP8].pointOn
		s.Line(point1=(DVct[0], DVct[2]+Pad_thk), point2=(DVct2[0], DVct2[2]+Pad_thk))
		DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP5].pointOn
		s.Line(point1=(DVct2[0], DVct2[2]+Pad_thk), point2=(DVct[0], DVct[2]+Pad_thk))
		PlatePath=np.vstack((PlatePath,[0.0315*xloc,DVct2[0], DVct2[2]+Pad_thk],[0.0315*xloc,DVct[0], DVct[2]+Pad_thk]))
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP8].pointOn
		s.Line(point1=(DVct[0], DVct[2]+Pad_thk), point2=(DVct2[0], DVct2[2]+Pad_thk))
		DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP5].pointOn
		s.Line(point1=(DVct2[0], DVct2[2]+Pad_thk), point2=(DVct[0], DVct[2]+Pad_thk))
		PlatePath=np.vstack((PlatePath,[0.0315*xloc,DVct2[0], DVct2[2]+Pad_thk],[0.0315*xloc,DVct[0], DVct[2]+Pad_thk]))					
		if xloc==begin:
			Axis_start=(DVct2[0], DVct2[2]+Pad_thk,0.0315*xloc)
			Axis_end=(DVct[0], DVct[2]+Pad_thk,0.0315*xloc)
			DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP5].pointOn
			if rowname==rownamelist[0]:
				DPstart=(0.5*(DVct[0]+DVct2[0]), 0.5*(DVct[2]+DVct2[2])+Pad_thk+0.001,0.0315*xloc)
		if xloc==finish-1:
			if rowname==rownamelist[-1]:
				DPend=(0.5*(DVct[0]+DVct2[0]), 0.5*(DVct[2]+DVct2[2])+Pad_thk+0.001,0.0315*xloc)
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP8].pointOn
		s.Line(point1=(DVct[0], DVct[2]+Pad_thk), point2=(DVct2[0], DVct2[2]+Pad_thk))
		DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP5].pointOn
		s.Line(point1=(DVct2[0], DVct2[2]+Pad_thk), point2=(DVct[0], DVct[2]+Pad_thk))
		PlatePath=np.vstack((PlatePath,[0.0315*xloc,DVct2[0], DVct2[2]+Pad_thk],[0.0315*xloc,DVct[0], DVct[2]+Pad_thk]))
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP8].pointOn
		s.Line(point1=(DVct[0], DVct[2]+Pad_thk), point2=(DVct2[0], DVct2[2]+Pad_thk))
		DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP5].pointOn
		s.Line(point1=(DVct2[0], DVct2[2]+Pad_thk), point2=(DVct[0], DVct[2]+Pad_thk))
		PlatePath=np.vstack((PlatePath,[0.0315*xloc,DVct2[0], DVct2[2]+Pad_thk],[0.0315*xloc,DVct[0], DVct[2]+Pad_thk]))
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP9].pointOn
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP8].pointOn
		DVct3=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP12].pointOn
		s.Line(point1=(DVct[0], DVct[2]+Pad_thk), point2=(-DVct3[0]+Pad, DVct2[2]+Pad_thk))
		s.Line(point1=(Dataxloc[-1,1], -(Dataxloc[-1,2])), point2=(-DVct3[0]+Pad, DVct2[2]+Pad_thk))
		if xloc==(begin+finish)/2:
			DVct1=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP13].pointOn
			DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP13].pointOn
			DVct3=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP13].pointOn
			DVct4=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP13].pointOn
			DVct5=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP13].pointOn
			A=[[DVct1[0],DVct1[2]+Pad_thk,0.0315*xloc],[DVct2[0],DVct2[2]+Pad_thk,0.0315*xloc],[DVct3[0],DVct3[2]+Pad_thk,0.0315*xloc],[DVct4[0],DVct4[2]+Pad_thk,0.0315*xloc],[DVct5[0],DVct5[2]+Pad_thk,0.0315*xloc]]
		if xloc==begin or xloc==finish-1:
			Path=np.hstack((Path,(11.85+xloc*0.0315,-DVct3[0]+Pad, DVct2[2]+Pad_thk)))
			Pathloc=np.vstack((Pathloc,Path))
		p = mdb.models['Model-1'].parts['Elastomer']
		if xloc==0:
			p.BaseWire(sketch=s)
			wire_count+=1
		else:
			p.Wire(sketchPlane=d2[DPP], sketchUpEdge=d2[DPE], sketchPlaneSide=SIDE1, 
				sketchOrientation=RIGHT, sketch=s)
			wire_count+=1
		s.unsetPrimaryObject()
		p = mdb.models['Model-1'].parts['Elastomer']
		del mdb.models['Model-1'].sketches['__profile__']
		p = mdb.models['Model-1'].parts['Elastomer']
		p = mdb.models['Model-1'].parts['Elastomer']
		e = p.edges
		f=p.faces
		edges = e.getByBoundingBox(-1000000,-1000000,0.0315*xloc-0.0001,1000000,1000000,0.0315*xloc+0.0001)
		p.Set(edges=edges, name='E%i'%xloc)
		exec('Dataxloc%i=Dataxloc'%xloc)
		exec('PlatePath%i=PlatePath'%xloc)
	p.SolidLoft(loftsections=[p.sets['E%i'%xloc].edges[:] for xloc in range(begin,finish)],globalSmoothing=OFF,keepInternalBoundaries=ON)
faces1=f.findAt((DPstart, ))
faces2=f.findAt((DPend, ))
p.Set(faces=faces1+faces2, name='StartEnd')	
for rowname in rownamelist:
	for locc in [eval('start%s'%rowname),eval('end%s'%rowname)-1]:
		edges = e.getByBoundingBox(-1000000,-1000000,0.0315*locc-0.0001,1000000,1000000,0.0315*locc+0.0001)
		p.Set(edges=edges, name='E%i'%locc)	
f=p.faces
for i in range(ND-1):	
	p.SolidLoft(loftsections=[p.sets['E%i'%xloc].edges[:] for xloc in [eval('end%s'%chr(65+i))-1,eval('start%s'%chr(65+i+1))]],globalSmoothing=ON)
for i in range(ND-1):	
	exec('LB=min(-Dataxloc%i[:,2])'%(eval('end%s'%chr(65+i))-1))
	if i==0:
		faces = f.getByBoundingBox(-1000000,LB-0.0001,0.0315*eval('end%s'%chr(65+i))-1-0.0001,1000000,1000000,0.0315*eval('start%s'%chr(65+i+1))+0.0001)
	else:
		faces += f.getByBoundingBox(-1000000,LB-0.0001,0.0315*eval('end%s'%chr(65+i))-1-0.0001,1000000,1000000,0.0315*eval('start%s'%chr(65+i+1))+0.0001)
p.Surface(side1Faces=faces, name='Preslocface')
a = mdb.models['Model-1'].rootAssembly
p=mdb.models['Model-1'].parts['Elastomer']
s1 = mdb.models['Model-1'].parts['Elastomer'].faces
p.Surface(side1Faces=[s1.findAt(((Presloc[i,0],Presloc[i,1],Presloc[i,2]), ),  ) for i in range(len(Presloc))], name='Elastomer_Surf')	
	
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Elastomer']
a.Instance(name='Elastomer', part=p, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Elastomer', ), axisPoint=[0,0,0], 
	axisDirection=(Axis_end[0]-Axis_start[0],0,0), angle=90+6.841917290255273E-02*180/pi)
a.rotate(instanceList=('Elastomer', ), axisPoint=[0,0,0], 
	axisDirection=(0,0,1), angle=180)	
for loc in range(3):
	Vct[loc]= 0.5*(mdb.models['Model-1'].rootAssembly.instances['Middle_PlateA'].datums[DP5].pointOn[loc]+mdb.models['Model-1'].rootAssembly.instances['Middle_PlateA'].datums[DP8].pointOn[loc])-mdb.models['Model-1'].rootAssembly.instances['Elastomer'].datums[DPloc].pointOn[loc]
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('Elastomer', ), vector=(Vct[0],Vct[1],Vct[2]))

DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2A'].datums[DP12].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2A'].datums[DP8].pointOn
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Elastomer'].faces
faces1 = f1.findAt(((DVct[0]-Pad,DVct[1],DVct2[2]+0.001), ))
faces2 = f1.findAt(((-DVct[0]+Pad,DVct[1],DVct2[2]+0.001), ))
faces3 = f1.findAt(((-DVct[0]+Pad,DVct[1]-0.0315*0.5*(startA-start),DVct2[2]+0.001), ))
faces4 = f1.findAt(((DVct[0]-Pad,DVct[1]-0.0315*0.5*(startA-start),DVct2[2]+0.001), ))	
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rownamelist[-1]].datums[DP11].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rownamelist[-1]].datums[DP7].pointOn
faces5 = f1.findAt(((-DVct[0]+Pad,DVct[1]+0.0315*0.5*(end-eval('end%s'%rownamelist[-1])),DVct2[2]+0.001), ))
faces6 = f1.findAt(((DVct[0]-Pad,DVct[1]+0.0315*0.5*(end-eval('end%s'%rownamelist[-1])),DVct2[2]+0.001), ))
region = a.Set(faces=faces1+faces2+faces3+faces4+faces5+faces6, name='Elast_ends')

for rowname in rownamelist:
	p = mdb.models['Model-1'].parts['Elastomer']
	s = p.faces 
	side1Faces = s.findAt(((eval('Setloc%s'%rowname)[0][1], eval('Setloc%s'%rowname)[0][2], eval('Setloc%s'%rowname)[0][0]), ))
	p.Surface(side1Faces=side1Faces, name='Surf_R2%s'%rowname)
	side1Faces = s.findAt(((eval('Setloc%s'%rowname)[1][1], eval('Setloc%s'%rowname)[1][2], eval('Setloc%s'%rowname)[1][0]), ))
	p.Surface(side1Faces=side1Faces, name='Surf_R%s'%rowname)
	side1Faces = s.findAt(((eval('Setloc%s'%rowname)[2][1], eval('Setloc%s'%rowname)[2][2], eval('Setloc%s'%rowname)[2][0]), ))
	p.Surface(side1Faces=side1Faces, name='Surf_M%s'%rowname)
	side1Faces = s.findAt(((eval('Setloc%s'%rowname)[3][1], eval('Setloc%s'%rowname)[3][2], eval('Setloc%s'%rowname)[3][0]), ))
	p.Surface(side1Faces=side1Faces, name='Surf_L%s'%rowname)
	side1Faces = s.findAt(((eval('Setloc%s'%rowname)[4][1], eval('Setloc%s'%rowname)[4][2], eval('Setloc%s'%rowname)[4][0]), ))
	p.Surface(side1Faces=side1Faces, name='Surf_L2%s'%rowname)	

# Create Material
print 'Creating the Materials'
# Create Aluminium
mdb.models['Model-1'].Material(name='Aluminium')
mdb.models['Model-1'].materials['Aluminium'].Elastic(table=((70E09, 0.3), ))
mdb.models['Model-1'].materials['Aluminium'].Density(table=((6500000.0, ), ))

#Create Elastomer
mdb.models['Model-1'].Material(name='Elastomer')
mdb.models['Model-1'].materials['Elastomer'].Hyperelastic(
    materialType=ISOTROPIC, testData=OFF, type=MOONEY_RIVLIN, 
    volumetricResponse=VOLUMETRIC_DATA, table=((0.3339E06,-0.000337E06, 0.0015828E-06), ))
mdb.models['Model-1'].materials['Elastomer'].Density(table=((1522000.0, ), ))

## Create Material##
mdb.models['Model-1'].Material(name='TT_matl')
mdb.models['Model-1'].materials['TT_matl'].Density(table=((6500000.0, ), ))
mdb.models['Model-1'].materials['TT_matl'].Expansion(table=((1.0, ), ))
mdb.models['Model-1'].materials['TT_matl'].expansion.setValues(
    type=ANISOTROPIC, table=((0.0, 0.0, 0.0, 0.0, 0.0, -0.1), ), zero=0.0)
mdb.models['Model-1'].materials['TT_matl'].Elastic(table=((83E09, 0.3), ))

mdb.models['Model-1'].Material(name='TT_matl1')
mdb.models['Model-1'].materials['TT_matl1'].Density(table=((6500000.0, ), ))
mdb.models['Model-1'].materials['TT_matl1'].Expansion(table=((1.0, ), ))
mdb.models['Model-1'].materials['TT_matl1'].expansion.setValues(
    type=ANISOTROPIC, table=((0.0, 0.0, 0.0, 0.0, 0.0, 0.1), ), zero=0.0)
mdb.models['Model-1'].materials['TT_matl1'].Elastic(table=((83E09, 0.3), ))

#Create/Assign Section
print 'Creating the Sections'
mdb.models['Model-1'].HomogeneousShellSection(name='Solid', 
	preIntegrate=OFF, material='Aluminium', thicknessType=UNIFORM, 
	thickness=1.0, thicknessField='', idealization=NO_IDEALIZATION, 
	poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
	useDensity=OFF, integrationRule=SIMPSON, numIntPts=5)

#Create/Assign Section
print 'Creating the Sections'
mdb.models['Model-1'].HomogeneousSolidSection(name='Stretchable', 
	material='Elastomer', thickness=None)	
print 'Assigning the Sections'

for rowname in rownamelist:
	for plate_num in [1,2,3,4,5]:
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		f = p.faces
		faces = f.getByBoundingBox(-eval('H%i%s'%(plate_num,rowname)),-eval('V%s'%rowname),-eval('G%s'%rowname),eval('H%i%s'%(plate_num,rowname)),eval('V%s'%rowname),eval('G%s'%rowname))
		region = p.Set(faces=faces, name='Set-1')
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		p.SectionAssignment(region=region, sectionName='Solid', offset=0.0, 
			offsetType=TOP_SURFACE, offsetField='', 
			thicknessAssignment=FROM_SECTION)

p = mdb.models['Model-1'].parts['Elastomer']
c = p.cells
cells=c[:]
#faces = f.getByBoundingBox(-3*H,-3V,-G,H,V,G)
region = p.Set(cells=cells, name='Set-1')
p = mdb.models['Model-1'].parts['Elastomer']
p.SectionAssignment(region=region, sectionName='Stretchable', offset=0.0, 
	offsetType=TOP_SURFACE, offsetField='', 
	thicknessAssignment=FROM_SECTION)

## Create Section TT##
mdb.models['Model-1'].HomogeneousSolidSection(name='TT_sect', material='TT_matl', 
    thickness=None)
mdb.models['Model-1'].HomogeneousSolidSection(name='TT_sect1', material='TT_matl1', 
    thickness=None)
p = mdb.models['Model-1'].parts['Torque_Tube']
v1, e = p.vertices, p.edges
DP=p.DatumCsysByThreePoints(point1=v1.findAt(coordinates=(ORTT, 0.0,-(0.5*TTL))), 
    name='Datum csys-1', coordSysType=CYLINDRICAL, origin=(0,0,-(0.5*TTL)), 
    point2=(ORTT,ORTT,-(0.5*TTL)))
DPCYS=DP.id
DP=p.DatumCsysByThreePoints(point1=v1.findAt(coordinates=(ORTT, 0.0,(0.5*TTL))), 
    name='Datum csys-2', coordSysType=CYLINDRICAL, origin=(0,0,(0.5*TTL)), 
    point2=(ORTT,ORTT,(0.5*TTL)))
DPCYS2=DP.id
p = mdb.models['Model-1'].parts['Torque_Tube']
c = p.cells
cells = c.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), 0), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), 
    -0.5*(IRTT+ORTT)*cos(pi/4), 0), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 0), ), ((
    0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4),0), ))
region = regionToolset.Region(cells=cells)
orientation = mdb.models['Model-1'].parts['Torque_Tube'].datums[DPCYS2]
mdb.models['Model-1'].parts['Torque_Tube'].MaterialOrientation(region=region, 
    orientationType=SYSTEM, axis=AXIS_3, localCsys=orientation, fieldName='', 
    additionalRotationType=ROTATION_NONE, angle=0.0, 
    additionalRotationField='', stackDirection=STACK_3)

p = mdb.models['Model-1'].parts['Torque_Tube']
c = p.cells
cells = c.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), 0), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), 
    -0.5*(IRTT+ORTT)*cos(pi/4), 0), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 0), ), ((
    0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4),0), ))
region = p.Set(cells=cells, name='TT')
p = mdb.models['Model-1'].parts['Torque_Tube']
p.SectionAssignment(region=region, sectionName='TT_sect', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)
	
#Create Wires
print 'Defining Wires'
for rowname in rownamelist:
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	v1 = a.instances['Middle_Plate%s'%rowname].vertices
	v2 = a.instances['Right_Plate%s'%rowname].vertices
	v3 = a.instances['Left_Plate%s'%rowname].vertices
	v4 = a.instances['Right_Plate_2%s'%rowname].vertices
	v5 = a.instances['Left_Plate_2%s'%rowname].vertices
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP2].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[eval('RP1%s'%rowname)], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP3].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[eval('RP2%s'%rowname)], v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP1].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[eval('RP3%s'%rowname)], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP4].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[eval('RP4%s'%rowname)], v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP2].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[eval('RP6%s'%rowname)], v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP1].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[eval('RP8%s'%rowname)], v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP4].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[eval('RP7%s'%rowname)], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP3].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[eval('RP5%s'%rowname)], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP4].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[eval('RP16%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)			
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP3].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[eval('RP14%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)				
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP2].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[eval('RP10%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)				
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP1].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[eval('RP12%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)				
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP3].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[eval('RP9%s'%rowname)], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)	
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP4].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[eval('RP11%s'%rowname)], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)	
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP2].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[eval('RP13%s'%rowname)], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)		
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP1].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[eval('RP15%s'%rowname)], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)			
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP2].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[eval('RP20%s'%rowname)], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP4].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[eval('RP23%s'%rowname)], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)

a = mdb.models['Model-1'].rootAssembly
e1 = a.edges 
edges1 = e1.getByBoundingBox(-LMD,-LMD,-LMD,LMD,LMD,LMD)
a.Set(edges=edges1, name='Wire-Beams')

a = mdb.models['Model-1'].rootAssembly
for wires in range(18*ND):
	a.suppressFeatures(('Wire-%i'% (wires+1), ))

for rowname in rownamelist:
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r11[eval('RP1%s'%rowname)], r11[eval('RP2%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)
	a = mdb.models['Model-1'].rootAssembly
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP10].pointOn
	e1 = a.edges
	edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
	a.Set(edges=edges1, name='Wire-1%s'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP9].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r11[eval('RP3%s'%rowname)], r11[eval('RP4%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)
	a = mdb.models['Model-1'].rootAssembly
	e1 = a.edges
	edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
	a.Set(edges=edges1, name='Wire-2%s'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP11].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r11[eval('RP5%s'%rowname)], r11[eval('RP6%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)
	a = mdb.models['Model-1'].rootAssembly
	e1 = a.edges
	edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
	a.Set(edges=edges1, name='Wire-3%s'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP12].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r11[eval('RP7%s'%rowname)], r11[eval('RP8%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)
	a = mdb.models['Model-1'].rootAssembly
	e1 = a.edges
	edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
	a.Set(edges=edges1, name='Wire-4%s'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP10].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r11[eval('RP9%s'%rowname)], r11[eval('RP10%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)
	a = mdb.models['Model-1'].rootAssembly
	e1 = a.edges
	edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
	a.Set(edges=edges1, name='Wire-5%s'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP9].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r11[eval('RP11%s'%rowname)], r11[eval('RP12%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)
	a = mdb.models['Model-1'].rootAssembly
	e1 = a.edges
	edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
	a.Set(edges=edges1, name='Wire-6%s'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP11].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r11[eval('RP13%s'%rowname)], r11[eval('RP14%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)
	a = mdb.models['Model-1'].rootAssembly
	e1 = a.edges
	edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
	a.Set(edges=edges1, name='Wire-7%s'%rowname)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP12].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r11[eval('RP15%s'%rowname)], r11[eval('RP16%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)
	a = mdb.models['Model-1'].rootAssembly
	e1 = a.edges
	edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
	a.Set(edges=edges1, name='Wire-8%s'%rowname)

	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP9].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r11[eval('RP21%s'%rowname)], r11[eval('RP27%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)
	a = mdb.models['Model-1'].rootAssembly
	e1 = a.edges
	edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
	a.Set(edges=edges1, name='Wire-9%s'%rowname)

	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP10].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r11[eval('RP20%s'%rowname)], r11[eval('RP26%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)
	a = mdb.models['Model-1'].rootAssembly
	e1 = a.edges
	edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
	a.Set(edges=edges1, name='Wire-10%s'%rowname)

	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP11].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r11[eval('RP22%s'%rowname)], r11[eval('RP28%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)
	a = mdb.models['Model-1'].rootAssembly
	e1 = a.edges
	edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
	a.Set(edges=edges1, name='Wire-11%s'%rowname)

	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP12].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r11[eval('RP23%s'%rowname)], r11[eval('RP29%s'%rowname)]), ), mergeType=IMPRINT, meshable=OFF)
	a = mdb.models['Model-1'].rootAssembly
	e1 = a.edges
	edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
	a.Set(edges=edges1, name='Wire-12%s'%rowname)
	a=mdb.models['Model-1'].rootAssembly
	a.SetByBoolean(name='MRH%s'%rowname, sets=(a.sets['Wire-1%s'%rowname], a.sets['Wire-2%s'%rowname], ))
	a=mdb.models['Model-1'].rootAssembly
	a.SetByBoolean(name='MLH%s'%rowname, sets=(a.sets['Wire-3%s'%rowname], a.sets['Wire-4%s'%rowname], ))
	a=mdb.models['Model-1'].rootAssembly
	a.SetByBoolean(name='RRH%s'%rowname, sets=(a.sets['Wire-5%s'%rowname], a.sets['Wire-6%s'%rowname], ))
	a=mdb.models['Model-1'].rootAssembly
	a.SetByBoolean(name='LLH%s'%rowname, sets=(a.sets['Wire-7%s'%rowname], a.sets['Wire-8%s'%rowname], ))

a = mdb.models['Model-1'].rootAssembly
for wires in range(68):
	a.resumeFeatures(('Wire-%i'% (wires+1), ))
	
#Create and Assign Connector Sections
print 'Create Connector Section'	
mdb.models['Model-1'].ConnectorSection(name='HINGE', assembledType=HINGE)

print 'Create Connector Section'	
mdb.models['Model-1'].ConnectorSection(name='BEAM', assembledType=BEAM)

print 'Define Datum Coordinate Systems'
for rowname in rownamelist:
	a = mdb.models['Model-1'].rootAssembly
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP10].pointOn
	a.DatumCsysByThreePoints(name='CSYS_MR%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-eval('V%s'%rowname), DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for Middle to Right 1
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP9].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_MR2%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-eval('V%s'%rowname), DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for Middle to R2
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP11].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_ML%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+eval('V%s'%rowname), DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for Middle to Left
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP12].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_ML2%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+eval('V%s'%rowname), DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for Middle to L2
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP10].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_RR%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-eval('V%s'%rowname), DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for  RHS Support
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP9].pointOn	
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_RR2%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-eval('V%s'%rowname), DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for RHS Support 2
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP11].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_LL%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+eval('V%s'%rowname), DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP12].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_LL2%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+eval('V%s'%rowname), DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP9].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_RS%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+eval('V%s'%rowname), DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP10].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_RS2%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+eval('V%s'%rowname), DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP11].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_LS%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+eval('V%s'%rowname), DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP12].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_LS2%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+eval('V%s'%rowname), DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2

print 'Assign Connector Sections'
for rowname in rownamelist:
	a = mdb.models['Model-1'].rootAssembly
	region=a.sets['Wire-1%s'%rowname]
	datum1 = a.datums[a.features['CSYS_MR%s'%rowname].id]
	csa = a.SectionAssignment(sectionName='HINGE', region=region)
	a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

	a = mdb.models['Model-1'].rootAssembly
	region=a.sets['Wire-2%s'%rowname]
	datum1 = a.datums[a.features['CSYS_MR2%s'%rowname].id]
	csa = a.SectionAssignment(sectionName='HINGE', region=region)
	a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)
	a = mdb.models['Model-1'].rootAssembly
	region=a.sets['Wire-3%s'%rowname]
	datum1 = a.datums[a.features['CSYS_ML%s'%rowname].id]
	csa = a.SectionAssignment(sectionName='HINGE', region=region)
	a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

	a = mdb.models['Model-1'].rootAssembly
	region=a.sets['Wire-4%s'%rowname]
	datum1 = a.datums[a.features['CSYS_ML2%s'%rowname].id]
	csa = a.SectionAssignment(sectionName='HINGE', region=region)
	a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

	a = mdb.models['Model-1'].rootAssembly
	region=a.sets['Wire-5%s'%rowname]
	datum1 = a.datums[a.features['CSYS_RR%s'%rowname].id]
	csa = a.SectionAssignment(sectionName='HINGE', region=region)
	a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

	a = mdb.models['Model-1'].rootAssembly
	region=a.sets['Wire-6%s'%rowname]
	datum1 = a.datums[a.features['CSYS_RR2%s'%rowname].id]
	csa = a.SectionAssignment(sectionName='HINGE', region=region)
	a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

	a = mdb.models['Model-1'].rootAssembly
	region=a.sets['Wire-7%s'%rowname]
	datum1 = a.datums[a.features['CSYS_LL%s'%rowname].id]
	csa = a.SectionAssignment(sectionName='HINGE', region=region)
	a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

	a = mdb.models['Model-1'].rootAssembly
	region=a.sets['Wire-8%s'%rowname]
	datum1 = a.datums[a.features['CSYS_LL2%s'%rowname].id]
	csa = a.SectionAssignment(sectionName='HINGE', region=region)
	a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

	a = mdb.models['Model-1'].rootAssembly
	region=a.sets['Wire-9%s'%rowname]
	datum1 = a.datums[a.features['CSYS_RS%s'%rowname].id]
	csa = a.SectionAssignment(sectionName='HINGE', region=region)
	a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

	a = mdb.models['Model-1'].rootAssembly
	region=a.sets['Wire-10%s'%rowname]
	datum1 = a.datums[a.features['CSYS_RS2%s'%rowname].id]
	csa = a.SectionAssignment(sectionName='HINGE', region=region)
	a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

	a = mdb.models['Model-1'].rootAssembly
	region=a.sets['Wire-11%s'%rowname]
	datum1 = a.datums[a.features['CSYS_LS%s'%rowname].id]
	csa = a.SectionAssignment(sectionName='HINGE', region=region)
	a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

	a = mdb.models['Model-1'].rootAssembly
	region=a.sets['Wire-12%s'%rowname]
	datum1 = a.datums[a.features['CSYS_LS2%s'%rowname].id]
	csa = a.SectionAssignment(sectionName='HINGE', region=region)
	a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-Beams']
csa = a.SectionAssignment(sectionName='BEAM', region=region)

#Mesh Parts
print 'Meshing the Part'
for rowname in rownamelist:
	for plate_num in [1,2,3,4,5]:	
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		p.seedPart(size=0.01, deviationFactor=0.1, minSizeFactor=0.1)		#G/6
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		p.generateMesh()

p = mdb.models['Model-1'].parts['Elastomer']
d2 = p.datums
pickedCells=p.cells[:]
Partition_xloc=np.linspace(start+1,startA-1,(startA-start-1))
for rowname in rownamelist:
	str=eval('start%s'%rowname)
	fin=eval('end%s'%rowname)
	Partition_xloc=np.append(Partition_xloc,np.linspace(str+1,fin-2,(fin-str-2)))
Partition_xloc=np.append(Partition_xloc,np.linspace(fin,end-2,(end-fin-1)))

for xloc in Partition_xloc:
	pickedCells=p.cells[:]
	mdb.models['Model-1'].parts['Elastomer'].PartitionCellByDatumPlane(cells=pickedCells,datumPlane=d2[eval('DPP%i'%xloc)])

v=a.instances['Elastomer'].vertices
v=p.vertices
p.Set(vertices=[v.findAt(((Curveloc[l,0],Curveloc[l,1],Curveloc[l,2]), ),  ) for l in range(len(Curveloc))], name='TIPNODE')
verts=a.instances['Elastomer'].sets['TIPNODE'].vertices[:]
a.Set(vertices=verts, name='TIPNODE')

elemType1 = mesh.ElemType(elemCode=C3D8H, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=C3D6H, elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D4H, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['Elastomer']
p.seedPart(size=0.01, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['Elastomer']
c = p.cells
cells = c[:]
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
	elemType3))
p.setMeshControls(regions=cells, elemShape=TET, technique=FREE)
p = mdb.models['Model-1'].parts['Elastomer']
p.generateMesh()

## Mesh##
elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD, 
    kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, 
    hourglassControl=DEFAULT, distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['Torque_Tube']
c = p.cells
cells = c.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), 0), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), 
    -0.5*(IRTT+ORTT)*cos(pi/4), 0), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 0), ), ((
    0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4),0), ))
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
    elemType3))
p = mdb.models['Model-1'].parts['Torque_Tube']
e = p.edges
pickedEdges = p.sets['Longitudinal_edges'].edges[:]
p.seedEdgeByNumber(edges=pickedEdges, number=15, constraint=FIXED)
p = mdb.models['Model-1'].parts['Torque_Tube']
e = p.edges
pickedEdges = p.sets['Radial_edges'].edges[:]
p.seedEdgeByNumber(edges=pickedEdges, number=3, constraint=FIXED)

p = mdb.models['Model-1'].parts['Torque_Tube']
e = p.edges
pickedEdges = p.sets['Circumferential_edges'].edges[:]
p.seedEdgeByNumber(edges=pickedEdges, number=4, constraint=FIXED)
p = mdb.models['Model-1'].parts['Torque_Tube']
p.generateMesh()

#Assemble TTs

p = mdb.models['Model-1'].parts['Torque_Tube']
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
for rowname in rownamelist:
	a.Instance(name='Torque_Tube_3%s'%rowname, part=p, dependent=ON)


## Constraints##
for rowname in rownamelist:
	a = mdb.models['Model-1'].rootAssembly
	region1=a.sets['RP_Z1-%s'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	f1 = a.instances['Torque_Tube_1%s'%rowname].faces
	faces1 = f1.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*TTL)), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 
		-(0.5*TTL)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*TTL)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*TTL)), ))
	region2=a.Set(faces=faces1, name='s_Set-1-%s'%rowname)
	mdb.models['Model-1'].Coupling(name='TT1_KC-%s'%rowname, controlPoint=region1, 
		surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
		localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

	a = mdb.models['Model-1'].rootAssembly
	region1=a.sets['RP_Z1+%s'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	e1 = a.instances['Torque_Tube_1%s'%rowname].edges
	orientation = mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_1%s'%rowname].datums[DPCYS2]
	edges1 = e1.findAt(((-ORTT*cos(pi/4), ORTT*cos(pi/4), (0.5*TTL)), ), ((-ORTT*cos(pi/4), -ORTT*cos(pi/4), 
		(0.5*TTL)), ), ((ORTT*cos(pi/4), ORTT*cos(pi/4), (0.5*TTL)), ), ((ORTT*cos(pi/4), -ORTT*cos(pi/4), (0.5*TTL)), ))
	region2=a.Set(edges=edges1, name='s_Set-1+%s'%rowname)
	mdb.models['Model-1'].Coupling(name='TT1_KC+%s'%rowname, controlPoint=region1, 
		surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
		localCsys=orientation, u1=ON, u2=ON, u3=OFF, ur1=ON, ur2=ON, ur3=ON)

	a = mdb.models['Model-1'].rootAssembly
	region1=a.sets['RP_Z3-%s'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	orientation = mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_3%s'%rowname].datums[DPCYS]
	f1 = a.instances['Torque_Tube_3%s'%rowname].faces
	faces1 = f1.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*TTL)), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 
		-(0.5*TTL)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*TTL)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*TTL)), ))
	region2=a.Set(faces=faces1,name='s_Set-3-%s'%rowname)
	mdb.models['Model-1'].Coupling(name='TT3_KC-%s'%rowname, controlPoint=region1, 
		surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
		localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
		
	a = mdb.models['Model-1'].rootAssembly
	region1=a.sets['RP_Z3+%s'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	orientation = mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_3%s'%rowname].datums[DPCYS2]
	e1 = a.instances['Torque_Tube_3%s'%rowname].edges
	edges1 = e1.findAt(((-ORTT*cos(pi/4), ORTT*cos(pi/4), (0.5*TTL)), ), ((-ORTT*cos(pi/4), -ORTT*cos(pi/4), 
		(0.5*TTL)), ), ((ORTT*cos(pi/4), ORTT*cos(pi/4), (0.5*TTL)), ), ((ORTT*cos(pi/4), -ORTT*cos(pi/4), (0.5*TTL)), ))
	region2=a.Set(edges=edges1, name='s_Set-3+%s'%rowname)
	mdb.models['Model-1'].Coupling(name='TT3_KC+%s'%rowname, controlPoint=region1, 
		surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
		localCsys=orientation, u1=ON, u2=ON, u3=OFF, ur1=ON, ur2=ON, ur3=ON)

	a = mdb.models['Model-1'].rootAssembly
	a.rotate(instanceList=('Torque_Tube_1%s'%rowname,'Torque_Tube_3%s'%rowname, ), axisPoint=(-ORTT,0, 0), 
		axisDirection=(ORTT,0, 0), angle=(90))		
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP12].pointOn
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_1%s'%rowname].datums[DP_TT1].pointOn
	a.translate(instanceList=('Torque_Tube_1%s'%rowname, ), vector=(DVct[0]-DVct2[0], DVct[1]-DVct2[1], DVct[2]-DVct2[2]))
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP12].pointOn

	DVct2=mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_3%s'%rowname].datums[DP_TT2].pointOn
	a.rotate(instanceList=('Torque_Tube_3%s'%rowname, ), axisPoint=(0,0,-ORTT), 
	axisDirection=(0,0,ORTT), angle=(180))	
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP10].pointOn	

	a.translate(instanceList=('Torque_Tube_3%s'%rowname, ), vector=(DVct[0]-DVct2[0], DVct[1]-DVct2[1], DVct[2]-DVct2[2]))
for rowname in rownamelist:
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Right_Plate%s'%rowname].sets['Set-1']
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refPoints1=(r1[eval('RP17%s'%rowname)], )
	region1=regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].RigidBody(name='RB_R%s'%rowname, refPointRegion=region1, bodyRegion=region2)

	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Left_Plate%s'%rowname].sets['Set-1']
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refPoints1=(r1[eval('RP18%s'%rowname)], )
	region1=regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].RigidBody(name='RB_L%s'%rowname, refPointRegion=region1, bodyRegion=region2)

	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Middle_Plate%s'%rowname].sets['Set-1']
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refPoints1=(r1[eval('RP19%s'%rowname)], )
	region1=regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].RigidBody(name='RB_M%s'%rowname, refPointRegion=region1, bodyRegion=region2)

	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Right_Plate_2%s'%rowname].sets['Set-1']
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refPoints1=(r1[eval('RP24%s'%rowname)], )
	region1=regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].RigidBody(name='RB_R_2%s'%rowname, refPointRegion=region1, bodyRegion=region2)

	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Left_Plate_2%s'%rowname].sets['Set-1']
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refPoints1=(r1[eval('RP25%s'%rowname)], )
	region1=regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].RigidBody(name='RB_L_2%s'%rowname, refPointRegion=region1, bodyRegion=region2)
for rowname in rownamelist:
	for plate_num in [1,2,3,4,5]:
		p = mdb.models['Model-1'].parts['Plate%i%s'%(plate_num,rowname)]
		s = p.elements
		side2Elements = s[:]
		p.Surface(side2Elements=side2Elements, name='Surf-Bottom')

for rowname in rownamelist:
	a = mdb.models['Model-1'].rootAssembly
	p=mdb.models['Model-1'].parts['Elastomer']
	region1=a.instances['Elastomer'].surfaces['Surf_L2%s'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Left_Plate_2%s'%rowname].surfaces['Surf-Bottom']
	mdb.models['Model-1'].Tie(name='Tie_L2%s'%rowname, master=region2, slave=region1, 
		positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

	a = mdb.models['Model-1'].rootAssembly
	p=mdb.models['Model-1'].parts['Elastomer']
	region1=a.instances['Elastomer'].surfaces['Surf_L%s'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Left_Plate%s'%rowname].surfaces['Surf-Bottom']
	mdb.models['Model-1'].Tie(name='Tie_L%s'%rowname, master=region2, slave=region1, 
		positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
			
	a = mdb.models['Model-1'].rootAssembly
	p=mdb.models['Model-1'].parts['Elastomer']
	region1=a.instances['Elastomer'].surfaces['Surf_M%s'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Middle_Plate%s'%rowname].surfaces['Surf-Bottom']
	mdb.models['Model-1'].Tie(name='Tie_M%s'%rowname, master=region2, slave=region1, 
		positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
		
	a = mdb.models['Model-1'].rootAssembly
	p=mdb.models['Model-1'].parts['Elastomer']
	region1=a.instances['Elastomer'].surfaces['Surf_R%s'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Right_Plate%s'%rowname].surfaces['Surf-Bottom']
	mdb.models['Model-1'].Tie(name='Tie_R%s'%rowname, master=region2, slave=region1, 
		positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
		
	a = mdb.models['Model-1'].rootAssembly
	p=mdb.models['Model-1'].parts['Elastomer']
	region1=a.instances['Elastomer'].surfaces['Surf_R2%s'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Right_Plate_2%s'%rowname].surfaces['Surf-Bottom']
	mdb.models['Model-1'].Tie(name='Tie_R2%s'%rowname, master=region2, slave=region1, 
		positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

#Define Steps
print 'Defining the Steps'

mdb.models['Model-1'].ImplicitDynamicsStep(name='RBM', previous='Initial', 
    application=QUASI_STATIC, nohaf=OFF, amplitude=RAMP, alpha=DEFAULT, 
    initialConditions=OFF,timePeriod=10.0)
mdb.models['Model-1'].steps['RBM'].setValues(nlgeom=ON)	
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'E', 'PE', 'LE', 'U', 'RF', 'CF','NT','COORD','THE','EE'))
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(
    timeInterval=1.0)
	
#Define BCs				
print 'Defining all BCs'
region = a.sets['Elast_ends']
mdb.models['Model-1'].EncastreBC(name='Fix_Elast', createStepName='Initial', 
	region=region, localCsys=None)
p=mdb.models['Model-1'].parts['Elastomer']
f2=p.faces
faces7=p.sets['StartEnd'].faces[:]
region = a.instances['Elastomer'].sets['StartEnd']
mdb.models['Model-1'].EncastreBC(name='Fix_Elast2', createStepName='Initial', 
	region=region, localCsys=None)
		
for rowname in rownamelist:
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refPoints1=(r1[eval('RP27%s'%rowname)], r1[eval('RP26%s'%rowname)], r1[eval('RP28%s'%rowname)], r1[eval('RP29%s'%rowname)], )
	region = regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].DisplacementBC(name='Fixed hinge Support%s'%rowname, createStepName='Initial', 
		region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
		amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
		localCsys=None)
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refPoints1=(r1[eval('RP21%s'%rowname)], r1[eval('RP22%s'%rowname)], )
	region = regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].DisplacementBC(name='Fixed TT%s'%rowname, createStepName='Initial', 
		region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
		amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
		localCsys=None)


## BCs##
mdb.models['Model-1'].SmoothStepAmplitude(name='Amp-1', timeSpan=STEP, data=((
    0.0, 0.0), (2.0, 0.2), (4.0, 0.4), (6.0, 0.6), (8.0, 0.8), (10.0, 1.0)))
	
for rowname in rownamelist:	
	datum = mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_1%s'%rowname].datums[DPCYS]
	mdb.models['Model-1'].ExpressionField(name='AnalyticalField-1%s'%rowname, 
		localCsys=datum, description='', expression='R/%f' % ORTT)
	a = mdb.models['Model-1'].rootAssembly
	region = a.instances['Torque_Tube_1%s'%rowname].sets['TT']
	mdb.models['Model-1'].Temperature(name='Predefined Field-1%s'%rowname, 
		createStepName='Initial', region=region, distributionType=UNIFORM, 
		crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(0.0, ))
	mdb.models['Model-1'].predefinedFields['Predefined Field-1%s'%rowname].setValues(
		distributionType=FIELD, field='AnalyticalField-1%s'%rowname)
	mdb.models['Model-1'].predefinedFields['Predefined Field-1%s'%rowname].setValuesInStep(
		stepName='RBM', amplitude='Amp-1', magnitudes=(eval('T1%s'%rowname), ))
	datum = mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_3%s'%rowname].datums[DPCYS]
	mdb.models['Model-1'].ExpressionField(name='AnalyticalField-3%s'%rowname, 
		localCsys=datum, description='', expression='R/%f' % ORTT)
	a = mdb.models['Model-1'].rootAssembly
	region = a.instances['Torque_Tube_3%s'%rowname].sets['TT']
	mdb.models['Model-1'].Temperature(name='Predefined Field-3%s'%rowname, 
		createStepName='Initial', region=region, distributionType=UNIFORM, 
		crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(0.0, ))
	mdb.models['Model-1'].predefinedFields['Predefined Field-3%s'%rowname].setValues(
		distributionType=FIELD, field='AnalyticalField-3%s'%rowname)
	mdb.models['Model-1'].predefinedFields['Predefined Field-3%s'%rowname].setValuesInStep(
		stepName='RBM', amplitude='Amp-1', magnitudes=(eval('T3%s'%rowname), ))	
	
# Pressure Load

a = mdb.models['Model-1'].rootAssembly
p=mdb.models['Model-1'].parts['Elastomer']
s1 = mdb.models['Model-1'].parts['Elastomer'].faces
region1 = a.instances['Elastomer'].surfaces['Elastomer_Surf']
mdb.models['Model-1'].Pressure(name='Aero_Pressure', createStepName='RBM', 
	region=region1, distributionType=UNIFORM, field='', magnitude=13750.0, 
	amplitude=UNSET)
region2 = a.instances['Elastomer'].surfaces['Preslocface']
mdb.models['Model-1'].Pressure(name='Aero_Pressure2', createStepName='RBM', 
	region=region2, distributionType=UNIFORM, field='', magnitude=13750.0, amplitude=UNSET)		

#####################################
### Creation/Execution of the Job ###
#####################################
print 'Creating/Running Job'

ModelName='Model-1'

mdb.Job(name=ModelName, model=ModelName, description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, 
    userSubroutine='', 
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
    numGPUs=0)
job=mdb.jobs[ModelName]
# delete lock file, which for some reason tends to hang around, if it exists
if os.access('%s.lck'%ModelName,os.F_OK):
	os.remove('%s.lck'%ModelName)	
# Run the job, then process the results.
if test==1:
	stahp
job.submit()
job.waitForCompletion()
print 'Completed job'
[initial,final]=getResults2(ModelName)
	
# Code for writing results
# DataFile = open('PostData.txt','a')
# DataFile.write('%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f\n' % (IRTT,ORTT,T1,T2,T3,G,eval('HNG%s'%rowname),
# phiL2,phiR2,tipDisp[2],tipDisp[12],tipDisp[3],tipDisp[13],tipDisp[0],tipDisp[10],tipDisp[1],tipDisp[11],tipDisp[8],tipDisp[18],tipDisp[9],tipDisp[19],tipDisp[6],tipDisp[16],tipDisp[7],tipDisp[17],
# tipDisp[4],tipDisp[14],tipDisp[5],tipDisp[15],tipDisp[20],TT1_THE,TT2_THE,TT3_THE,TT1_EE,TT2_EE,TT3_EE,ELAST_EE,TT1_S,TT2_S,TT3_S,ELAST_S, ))
# DataFile.close()
# DataFile = open('PostDataInitial.txt','a')
# DataFile.write('%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f\n' % (IRTT,ORTT,T1,T2,T3,G,eval('HNG%s'%rowname),phiL2,phiR2,initial[2],initial[12],initial[3],initial[13],initial[0],initial[10],initial[1],initial[11],initial[8],initial[18],initial[9],initial[19],initial[6],initial[16],initial[7],initial[17],initial[4],initial[14],initial[5],initial[15],tipDisp[20],TT1_THE,TT2_THE,TT3_THE,TT1_EE,TT2_EE,TT3_EE,ELAST_EE,TT1_S,TT2_S,TT3_S,ELAST_S, ))
# DataFile.close()
# DataFile = open('PostDataFinal.txt','a')
# DataFile.write('%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,' % (IRTT,ORTT,T1,T2,T3,G,eval('HNG%s'%rowname),phiL2,phiR2,final[2],final[12],final[3],final[13],final[0],final[10],final[1],final[11],final[8],final[18],final[9],final[19],final[6],final[16],final[7],final[17],final[4],final[14],final[5],final[15],tipDisp[20],TT1_THE,TT2_THE,TT3_THE,TT1_EE,TT2_EE,TT3_EE,ELAST_EE,TT1_S,TT2_S,TT3_S,ELAST_S, ))
# DataFile.close()
# print 'DONE!!'
