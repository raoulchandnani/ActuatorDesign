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
from math import atan, sin, cos, tan,atan2,pi
import numpy as np
from Post_P_Script import getResults,getResults2,find_strain,find_stress
##########################

# Variables for Plates & Assembly
H=0.55; # Horizontal plate length 
V=0.055; # Vertical plate length
G=0.025/4;   # Gap length
HNG=-0.030/4; # Hinge height
D=0.015/4; #Hinge distance from edge
elast_thk=-0.005*6
Total = 0.28
Pad=0.01

##Angles
phiL2=16.3; 
phiR2=16.3;
phiL=3.15;
phiR=3.15;
phiM=0;

# Variables for Elastomer
with open("Aero_Data_I.txt", "r") as a:
    Datalow= [[float(num) for num in line.split(' ')] for line in a]
Datalow=np.array(Datalow) 
start=83
end=86
Elastfull=0

# Variables for TT
Ti=280 	# Initial temperature of all TTs
T1=520	# Final temperature of all TT1
T2=520	# Final temperature of all TT2
ORTT=0.008
IRTT=0.004
NumLayer = 5 #Number of layers in TT for Max transformation strain variation

# Variables for Plates & Assembly

Hnet=(Total)/(2*cos(phiL2*pi/180)+2*cos(phiL*pi/180)+1);
H=Hnet-G;
H1=H
H2=H
H3=H
H4=H
H5=H
seedsize=G/4
test=1
##########################


### Write data file column headings

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)


###Sketch Geometry and Create Parts
print 'Sketching/Creating the part'
Mdb()

#### TOrque tube ###
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
p.BaseSolidExtrude(sketch=s, depth=0.5*(V-2*D))
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
    (0.5*V-D))), point3=v1.findAt(coordinates=(ORTT, 0.0, -(0.5*V-D))), cells=pickedCells, 
    point1=p.InterestingPoint(edge=e.findAt(coordinates=(0,-ORTT, (0.5*V-D))), 
    rule=MIDDLE))
p = mdb.models['Model-1'].parts['Torque_Tube']
c = p.cells
pickedCells = c.findAt(((0, 0.5*(IRTT+ORTT), 0), ), ((0, 
     -0.5*(IRTT+ORTT), 0), ))
v2, e1, d2 = p.vertices, p.edges, p.datums
p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=p.InterestingPoint(
    edge=e1.findAt(coordinates=(0, ORTT, (0.5*V-D))), rule=MIDDLE), 
    point2=p.InterestingPoint(edge=e1.findAt(coordinates=(0, -ORTT, 
    (0.5*V-D))), rule=MIDDLE), point3=p.InterestingPoint(edge=e1.findAt(
    coordinates=(0, ORTT, -(0.5*V-D))), rule=MIDDLE))
p = mdb.models['Model-1'].parts['Torque_Tube']
DP=p.DatumPointByCoordinate(coords=(0, 0, (0.5*V-D)))
DP_TT1=DP.id
p = mdb.models['Model-1'].parts['Torque_Tube']
DP=p.DatumPointByCoordinate(coords=(0, 0, -(0.5*V-D)))
DP_TT2=DP.id

## Create sets##
p = mdb.models['Model-1'].parts['Torque_Tube']
e = p.edges
edges = e.findAt(((0.0, ORTT, 0), ), ((0.0, -ORTT, 0), ), ((-ORTT, 0.0, 
    0), ), ((ORTT, 0.0, 0), ))
p.Set(edges=edges, name='Longitudinal_edges')

p = mdb.models['Model-1'].parts['Torque_Tube']
e = p.edges
edges = e.findAt(((0.0, 0.5*(IRTT+ORTT), -(0.5*V-D)), ), ((0.0, 0.5*(IRTT+ORTT), (0.5*V-D)), ), ((0.0, 
    -0.5*(IRTT+ORTT), (0.5*V-D)), ), ((0.0, -0.5*(IRTT+ORTT), -(0.5*V-D)), ), ((0.5*(IRTT+ORTT), 0.0, (0.5*V-D)), ), ((
    -0.5*(IRTT+ORTT), 0.0, -(0.5*V-D)), ), ((0.5*(IRTT+ORTT), 0.0, -(0.5*V-D)), ), ((-0.5*(IRTT+ORTT), 0.0, (0.5*V-D)), ))
p.Set(edges=edges, name='Radial_edges')

p = mdb.models['Model-1'].parts['Torque_Tube']
e = p.edges
edges = e.findAt(((ORTT*cos(pi/4), ORTT*sin(pi/4), -(0.5*V-D)), ), ((-ORTT*cos(pi/4), ORTT*sin(pi/4), 
    (0.5*V-D)), ), ((ORTT*cos(pi/4), ORTT*sin(pi/4), (0.5*V-D)), ), ((-ORTT*cos(pi/4), -ORTT*sin(pi/4), -(0.5*V-D)), 
    ), ((ORTT*cos(pi/4), -ORTT*sin(pi/4), (0.5*V-D)), ), ((-ORTT*cos(pi/4), -ORTT*sin(pi/4), (0.5*V-D)), ), (
    (-ORTT*cos(pi/4), ORTT*sin(pi/4), -(0.5*V-D)), ), ((ORTT*cos(pi/4), -ORTT*sin(pi/4), -(0.5*V-D)), ))
p.Set(edges=edges, name='Circumferential_edges')
f = p.faces
faces1 = f.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 
    -(0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ))
p.Set(faces=faces1, name='Rear_Faces')

faces1 = f.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), (0.5*V-D)), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 
    (0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), (0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), (0.5*V-D)), ))
p.Set(faces=faces1, name='Front_Faces')

## Partition the geometry for Max transformation strain variation##

p = mdb.models['Model-1'].parts['Torque_Tube']
f1, e1, d2 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f1.findAt(coordinates=(0.5*(IRTT+ORTT)*cos(pi/4), 
    0.5*(IRTT+ORTT)*cos(pi/4), (0.5*V-D))), sketchUpEdge=e1.findAt(coordinates=(ORTT*cos(pi/4), 
    ORTT*cos(pi/4), (0.5*V-D))), sketchPlaneSide=SIDE1, origin=(0, 0, 
    (0.5*V-D)))
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=0.052, gridSpacing=0.001, transform=t)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Torque_Tube']
p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
for n in range(NumLayer-1):
	s1.CircleByCenterPerimeter(center=(0, 0), point1=(IRTT + (n+1)*(ORTT-IRTT)/(NumLayer), 0))
p = mdb.models['Model-1'].parts['Torque_Tube']
c = p.cells
f, e, d1 = p.faces, p.edges, p.datums
p.PartitionCellBySketch(sketchPlane=f.findAt(coordinates=(0.5*(IRTT+ORTT)*cos(pi/4), 
    0.5*(IRTT+ORTT)*cos(pi/4), (0.5*V-D))), sketchUpEdge=e1.findAt(coordinates=(ORTT*cos(pi/4), 
    ORTT*cos(pi/4), (0.5*V-D))), cells=c[:], sketch=s1)
for n in range(NumLayer-1):
	pickedEdges = [e.findAt(coordinates=((IRTT + (n+1)*(ORTT-IRTT)/(NumLayer))*cos(pi/4),(IRTT + (n+1)*(ORTT-IRTT)/(NumLayer))*cos(pi/4), (0.5*V-D))), 
    e.findAt(coordinates=((IRTT + (n+1)*(ORTT-IRTT)/(NumLayer))*cos(pi/4),-(IRTT + (n+1)*(ORTT-IRTT)/(NumLayer))*cos(pi/4), (0.5*V-D))), 
	e.findAt(coordinates=(-(IRTT + (n+1)*(ORTT-IRTT)/(NumLayer))*cos(pi/4),(IRTT + (n+1)*(ORTT-IRTT)/(NumLayer))*cos(pi/4), (0.5*V-D))), 
	e.findAt(coordinates=(-(IRTT + (n+1)*(ORTT-IRTT)/(NumLayer))*cos(pi/4),-(IRTT + (n+1)*(ORTT-IRTT)/(NumLayer))*cos(pi/4), (0.5*V-D)))]
	pickedEdges = tuple(pickedEdges)
	p.PartitionCellByExtrudeEdge(line=e.findAt(coordinates=(ORTT, 0.0, 0.25*(V-D))), 
		cells=c[:], edges=pickedEdges, sense=FORWARD)
s1.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']

#### Plate ###
for plate_num in [1,2,3,4,5]:
	s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
		sheetSize=200.0)
	g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
	s.setPrimaryObject(option=STANDALONE)
	s.rectangle(point1=((-0.5*eval('H%i'%plate_num)), (-0.5*V)), point2=(0.5*eval('H%i'%plate_num),0.5*V))
	p = mdb.models['Model-1'].Part(name='Plate%i'%plate_num, dimensionality=THREE_D, 
		type=DEFORMABLE_BODY)
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	p.BaseShell(sketch=s)
	s.unsetPrimaryObject()
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	del mdb.models['Model-1'].sketches['__profile__']
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	s = p.faces
	side1Faces = s.findAt(((0, 0, 0.0), ))
	p.Surface(side2Faces=side1Faces, name='Surf-Bottom')
	
	#Defining the face partitions
	print 'Partitioning part'
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	f1, e2, d2 = p.faces, p.edges, p.datums
	t = p.MakeSketchTransform(sketchPlane=f1.findAt(coordinates=(0, 
		0, 0.0), normal=(0.0, 0.0, 1.0)), sketchUpEdge=e2.findAt(
		coordinates=(0.5*eval('H%i'%plate_num),0, 0.0)), sketchPlaneSide=SIDE1, origin=(0, 0,0.0))
	s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
		sheetSize=459.88, gridSpacing=11.49, transform=t)
	g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
	s.setPrimaryObject(option=SUPERIMPOSE)
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

	s.rectangle(point1=(-(0.5*eval('H%i'%plate_num)-D), -(0.5*V-D)), point2=((0.5*eval('H%i'%plate_num)-D), (0.5*V-D)))
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	f = p.faces
	pickedFaces = f.findAt(((0,0, 0.0), ))
	e, d1 = p.edges, p.datums
	p.PartitionFaceBySketch(sketchUpEdge=e.findAt(coordinates=(0.5*eval('H%i'%plate_num),0, 0.0)), 
		faces=pickedFaces, sketch=s)
	s.unsetPrimaryObject()
	del mdb.models['Model-1'].sketches['__profile__']

	###Datum Points
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	DP=p.DatumPointByCoordinate(coords=((0.5*eval('H%i'%plate_num)-D), 0.5*V-D, 0))
	DP2=DP.id
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	DP=p.DatumPointByCoordinate(coords=((0.5*eval('H%i'%plate_num)-D), -0.5*V+D, 0))
	DP1=DP.id
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	DP=p.DatumPointByCoordinate(coords=(-(0.5*eval('H%i'%plate_num)-D), 0.5*V-D, 0))
	DP3=DP.id
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	DP=p.DatumPointByCoordinate(coords=(-(0.5*eval('H%i'%plate_num)-D), -0.5*V+D, 0))
	DP4=DP.id
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	DP=p.DatumPointByCoordinate(coords=((0.5*eval('H%i'%plate_num)), 0.5*V, 0))
	DP6=DP.id
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	DP=p.DatumPointByCoordinate(coords=((0.5*eval('H%i'%plate_num)), -0.5*V, 0))
	DP5=DP.id
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	DP=p.DatumPointByCoordinate(coords=(-(0.5*eval('H%i'%plate_num)), 0.5*V, 0))
	DP7=DP.id
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	DP=p.DatumPointByCoordinate(coords=(-(0.5*eval('H%i'%plate_num)), -0.5*V, 0))
	DP8=DP.id
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	DP=p.DatumPointByCoordinate(coords=((0.5*(eval('H%i'%plate_num)+G)), 0.5*V-D, HNG))
	DP10=DP.id
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	DP=p.DatumPointByCoordinate(coords=((0.5*(eval('H%i'%plate_num)+G)), -0.5*V+D, HNG))
	DP9=DP.id
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	DP=p.DatumPointByCoordinate(coords=(-(0.5*(eval('H%i'%plate_num)+G)), 0.5*V-D, HNG))
	DP11=DP.id
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	DP=p.DatumPointByCoordinate(coords=(-(0.5*(eval('H%i'%plate_num)+G)), -0.5*V+D, HNG))
	DP12=DP.id
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	DP=p.DatumPointByCoordinate(coords=(0, 0, 0))
	DP13=DP.id


#### Assembly ###
print 'Placing Parts in Space'
#Create Instances here

a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Plate1']
a.Instance(name='Left_Plate_2', part=p, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Plate2']
a.Instance(name='Left_Plate', part=p, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Plate3']
a.Instance(name='Middle_Plate', part=p, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Plate4']
a.Instance(name='Right_Plate', part=p, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Plate5']
a.Instance(name='Right_Plate_2', part=p, dependent=ON)

#Position plates
Vct=[0,0,0]
DVct=[-0.5*Total,-0.5*V+D,0]
for loc in range(3):
	Vct[loc]= DVct[loc] - mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn[loc]
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('Left_Plate_2', ), vector=(Vct[0],Vct[1], Vct[2]))
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Left_Plate_2', ), axisPoint=(DVct[0],DVct[1], DVct[2]), 
	axisDirection=(0.0, -(2*V), 0.0), angle=(phiL2))	
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP9].pointOn
for loc in range(3):
	Vct[loc]= DVct[loc] - mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP12].pointOn[loc]
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('Left_Plate', ), vector=(Vct[0],Vct[1], Vct[2]))
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Left_Plate', ), axisPoint=(DVct[0],DVct[1], DVct[2]), 
	axisDirection=(0.0, -(2*V), 0.0), angle=(phiL))		
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP9].pointOn
for loc in range(3):
	Vct[loc]= DVct[loc] - mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP12].pointOn[loc]
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('Middle_Plate', ), vector=(Vct[0],Vct[1], Vct[2]))
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Middle_Plate', ), axisPoint=(DVct[0],DVct[1], DVct[2]), 
	axisDirection=(0.0, -(2*V), 0.0), angle=(phiM))			

DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP9].pointOn
for loc in range(3):
	Vct[loc]= DVct[loc] - mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP12].pointOn[loc]
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('Right_Plate', ), vector=(Vct[0],Vct[1], Vct[2]))
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Right_Plate', ), axisPoint=(DVct[0],DVct[1], DVct[2]), 
	axisDirection=(0.0, -(2*V), 0.0), angle=(-phiR))	

DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP9].pointOn
for loc in range(3):
	Vct[loc]= DVct[loc] - mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP12].pointOn[loc]
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('Right_Plate_2', ), vector=(Vct[0],Vct[1], Vct[2]))
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Right_Plate_2', ), axisPoint=(DVct[0],DVct[1], DVct[2]), 
	axisDirection=(0.0, (2*V), 0.0), angle=(phiR2))		

p = mdb.models['Model-1'].parts['Torque_Tube']
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
a.Instance(name='Torque_Tube_1', part=p, dependent=ON)

#Create Reference Points
print 'Defining Reference Points'

## Reference points##

DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly													# RPs for Middle to Right 1
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#8,9
RP1=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP2=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly													# RPs for Middle to Right 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))													#10,11
RP3=RP.id															
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP4=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP6=RP.id
a = mdb.models['Model-1'].rootAssembly													# RPs for Middle to Left 1
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#12,13
RP5=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP7=RP.id
a = mdb.models['Model-1'].rootAssembly														# RPs for Middle to Left 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#14,15
RP8=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly															# RPs for Right to Right
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#16,17
RP9=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP10=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly															# RPs for Right to Right 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#18,19
RP11=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP12=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly														# RPs for Left to Left 
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#RP13,21
RP13=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP14=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly														# RPs for Left to Left  2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#22,23
RP15=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP16=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly														# RPs for Rigid Body Constraints								#24,25,26
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP17=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP18=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP19=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly																	# RPs for Hinge Support at endswith								#27,28,29,30
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP20=RP.id
a = mdb.models['Model-1'].rootAssembly																	# RPs for Hinge Support at endswith								#27,28,29,30
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP26=RP.id
r1 = a.referencePoints
refPoints1=(r1[RP20], )
a.Set(referencePoints=refPoints1, name='RP_Z3+')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP21=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP27=RP.id
r1 = a.referencePoints
refPoints1=(r1[RP21], )
a.Set(referencePoints=refPoints1, name='RP_Z3-')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP22=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP28=RP.id
r1 = a.referencePoints
refPoints1=(r1[RP22], )
a.Set(referencePoints=refPoints1, name='RP_Z1-')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP23=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP29=RP.id
r1 = a.referencePoints
refPoints1=(r1[RP23], )
a.Set(referencePoints=refPoints1, name='RP_Z1+')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP24=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP25=RP.id

## Create Elastomer
Presloc=[0,0,0]
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
	sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
p = mdb.models['Model-1'].Part(name='Elastomer', dimensionality=THREE_D, 
	type=DEFORMABLE_BODY)
d2 = p.datums
Pathloc=(0,0,0,0,0,0)
Setloc=np.zeros((5,3))
SetlocL=np.zeros((5,3))
SetlocT=np.zeros((5,3))
Curveloc=np.zeros((1,3))
for xloc in range(start,end):
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
	Dataxloc=[[0.000000 for x in range(7)] for y in range(len(Datalow))]
	Dataxloc=np.array(Dataxloc)
	for i in range(len(Datalow)):
		if abs(Datalow[i,0]-XLOC)<1E-06:
			Dataxloc[count]=Datalow[i]
			count=count+1
	Dataxloc=Dataxloc[~np.all(Dataxloc == 0, axis=1)]
	Pad_thk=-(min(Datalow[:,2])-min(Dataxloc[:,2]))-0.04#-(3.15-0.0315*xloc)*tan(6.841917290255273E-02)+0.04
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
	Pad_thk=DVct2[2]+max(Datalow[:,2])+(0.0315*100-0.0315*xloc)*tan(6.841917290255273E-02)-0.215
	j=0
	if xloc==start:	
		for i in range(len(Dataxloc)):
			Curveloc=np.vstack((Curveloc,[Dataxloc[i,1],-(Dataxloc[i,2]),0.0315*xloc]))
		Curveloc=Curveloc[~np.all(Curveloc == 0, axis=1)]	
	if xloc==(start+end)/2:
		for i in range(len(Dataxloc)-1):
			Presloc=np.vstack((Presloc,[0.7*(Dataxloc[i,1])+0.3*(Dataxloc[i+1,1]),0.7*(-(Dataxloc[i,2]))-0.3*(Dataxloc[i+1,2]),0.0315*xloc]))
			Presloc=np.vstack((Presloc,[0.3*(Dataxloc[i,1])+0.7*(Dataxloc[i+1,1]),0.3*(-(Dataxloc[i,2]))-0.7*(Dataxloc[i+1,2]),0.0315*xloc]))
			s.Line(point1=(Dataxloc[i,1], -(Dataxloc[i,2])), point2=(Dataxloc[i+1,1], -(Dataxloc[i+1,2])))				
	else:		
		for i in range(len(Dataxloc)-1):
			s.Spline(points=((Dataxloc[i,1], -(Dataxloc[i,2])), (Dataxloc[i+1,1], -(Dataxloc[i+1,2]))))				
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP8].pointOn
	s.Line(point1=(DVct[0]-Pad, DVct2[2]-Pad_thk), point2=(Dataxloc[0,1], -(Dataxloc[0,2])))
#	s.Line(point1=(DVct[0]-Pad, DVct2[2]-Pad_thk), point2=(DVct[0]-Pad, -(Dataxloc[0,2])))
	s.Line(point1=(DVct[0]-Pad, DVct2[2]-Pad_thk), point2=(DVct2[0], DVct2[2]-Pad_thk))
	if xloc==start or xloc==end-1:
		Path=(11.85+xloc*0.0315,DVct[0]-Pad, DVct2[2]-Pad_thk)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP5].pointOn
	s.Line(point1=(DVct2[0], DVct2[2]-Pad_thk), point2=(DVct[0], DVct[2]-Pad_thk))
	if xloc==start+1:
		SetlocT[0]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])-Pad_thk)
	if xloc==end-2:
		SetlocL[0]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])-Pad_thk)		
	if xloc==(start+end)/2:
		Setloc[0]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])-Pad_thk)
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP8].pointOn
	s.Line(point1=(DVct[0], DVct[2]-Pad_thk), point2=(DVct2[0], DVct2[2]-Pad_thk))
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP5].pointOn
	s.Line(point1=(DVct2[0], DVct2[2]-Pad_thk), point2=(DVct[0], DVct[2]-Pad_thk))
	if xloc==start+1:
		SetlocT[1]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])-Pad_thk)
	if xloc==end-2:
		SetlocL[1]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])-Pad_thk)		
	if xloc==(start+end)/2:
		Setloc[1]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])-Pad_thk)
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP8].pointOn
	s.Line(point1=(DVct[0], DVct[2]-Pad_thk), point2=(DVct2[0], DVct2[2]-Pad_thk))
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP5].pointOn
	s.Line(point1=(DVct2[0], DVct2[2]-Pad_thk), point2=(DVct[0], DVct[2]-Pad_thk))
	if xloc==start+1:
		SetlocT[2]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])-Pad_thk)
	if xloc==end-2:
		SetlocL[2]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])-Pad_thk)		
	if xloc==(start+end)/2:
		Setloc[2]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])-Pad_thk)			
	if xloc==start:
		Axis_start=(DVct2[0], DVct2[2]-Pad_thk,0.0315*xloc)
		Axis_end=(DVct[0], DVct[2]-Pad_thk,0.0315*xloc)
		Axis_1_avg=(0.5*(DVct[0]+DVct2[0]), 0.5*(DVct[2]+DVct2[2])-Pad_thk,0.0315*xloc)
	if xloc==end-1:
		Axis_2_avg=(0.5*(DVct[0]+DVct2[0]), 0.5*(DVct[2]+DVct2[2])-Pad_thk,0.0315*xloc)
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP8].pointOn
	s.Line(point1=(DVct[0], DVct[2]-Pad_thk), point2=(DVct2[0], DVct2[2]-Pad_thk))
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP5].pointOn
	s.Line(point1=(DVct2[0], DVct2[2]-Pad_thk), point2=(DVct[0], DVct[2]-Pad_thk))
	if xloc==start+1:
		SetlocT[3]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])-Pad_thk)
	if xloc==end-2:
		SetlocL[3]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])-Pad_thk)		
	if xloc==(start+end)/2:
		Setloc[3]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])-Pad_thk)
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP8].pointOn
	s.Line(point1=(DVct[0], DVct[2]-Pad_thk), point2=(DVct2[0], DVct2[2]-Pad_thk))
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP5].pointOn
	s.Line(point1=(DVct2[0], DVct2[2]-Pad_thk), point2=(DVct[0], DVct[2]-Pad_thk))
	if xloc==start+1:
		SetlocT[4]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])-Pad_thk)
	if xloc==end-2:
		SetlocL[4]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])-Pad_thk)		
	if xloc==(start+end)/2:
		Setloc[4]=(0.0315*xloc,0.5*(DVct[0]+DVct2[0]),0.5*(DVct[2]+DVct2[2])-Pad_thk)
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP9].pointOn
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP8].pointOn
	DVct3=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
	s.Line(point1=(DVct[0], DVct[2]-Pad_thk), point2=(-DVct3[0]+Pad, DVct2[2]-Pad_thk))
	s.Line(point1=(Dataxloc[-1,1], -(Dataxloc[-1,2])), point2=(-DVct3[0]+Pad, DVct2[2]-Pad_thk))
	if xloc==(start+end)/2:
		DVct1=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP13].pointOn
		DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP13].pointOn
		DVct3=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
		DVct4=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP13].pointOn
		DVct5=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP13].pointOn
		A=[[DVct1[0],DVct1[2]-Pad_thk,0.0315*xloc],[DVct2[0],DVct2[2]-Pad_thk,0.0315*xloc],[DVct3[0],DVct3[2]-Pad_thk,0.0315*xloc],[DVct4[0],DVct4[2]-Pad_thk,0.0315*xloc],[DVct5[0],DVct5[2]-Pad_thk,0.0315*xloc]]
	if xloc==start or xloc==end-1:
		Path=np.hstack((Path,(11.85+xloc*0.0315,-DVct3[0]+Pad, DVct2[2]-Pad_thk)))
		Pathloc=np.vstack((Pathloc,Path))
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
	del d2[DPE]

if Elastfull==1:
	for i in range(29):
		Pathh=[0,0,0]
		for j in range(i+1,len(Datalow),31):
			Pathh=np.vstack((Pathh,Datalow[j,0:3]))
		Pathh=Pathh[~np.all(Pathh == 0, axis=1)]
		p.Set(edges=e.getByBoundingBox(Pathh[0,1],0,0,Pathh[0,1],100,4), name='PathWire%i'%i)
		DP=p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=Pathh[0,1])
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
		# for k in range(len(Pathh)-1):
			# s.Spline(points=((-0.0315*(start+k), -Pathh[k,2]), (-0.0315*(start+k+1), -Pathh[k+1,2])))
		s.Spline(points=[(-0.0315*(start+k), -Pathh[k,2]) for k in range(len(Pathh))])
		W=p.Wire(sketchPlane=d2[DPP], sketchUpEdge=d2[DPE], sketchPlaneSide=SIDE1, 
			sketchOrientation=RIGHT, sketch=s)
		s.unsetPrimaryObject()
		p.Set(edges=e.getByBoundingBox(Pathh[0,1],0,0,Pathh[0,1],100,4), name='PathWire%i'%i)

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
s.unsetPrimaryObject()
p.Set(edges=[e.findAt(((Pathloc[1,1], Pathloc[1,2]+0.0315*(loca-0.5)*tan(6.841917290255273E-02), Pathloc[1,0]-11.85+0.0315*(loca-0.5)), ),  ) for loca in range(1,end-start)], name='PathWireL')

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
s.unsetPrimaryObject()
p.Set(edges=[e.findAt(((Pathloc[1,4], Pathloc[1,5]+0.0315*(loca-0.5)*tan(6.841917290255273E-02), Pathloc[1,3]-11.85+0.0315*(loca-0.5)), ),  ) for loca in range(1,end-start)], name='PathWireR')

if Elastfull==1:
	p.SolidLoft(loftsections=[p.sets['E%i'%xloc].edges[:] for xloc in range(start,end)],paths=[p.sets['PathWire%i'%locc].edges[:] for locc in range(29)],globalSmoothing=ON)
else:
	p.SolidLoft(loftsections=[p.sets['E%i'%xloc].edges[:] for xloc in range(start,end)],paths=(p.sets['PathWireL'].edges[:],p.sets['PathWireR'].edges[:],),globalSmoothing=ON)

v=p.vertices
p.Set(vertices=[v.findAt(((Curveloc[l,0],Curveloc[l,1],Curveloc[l,2]), ),  ) for l in range(len(Curveloc))], name='CurveTip')
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Elastomer']
a.Instance(name='Elastomer', part=p, dependent=ON)
v=a.instances['Elastomer'].vertices
a.Set(vertices=[v.findAt(((Curveloc[l,0],Curveloc[l,1],Curveloc[l,2]), ),  ) for l in range(len(Curveloc))], name='TIPNODE')
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Elastomer', ), axisPoint=Axis_start, 
	axisDirection=(Axis_end[0]-Axis_start[0],0,0), angle=90+6.841917290255273E-02*180/pi)	
for loc in range(3):
	Vct[loc]= mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP7].pointOn[loc]-Axis_start[loc]
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('Elastomer', ), vector=(Vct[0],Vct[1]+0.5*((0.0315*(end-1)-0.0315*start)/cos(6.841917290255273E-02)-V),Vct[2]))

#Defining the face partitions of elastomer
print 'Partitioning part'
p = mdb.models['Model-1'].parts['Elastomer']
f, e, d = p.faces, p.edges, p.datums
Vct2=[0,0,0]
for loc in range(3):
	Vct[loc]= 0.5*(Axis_end[loc]+Axis_start[loc])
	Vct2[loc]=0.5*(Axis_1_avg[loc]+Axis_2_avg[loc])
t = p.MakeSketchTransform(sketchPlane=f.findAt(coordinates=A[2]), sketchUpEdge=e.findAt(coordinates=Vct), 
	sketchPlaneSide=SIDE2, origin=Vct2)
#	0.0, 0.59661, 2.142
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=6.01, 
	gridSpacing=0.15, transform=t)
g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Elastomer']
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
s.Line(point1=(0.5*V, 0.15), point2=(0.5*V, -0.15))
s.Line(point1=(-0.5*V, 0.15), point2=(-0.5*V, -0.15))
p = mdb.models['Model-1'].parts['Elastomer']
f = p.faces
pickedFaces = f.findAt((A[0], ), (A[1], ), (A[2], ), (A[3], ), (A[4], ))
f1, e1, d2 = p.faces, p.edges, p.datums
p.PartitionFaceBySketchThruAll(sketchPlane=f1.findAt(coordinates=A[2]), sketchUpEdge=e1.findAt(coordinates=Vct), faces=pickedFaces, sketchPlaneSide=SIDE2, sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
for rowname in ['']:
	p = mdb.models['Model-1'].parts['Elastomer']
	s = p.faces 
	side1Faces = s.findAt(((eval('Setloc%s'%rowname)[0][1], eval('Setloc%s'%rowname)[0][2], eval('Setloc%s'%rowname)[0][0]), ))
	p.Surface(side1Faces=side1Faces, name='Surf_L2%s'%rowname)
	side1Faces = s.findAt(((eval('Setloc%s'%rowname)[1][1], eval('Setloc%s'%rowname)[1][2], eval('Setloc%s'%rowname)[1][0]), ))
	p.Surface(side1Faces=side1Faces, name='Surf_L%s'%rowname)
	side1Faces = s.findAt(((eval('Setloc%s'%rowname)[2][1], eval('Setloc%s'%rowname)[2][2], eval('Setloc%s'%rowname)[2][0]), ))
	p.Surface(side1Faces=side1Faces, name='Surf_M%s'%rowname)
	side1Faces = s.findAt(((eval('Setloc%s'%rowname)[3][1], eval('Setloc%s'%rowname)[3][2], eval('Setloc%s'%rowname)[3][0]), ))
	p.Surface(side1Faces=side1Faces, name='Surf_R%s'%rowname)
	side1Faces = s.findAt(((eval('Setloc%s'%rowname)[4][1], eval('Setloc%s'%rowname)[4][2], eval('Setloc%s'%rowname)[4][0]), ))
	p.Surface(side1Faces=side1Faces, name='Surf_R2%s'%rowname)	



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
p = mdb.models['Model-1'].parts['Torque_Tube']
v1, e = p.vertices, p.edges
DP=p.DatumCsysByThreePoints(point1=v1.findAt(coordinates=(ORTT, 0.0,-(0.5*V-D))), 
    name='Datum csys-1', coordSysType=CYLINDRICAL, origin=(0,0,-(0.5*V-D)), 
    point2=(ORTT,ORTT,-(0.5*V-D)))
DPCYS=DP.id
DP=p.DatumCsysByThreePoints(point1=v1.findAt(coordinates=(ORTT, 0.0,(0.5*V-D))), 
    name='Datum csys-2', coordSysType=CYLINDRICAL, origin=(0,0,(0.5*V-D)), 
    point2=(ORTT,ORTT,(0.5*V-D)))
DPCYS2=DP.id
for n in range(NumLayer):
	nm = 'SMA' + str(n)
	mdb.models['Model-1'].Material(name=nm)
	mdb.models['Model-1'].materials[nm].Density(table=((6500000.0, ), ))
	mdb.models['Model-1'].materials[nm].Depvar(n=31)
	hsat = (IRTT + (n+0.5)*(ORTT-IRTT)/NumLayer)*0.048202422/(IRTT + (NumLayer-0.5)*(ORTT-IRTT)/NumLayer)
	mdb.models['Model-1'].materials[nm].UserMaterial(mechanicalConstants=(6500.00
	,6.62E+10,2.56E+10,0.33,300,0.00000000,0.00000000,1.50E+08,9017722.07,8342551.97
	,470.00,470.00,285.000,310.000,315.000,330.000,hsat,1.29E-08,5.00E+07
	,0.00E+00,0.00E+00,0.00E+00,2.50E+07,5.00E+07,0.9500,0.9500,0.9500,0.9500
	,1.00E-06,1.00E-09,0.9999,6,0,0,1,0,0,0))
																				
	mdb.models['Model-1'].materials[nm].Expansion(table=((0.0, ), ))
	mdb.models['Model-1'].materials[nm].Conductivity(table=((10.0, ), ))
																			 
	mdb.models['Model-1'].materials[nm].SpecificHeat(table=((0.32, ), ))
	nms ='SMA_Section' + str(n)
	mdb.models['Model-1'].HomogeneousSolidSection(name=nms, material=nm,                                                                
		thickness=None)
	c = p.cells
	cells = c.findAt(((-(IRTT + (n+0.5)*(ORTT-IRTT)/NumLayer)*cos(pi/4), (IRTT + (n+0.5)*(ORTT-IRTT)/NumLayer)*cos(pi/4), 0), ), ((-(IRTT + (n+0.5)*(ORTT-IRTT)/NumLayer)*cos(pi/4), 
		-(IRTT + (n+0.5)*(ORTT-IRTT)/NumLayer)*cos(pi/4), 0), ), (((IRTT + (n+0.5)*(ORTT-IRTT)/NumLayer)*cos(pi/4), -(IRTT + (n+0.5)*(ORTT-IRTT)/NumLayer)*cos(pi/4), 0), ), ((
		(IRTT + (n+0.5)*(ORTT-IRTT)/NumLayer)*cos(pi/4), (IRTT + (n+0.5)*(ORTT-IRTT)/NumLayer)*cos(pi/4),0), ))
	region = regionToolset.Region(cells=cells)
	orientation = mdb.models['Model-1'].parts['Torque_Tube'].datums[DPCYS2]
	mdb.models['Model-1'].parts['Torque_Tube'].MaterialOrientation(region=region, 
		orientationType=SYSTEM, axis=AXIS_3, localCsys=orientation, fieldName='', 
		additionalRotationType=ROTATION_NONE, angle=0.0, 
		additionalRotationField='', stackDirection=STACK_3)	
	region = p.Set(cells=cells, name='TT' + str(n))	
	p.SectionAssignment(region=region, sectionName=nms, offset=0.0, 
		offsetType=MIDDLE_SURFACE, offsetField='', 
		thicknessAssignment=FROM_SECTION)
p.SetByBoolean(name='TT', sets = tuple([p.sets['TT' + str(n)] for n in range(NumLayer)]))

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

for plate_num in [1,2,3,4,5]:	
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	f = p.faces
	faces = f.getByBoundingBox(-H,-V,-G,H,V,G)
	region = p.Set(faces=faces, name='Set-1')
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
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

	
#Create Wires
print 'Defining Wires'

a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
v1 = a.instances['Middle_Plate'].vertices
v2 = a.instances['Right_Plate'].vertices
v3 = a.instances['Left_Plate'].vertices
v4 = a.instances['Right_Plate_2'].vertices
v5 = a.instances['Left_Plate_2'].vertices

DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP1], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP2], v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP3], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP4], v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP6], v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP8], v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP7], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP5], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[RP16]), ), mergeType=IMPRINT, meshable=OFF)			
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[RP14]), ), mergeType=IMPRINT, meshable=OFF)				
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[RP10]), ), mergeType=IMPRINT, meshable=OFF)				
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[RP12]), ), mergeType=IMPRINT, meshable=OFF)				
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP9], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)	
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP11], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)	
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP13], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)		
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP15], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)		
# DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP1].pointOn
# a = mdb.models['Model-1'].rootAssembly
# r11 = a.referencePoints
# a.WirePolyLine(points=((r1[RP21], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)	
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP20], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
# DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP3].pointOn
# a = mdb.models['Model-1'].rootAssembly
# r11 = a.referencePoints
# a.WirePolyLine(points=((r1[RP22], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP23], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)

a = mdb.models['Model-1'].rootAssembly
e1 = a.edges 
edges1 = e1.getByBoundingBox(-4*H,-4*V,-2*H,4*H,4*V,2*H)
a.Set(edges=edges1, name='Wire-Beams')

a = mdb.models['Model-1'].rootAssembly
a.suppressFeatures(('Wire-1', 'Wire-2', 'Wire-3', 'Wire-4', 'Wire-5', 'Wire-6', 
	'Wire-7', 'Wire-8', 'Wire-9', 'Wire-10','Wire-11','Wire-12','Wire-13','Wire-14','Wire-15','Wire-16','Wire-17','Wire-18','Wire-19','Wire-20', ))

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP1], r11[RP2]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP10].pointOn
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-1')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP3], r11[RP4]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-2')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP5], r11[RP6]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-3')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP7], r11[RP8]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-4')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP9], r11[RP10]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-5')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP11], r11[RP12]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-6')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP13], r11[RP14]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-7')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP15], r11[RP16]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-8')

DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP21], r11[RP27]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-9')

DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP20], r11[RP26]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-10')

DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP22], r11[RP28]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-11')

DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP23], r11[RP29]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-12')



a=mdb.models['Model-1'].rootAssembly
a.SetByBoolean(name='MRH', sets=(a.sets['Wire-1'], a.sets['Wire-2'], ))
a=mdb.models['Model-1'].rootAssembly
a.SetByBoolean(name='MLH', sets=(a.sets['Wire-3'], a.sets['Wire-4'], ))
a=mdb.models['Model-1'].rootAssembly
a.SetByBoolean(name='RRH', sets=(a.sets['Wire-5'], a.sets['Wire-6'], ))
a=mdb.models['Model-1'].rootAssembly
a.SetByBoolean(name='LLH', sets=(a.sets['Wire-7'], a.sets['Wire-8'], ))

a = mdb.models['Model-1'].rootAssembly
a.resumeFeatures(('Wire-1', 'Wire-2', 'Wire-3', 'Wire-4', 'Wire-5', 'Wire-6', 
	'Wire-7', 'Wire-8', 'Wire-9', 'Wire-10','Wire-11','Wire-12','Wire-13','Wire-14','Wire-15','Wire-16','Wire-17','Wire-18','Wire-19','Wire-20', ))

#Create and Assign Connector Sections
print 'Create Connector Section'	
mdb.models['Model-1'].ConnectorSection(name='HINGE', assembledType=HINGE)

print 'Create Connector Section'	
mdb.models['Model-1'].ConnectorSection(name='BEAM', assembledType=BEAM)

print 'Define Datum Coordinate Systems'
a = mdb.models['Model-1'].rootAssembly
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP10].pointOn
a.DatumCsysByThreePoints(name='CSYS_MR', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for Middle to Right 1
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_MR2', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for Middle to R2
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_ML', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for Middle to Left
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_ML2', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for Middle to L2
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_RR', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for  RHS Support
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP9].pointOn	
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_RR2', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for RHS Support 2
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_LL', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_LL2', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2

DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_RS', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_RS2', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_LS', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_LS2', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2

print 'Assign Connector Sections'

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-1']
datum1 = a.datums[a.features['CSYS_MR'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-2']
datum1 = a.datums[a.features['CSYS_MR2'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)
a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-3']
datum1 = a.datums[a.features['CSYS_ML'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-4']
datum1 = a.datums[a.features['CSYS_ML2'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-5']
datum1 = a.datums[a.features['CSYS_RR'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-6']
datum1 = a.datums[a.features['CSYS_RR2'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-7']
datum1 = a.datums[a.features['CSYS_LL'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-8']
datum1 = a.datums[a.features['CSYS_LL2'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-9']
datum1 = a.datums[a.features['CSYS_RS'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-10']
datum1 = a.datums[a.features['CSYS_RS2'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-11']
datum1 = a.datums[a.features['CSYS_LS'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-12']
datum1 = a.datums[a.features['CSYS_LS2'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-Beams']
csa = a.SectionAssignment(sectionName='BEAM', region=region)

#Mesh Parts
print 'Meshing the Part'
for plate_num in [1,2,3,4,5]:	
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	p.seedPart(size=0.02, deviationFactor=0.1, minSizeFactor=0.1)		#G/6
	p = mdb.models['Model-1'].parts['Plate%i'%plate_num]
	p.generateMesh()

elemType1 = mesh.ElemType(elemCode=C3D8H, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['Elastomer']
p.seedPart(size=0.01, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['Elastomer']
c = p.cells
cells = c[:]
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
    elemType3))
p = mdb.models['Model-1'].parts['Elastomer']
p.generateMesh()
## Mesh##
elemType1 = mesh.ElemType(elemCode=C3D8, elemLibrary=STANDARD, 
    kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, 
    hourglassControl=DEFAULT, distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['Torque_Tube']
c = p.cells
cells = c[:]
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
    elemType3))
# p = mdb.models['Model-1'].parts['Torque_Tube']
# e = p.edges
# pickedEdges = p.sets['Longitudinal_edges'].edges[:]
# p.seedEdgeByNumber(edges=pickedEdges, number=30, constraint=FIXED)
p = mdb.models['Model-1'].parts['Torque_Tube']
e = p.edges
pickedEdges = p.sets['Radial_edges'].edges[:]
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FIXED)

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
a.Instance(name='Torque_Tube_2', part=p, dependent=ON)

## Constraints##
a = mdb.models['Model-1'].rootAssembly
region1=a.sets['RP_Z1-']
a = mdb.models['Model-1'].rootAssembly
region2=a.instances['Torque_Tube_1'].sets['Rear_Faces']
mdb.models['Model-1'].Coupling(name='TT1_KC-', controlPoint=region1, 
    surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
    localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

a = mdb.models['Model-1'].rootAssembly
region1=a.sets['RP_Z1+']
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Torque_Tube_1'].faces
orientation = mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_1'].datums[DPCYS2]
region2=a.instances['Torque_Tube_1'].sets['Front_Faces']
mdb.models['Model-1'].Coupling(name='TT1_KC+', controlPoint=region1, 
    surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
    localCsys=orientation, u1=ON, u2=ON, u3=OFF, ur1=ON, ur2=ON, ur3=ON)

a = mdb.models['Model-1'].rootAssembly
region1=a.sets['RP_Z3-']
a = mdb.models['Model-1'].rootAssembly
orientation = mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_2'].datums[DPCYS]
region2=a.instances['Torque_Tube_2'].sets['Rear_Faces']
mdb.models['Model-1'].Coupling(name='TT2_KC-', controlPoint=region1, 
    surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
    localCsys=orientation, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
	
a = mdb.models['Model-1'].rootAssembly
region1=a.sets['RP_Z3+']
a = mdb.models['Model-1'].rootAssembly
orientation = mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_2'].datums[DPCYS2]
region2=a.instances['Torque_Tube_2'].sets['Front_Faces']
mdb.models['Model-1'].Coupling(name='TT2_KC+', controlPoint=region1, 
    surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
    localCsys=orientation, u1=ON, u2=ON, u3=OFF, ur1=ON, ur2=ON, ur3=ON)

a = mdb.models['Model-1'].rootAssembly

a.rotate(instanceList=('Torque_Tube_1','Torque_Tube_2', ), axisPoint=(-ORTT,0, 0), 
    axisDirection=(ORTT,0, 0), angle=(90))		
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_1'].datums[DP_TT1].pointOn
a.translate(instanceList=('Torque_Tube_1', ), vector=(DVct[0]-DVct2[0], DVct[1]-DVct2[1], DVct[2]-DVct2[2]))
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP12].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_2'].datums[DP_TT1].pointOn
a.rotate(instanceList=('Torque_Tube_2', ), axisPoint=(0,0,-ORTT), 
axisDirection=(0,0,ORTT), angle=(180))	
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP9].pointOn	

a.translate(instanceList=('Torque_Tube_2', ), vector=(DVct[0]-DVct2[0], DVct[1]-DVct2[1], DVct[2]-DVct2[2]))

a = mdb.models['Model-1'].rootAssembly
region2=a.instances['Right_Plate'].sets['Set-1']
a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[RP17], )
region1=regionToolset.Region(referencePoints=refPoints1)
mdb.models['Model-1'].RigidBody(name='RB_R', refPointRegion=region1, bodyRegion=region2)

a = mdb.models['Model-1'].rootAssembly
region2=a.instances['Left_Plate'].sets['Set-1']
a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[RP18], )
region1=regionToolset.Region(referencePoints=refPoints1)
mdb.models['Model-1'].RigidBody(name='RB_L', refPointRegion=region1, bodyRegion=region2)

a = mdb.models['Model-1'].rootAssembly
region2=a.instances['Middle_Plate'].sets['Set-1']
a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[RP19], )
region1=regionToolset.Region(referencePoints=refPoints1)
mdb.models['Model-1'].RigidBody(name='RB_M', refPointRegion=region1, bodyRegion=region2)

a = mdb.models['Model-1'].rootAssembly
region2=a.instances['Right_Plate_2'].sets['Set-1']
a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[RP24], )
region1=regionToolset.Region(referencePoints=refPoints1)
mdb.models['Model-1'].RigidBody(name='RB_R_2', refPointRegion=region1, bodyRegion=region2)

a = mdb.models['Model-1'].rootAssembly
region2=a.instances['Left_Plate_2'].sets['Set-1']
a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[RP25], )
region1=regionToolset.Region(referencePoints=refPoints1)
mdb.models['Model-1'].RigidBody(name='RB_L_2', refPointRegion=region1, bodyRegion=region2)

for rowname in ['']:
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

# mdb.models['Model-1'].CoupledTempDisplacementStep(name='Temp-Disp', 
    # previous='Initial', timePeriod=20.0, maxNumInc=1000000, initialInc=0.1, 
    # minInc=1e-08, maxInc=0.1, deltmx=100.0)
# mdb.models['Model-1'].steps['Temp-Disp'].setValues(nlgeom=ON)
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'E', 'PE', 'LE', 'U', 'RF', 'CF','NT','COORD','SDV'))
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(
    timeInterval=1.0)
	
#Define BCs				
print 'Defining all BCs'
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP8].pointOn
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Elastomer'].faces 
edges1 = e1.findAt(((DVct[0]-Pad,DVct[1],DVct2[2]+0.001), ))
e2 = a.instances['Elastomer'].faces
edges2 = e1.findAt(((-DVct[0]+Pad,DVct[1],DVct2[2]+0.001), ))
region = a.Set(faces=edges1+edges2, name='Set-43')
mdb.models['Model-1'].EncastreBC(name='Fix_Elast', createStepName='Initial', 
	region=region, localCsys=None)

a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[RP27], r1[RP26], r1[RP28], r1[RP29], )
region = regionToolset.Region(referencePoints=refPoints1)
mdb.models['Model-1'].DisplacementBC(name='Fixed hinge Support', createStepName='Initial', 
	region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
	amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
	localCsys=None)
a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[RP21], r1[RP22], )
region = regionToolset.Region(referencePoints=refPoints1)
mdb.models['Model-1'].DisplacementBC(name='Fixed TT1', createStepName='Initial', 
    region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)	

## BCs##
mdb.models['Model-1'].SmoothStepAmplitude(name='Amp-1', timeSpan=STEP, data=((
    0.0, (Ti)), (2.0, (Ti+0.2*(T1-Ti))), (4.0, (Ti+0.4*(T1-Ti))), (6.0, (Ti+0.6*(T1-Ti))), (8.0, (Ti+0.8*(T1-Ti))), (10.0, (Ti+1.0*(T1-Ti)))))
mdb.models['Model-1'].SmoothStepAmplitude(name='Amp-3', timeSpan=STEP, data=((
    0.0, (Ti)), (2.0, (Ti+0.2*(T2-Ti))), (4.0, (Ti+0.4*(T2-Ti))), (6.0, (Ti+0.6*(T2-Ti))), (8.0, (Ti+0.8*(T2-Ti))), (10.0, (Ti+1.0*(T2-Ti)))))
a = mdb.models['Model-1'].rootAssembly
region = a.instances['Torque_Tube_1'].sets['TT']
mdb.models['Model-1'].Temperature(name='Predefined Field-1', 
    createStepName='Initial', region=region, distributionType=UNIFORM, 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(Ti, ))
mdb.models['Model-1'].predefinedFields['Predefined Field-1'].setValuesInStep(
    stepName='RBM', amplitude='Amp-1', magnitudes=(1, ))
a = mdb.models['Model-1'].rootAssembly
region = a.instances['Torque_Tube_2'].sets['TT']
mdb.models['Model-1'].Temperature(name='Predefined Field-3', 
    createStepName='Initial', region=region, distributionType=UNIFORM, 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(Ti, ))	
mdb.models['Model-1'].predefinedFields['Predefined Field-3'].setValuesInStep(
    stepName='RBM', amplitude='Amp-3', magnitudes=(1, ))	
	
# Pressure Load
a = mdb.models['Model-1'].rootAssembly
p=mdb.models['Model-1'].parts['Elastomer']
s1 = mdb.models['Model-1'].parts['Elastomer'].faces
p.Surface(side1Faces=[s1.findAt(((Presloc[i,0],Presloc[i,1],Presloc[i,2]), ),  ) for i in range(len(Presloc))], name='Elastomer_Surf')
#	side1Faces1 = [s1.findAt(((Presloc[i,0],Presloc[i,1],Presloc[i,2]), ),  ) for i in range(len(Presloc)]
region = a.instances['Elastomer'].surfaces['Elastomer_Surf']
mdb.models['Model-1'].Pressure(name='Aero_Pressure', createStepName='RBM', 
	region=region, distributionType=UNIFORM, field='', magnitude=13750.0, 
	amplitude=UNSET)
	
#Define Sets
print 'Defining Sets'
a = mdb.models['Model-1'].rootAssembly
a.regenerate()
e1 = a.instances['Torque_Tube_1'].elements
elements1 = e1[:]
e3 = a.instances['Torque_Tube_2'].elements
elements3 = e3[:]
e4 = a.instances['Elastomer'].elements
elements4 = e4[:]
a.Set(elements=elements1+elements3+elements4, name='ALL_PART')	
	
a.Set(elements=elements1, name='TT1_ELEM')
a.Set(elements=elements3, name='TT2_ELEM')
a.Set(elements=elements4, name='ELAST_ELEM')

a = mdb.models['Model-1'].rootAssembly
# v5 = a.instances['Left_Plate_2'].vertices
# DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP8].pointOn
# DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP5].pointOn
# verts5 = v5.findAt(((DVct[0], DVct[1],DVct[2]), ), ((DVct2[0], DVct2[1],DVct2[2]), ))
# v1 = a.instances['Left_Plate'].vertices
# DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP8].pointOn
# DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP5].pointOn
# verts1 = v1.findAt(((DVct[0], DVct[1],DVct[2]), ), ((DVct2[0], DVct2[1],DVct2[2]), ))
# v3 = a.instances['Middle_Plate'].vertices
# DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP8].pointOn
# DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP5].pointOn
# verts3 = v3.findAt(((DVct[0], DVct[1],DVct[2]), ), ((DVct2[0], DVct2[1],DVct2[2]), ))
# v4 = a.instances['Right_Plate'].vertices
# DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP8].pointOn
# DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP5].pointOn
# verts4 = v4.findAt(((DVct[0], DVct[1],DVct[2]), ), ((DVct2[0], DVct2[1],DVct2[2]), ))
# a = mdb.models['Model-1'].rootAssembly
# v6 = a.instances['Right_Plate_2'].vertices
# DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP8].pointOn
# DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP5].pointOn
# verts6 = v6.findAt(((DVct[0], DVct[1],DVct[2]), ), ((DVct2[0], DVct2[1],DVct2[2]), ))
# a.Set(vertices=verts5+verts1+verts3+verts4+verts6, name='TIPNODE')

# a = mdb.models['Model-1'].rootAssembly
# e1 = a.instances['Elastomer'].edges
# DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP8].pointOn
# edges1 = e1.findAt(((DVct[0], DVct[1], DVct[2]-elast_thk), ))
# a.Set(edges=edges1, name='TIPNODE')

# x=np.linspace(-0.5,0.5,100)
# DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP8].pointOn
# y=DVct[2]
# p = mdb.models['Model-1'].parts['Elastomer']
# for i in range(100):
	# p.DatumPointByCoordinate(coords=(x[i], y-elast_thk, 0))
	
# elast=getResults2('Model-1')
# np.savetxt("finalz.csv", finalz, delimiter=",")

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
    userSubroutine='C:\\Users\\raoulchandnani\\WorkingDirectory\\UMAT_3D_Coupled_ML_IP_Original.for', 
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
# tipDisp=[0 for x in range(21)]
# # if job.status==COMPLETED:
	# # tipDisp = getResults(ModelName)
# # else:
	# # tipDisp=[9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999]
# tipDisp = getResults(ModelName)
[initial,final]=getResults2(ModelName)
TT1_THE= find_strain(ModelName, 'RBM','THE',SetName = 'TT1_ELEM', object = 'assembly', object_name = None, strain_metric = 'principal')
TT2_THE= find_strain(ModelName, 'RBM','THE',SetName = 'TT2_ELEM', object = 'assembly', object_name = None, strain_metric = 'principal')
TT1_EE= find_strain(ModelName, 'RBM','EE',SetName = 'TT1_ELEM', object = 'assembly', object_name = None, strain_metric = 'Shear')
TT2_EE= find_strain(ModelName, 'RBM','EE',SetName = 'TT2_ELEM', object = 'assembly', object_name = None, strain_metric = 'Shear')
ELAST_EE= find_strain(ModelName, 'RBM','EE',SetName = 'TT2_ELEM', object = 'assembly', object_name = None, strain_metric = 'principal')
TT1_S=find_stress(ModelName,'TT1_ELEM')
TT2_S=find_stress(ModelName,'TT2_ELEM')
ELAST_S=find_stress(ModelName,'ELAST_ELEM')

A=(final[2],final[12],final[3],final[13],final[0],final[10],final[1],final[11],final[8],final[18],final[9],final[19],final[6],final[16],final[7],final[17],final[4],final[14],final[5],final[15])
x=[0 for ind in range(12)]
x[0]=-0.5*Total
x[1:11]=A[0:19:2]
x[11]=0.5*Total
z=[0 for ind in range(12)]
z[0]=initial[12]
z[1:10]=A[1:20:2]
z[11]=initial[15]
fileobject = open('temp.txt','wb')
for xloc in x:
	fileobject.write('%.4f\n' % xloc)
for zloc in z:
	fileobject.write('%.4f\n' % zloc)
fileobject.close()
	
DataFile = open('PostData.txt','a')
DataFile.write('%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f\n' % (IRTT,ORTT,T1,T2,G,HNG,
phiL2,phiR2,tipDisp[2],tipDisp[12],tipDisp[3],tipDisp[13],tipDisp[0],tipDisp[10],tipDisp[1],tipDisp[11],tipDisp[8],tipDisp[18],tipDisp[9],tipDisp[19],tipDisp[6],tipDisp[16],tipDisp[7],tipDisp[17],
tipDisp[4],tipDisp[14],tipDisp[5],tipDisp[15],tipDisp[20],TT1_THE,TT2_THE,TT1_EE,TT2_EE,ELAST_EE,TT1_S,TT2_S,ELAST_S, ))
DataFile.close()
DataFile = open('PostDataInitial.txt','a')
DataFile.write('%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f\n' % (IRTT,ORTT,T1,T2,G,HNG,phiL2,phiR2,initial[2],initial[12],initial[3],initial[13],initial[0],initial[10],initial[1],initial[11],initial[8],initial[18],initial[9],initial[19],initial[6],initial[16],initial[7],initial[17],initial[4],initial[14],initial[5],initial[15],tipDisp[20],TT1_THE,TT2_THE,TT1_EE,TT2_EE,ELAST_EE,TT1_S,TT2_S,ELAST_S, ))
DataFile.close()
DataFile = open('PostDataFinal.txt','a')
DataFile.write('%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,' % (IRTT,ORTT,T1,T2,G,HNG,phiL2,phiR2,final[2],final[12],final[3],final[13],final[0],final[10],final[1],final[11],final[8],final[18],final[9],final[19],final[6],final[16],final[7],final[17],final[4],final[14],final[5],final[15],tipDisp[20],TT1_THE,TT2_THE,TT1_EE,TT2_EE,ELAST_EE,TT1_S,TT2_S,ELAST_S, ))
DataFile.close()
print 'DONE!!'