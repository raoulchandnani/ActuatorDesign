from abaqus import *
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
from Post_P_Script import getResults
##########################
# Variables
H=0.550; # Horizontal plate length 
V=0.1275; # Vertical plate length
G=0.030;   # Gap length
HNG=0.030; # Hinge height
D=0.015; #Hinge distance from edge
phiL=0.2; #initial angle
phiR=0.2;
# theta_m=0.0;
# theta_l=-0.4;
# theta_l2=0.4;
# theta_ml=0.0;
# theta_ll=-0.3;
# theta_l2l=0.3;
# theta_mt=0.0;
# theta_lt=-0.3;
# theta_l2t=0.3;
elast_thk=0.005
Total = 1.000
seedsize=G/4
LD=0.1575
theta=asin((sin(phiR)-sin(phiL))/3)
Hnet=(Total)/(3*cos(theta)+cos(phiL)+cos(phiR));
H=Hnet-G;
test=1
Case=1
# Variables for TT
T1=0.49
T2=0.49
T3=0.22
# T1T=0.35
# T2T=0.39
# T3T=0.26
# T1L=0.35
# T2L=0.39
# T3L=0.26
T1T=T1
T1L=T1
T2T=T2
T2L=T2
T3T=T3
T3L=T3
T12=0.76
T22=0.76
T32=0.28
ORTT=0.02
IRTT=0.01
##########################


### Write data file column headings

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
DataFile = open('PostData.txt','w')
DataFile.write('X1		Y1		X2		Y2		X3		Y3		X4		Y4		X5		Y5		X6		Y6		X7		Y7		X8		Y8		X9		Y9		X10		Y10		MaxVM\n')
DataFile.close()

###Sketch Geometry and Create Parts
print 'Sketching/Creating the part'
Mdb()
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
	sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=((-0.5*H), (-0.5*V)), point2=(0.5*H,0.5*V))
p = mdb.models['Model-1'].Part(name='Plate', dimensionality=THREE_D, 
	type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Plate']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Plate']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']

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
	
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
	sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=((-0.5*H), (-0.5*V)), point2=(0.5*H,0.5*V))
p = mdb.models['Model-1'].Part(name='Plate', dimensionality=THREE_D, 
	type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Plate']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Plate']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']

#Defining the face partitions
print 'Partitioning part'
p = mdb.models['Model-1'].parts['Plate']
f1, e2, d2 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f1.findAt(coordinates=(0, 
	0, 0.0), normal=(0.0, 0.0, 1.0)), sketchUpEdge=e2.findAt(
	coordinates=(0.5*H,0, 0.0)), sketchPlaneSide=SIDE1, origin=(0, 0,0.0))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
	sheetSize=459.88, gridSpacing=11.49, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Plate']
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

#s.Line(point1=(0.0, 0.5*V), point2=(0.0, -0.5*V))
#s.Line(point1=(-0.5*H, 0.0), point2=(0.5*H, 0.0))
s.rectangle(point1=(-(0.5*H-D), -(0.5*V-D)), point2=((0.5*H-D), (0.5*V-D)))
p = mdb.models['Model-1'].parts['Plate']
f = p.faces
pickedFaces = f.findAt(((0,0, 0.0), ))
e, d1 = p.edges, p.datums
p.PartitionFaceBySketch(sketchUpEdge=e.findAt(coordinates=(0.5*H,0, 0.0)), 
	faces=pickedFaces, sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Plate']
s = p.faces
side1Faces = s.findAt(((0, 0, 0.0), ))
# side1Faces = s.findAt(((0, 0, 0.0), ), (((0.5*H-0.5*D), 0.25*V, 0.0), ), (((0.5*H-0.5*D), -0.25*V, 0.0), ), ((
    # -(0.5*H-0.5*D), 0.25*V, 0.0), ), ((-(0.5*H-0.5*D), -0.25*V, 0.0), ))
p.Surface(side2Faces=side1Faces, name='Surf-Bottom')

###Datum Points
p = mdb.models['Model-1'].parts['Plate']
DP=p.DatumPointByCoordinate(coords=((0.5*H-D), 0.5*V-D, 0))
DP2=DP.id
p = mdb.models['Model-1'].parts['Plate']
DP=p.DatumPointByCoordinate(coords=((0.5*H-D), -0.5*V+D, 0))
DP1=DP.id
p = mdb.models['Model-1'].parts['Plate']
DP=p.DatumPointByCoordinate(coords=(-(0.5*H-D), 0.5*V-D, 0))
DP3=DP.id
p = mdb.models['Model-1'].parts['Plate']
DP=p.DatumPointByCoordinate(coords=(-(0.5*H-D), -0.5*V+D, 0))
DP4=DP.id
p = mdb.models['Model-1'].parts['Plate']
DP=p.DatumPointByCoordinate(coords=((0.5*H), 0.5*V, 0))
DP6=DP.id
p = mdb.models['Model-1'].parts['Plate']
DP=p.DatumPointByCoordinate(coords=((0.5*H), -0.5*V, 0))
DP5=DP.id
p = mdb.models['Model-1'].parts['Plate']
DP=p.DatumPointByCoordinate(coords=(-(0.5*H), 0.5*V, 0))
DP7=DP.id
p = mdb.models['Model-1'].parts['Plate']
DP=p.DatumPointByCoordinate(coords=(-(0.5*H), -0.5*V, 0))
DP8=DP.id
p = mdb.models['Model-1'].parts['Plate']
DP=p.DatumPointByCoordinate(coords=((0.5*(H+G)), 0.5*V-D, HNG))
DP10=DP.id
p = mdb.models['Model-1'].parts['Plate']
DP=p.DatumPointByCoordinate(coords=((0.5*(H+G)), -0.5*V+D, HNG))
DP9=DP.id
p = mdb.models['Model-1'].parts['Plate']
DP=p.DatumPointByCoordinate(coords=(-(0.5*(H+G)), 0.5*V-D, HNG))
DP11=DP.id
p = mdb.models['Model-1'].parts['Plate']
DP=p.DatumPointByCoordinate(coords=(-(0.5*(H+G)), -0.5*V+D, HNG))
DP12=DP.id
p = mdb.models['Model-1'].parts['Plate']
DP=p.DatumPointByCoordinate(coords=(0, 0, 0))
DP13=DP.id
#Assemble Parts
print 'Placing Parts in Space'
#Create Instances here

a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Plate']
a.Instance(name='Middle_Plate', part=p, dependent=ON)

a = mdb.models['Model-1'].rootAssembly
a.LinearInstancePattern(instanceList=('Middle_Plate', ), direction1=(1.0, 0.0, 
	0.0), direction2=(0.0, 1.0, 0.0), number1=3, number2=1, spacing1=(H+G), 
	spacing2=(V+G))
a = mdb.models['Model-1'].rootAssembly
a.LinearInstancePattern(instanceList=('Middle_Plate', ), direction1=(-1.0, 0.0, 
	0.0), direction2=(0.0, 1.0, 0.0), number1=3, number2=1, spacing1=(H+G), 
	spacing2=(V+G))

mdb.models['Model-1'].rootAssembly.features.changeKey(
	fromName='Middle_Plate-lin-2-1', toName='Right_Plate')
mdb.models['Model-1'].rootAssembly.features.changeKey(
	fromName='Middle_Plate-lin-2-1-1', toName='Left_Plate')
mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='Middle_Plate-lin-3-1', 
	toName='Right_Plate_2')
mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='Middle_Plate-lin-3-1-1', 
	toName='Left_Plate_2')
Vct=[0,0,0]
DVct=[-0.5*Total,-0.5*V+D,0]
for loc in range(3):
	Vct[loc]= DVct[loc] - mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn[loc]
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('Left_Plate_2', ), vector=(Vct[0],Vct[1], Vct[2]))
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Left_Plate_2', ), axisPoint=(DVct[0],DVct[1], DVct[2]), 
	axisDirection=(0.0, -(2*V), 0.0), angle=(phiL*180/pi))	
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP9].pointOn
for loc in range(3):
	Vct[loc]= DVct[loc] - mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP12].pointOn[loc]
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('Left_Plate', ), vector=(Vct[0],Vct[1], Vct[2]))
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Left_Plate', ), axisPoint=(DVct[0],DVct[1], DVct[2]), 
	axisDirection=(0.0, -(2*V), 0.0), angle=(theta*180/pi))		
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP9].pointOn
for loc in range(3):
	Vct[loc]= DVct[loc] - mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP12].pointOn[loc]
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('Middle_Plate', ), vector=(Vct[0],Vct[1], Vct[2]))
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Middle_Plate', ), axisPoint=(DVct[0],DVct[1], DVct[2]), 
	axisDirection=(0.0, -(2*V), 0.0), angle=(theta*180/pi))			
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP9].pointOn
for loc in range(3):
	Vct[loc]= DVct[loc] - mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP12].pointOn[loc]
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('Right_Plate', ), vector=(Vct[0],Vct[1], Vct[2]))
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Right_Plate', ), axisPoint=(DVct[0],DVct[1], DVct[2]), 
	axisDirection=(0.0, -(2*V), 0.0), angle=(theta*180/pi))	

DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP9].pointOn
for loc in range(3):
	Vct[loc]= DVct[loc] - mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP12].pointOn[loc]
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('Right_Plate_2', ), vector=(Vct[0],Vct[1], Vct[2]))
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Right_Plate_2', ), axisPoint=(DVct[0],DVct[1], DVct[2]), 
	axisDirection=(0.0, (2*V), 0.0), angle=(phiR*180/pi))		

a = mdb.models['Model-1'].rootAssembly
a.LinearInstancePattern(instanceList=('Middle_Plate','Right_Plate','Right_Plate_2','Left_Plate','Left_Plate_2', ), direction1=(1.0, 0.0, 
	0.0), direction2=(0.0, 1.0, 0.0), number1=1, number2=2, spacing1=(LD), 
	spacing2=(LD))
a = mdb.models['Model-1'].rootAssembly
a.LinearInstancePattern(instanceList=('Middle_Plate','Right_Plate','Right_Plate_2','Left_Plate','Left_Plate_2', ), direction1=(1.0, 0.0, 
	0.0), direction2=(0.0, -1.0, 0.0), number1=1, number2=2, spacing1=(LD), 
	spacing2=(LD))

for rowname in ['Middle_Plate','Right_Plate','Right_Plate_2','Left_Plate','Left_Plate_2']:
	mdb.models['Model-1'].rootAssembly.features.changeKey(
	fromName='%s-lin-1-2'%rowname, toName='%sT'%rowname)
for rowname in ['Middle_Plate','Right_Plate','Right_Plate_2','Left_Plate','Left_Plate_2']:
	mdb.models['Model-1'].rootAssembly.features.changeKey(
	fromName='%s-lin-1-2-1'%rowname, toName='%sL'%rowname)

p = mdb.models['Model-1'].parts['Torque_Tube']
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
a.Instance(name='Torque_Tube_1', part=p, dependent=ON)
a.Instance(name='Torque_Tube_1T', part=p, dependent=ON)
a.Instance(name='Torque_Tube_1L', part=p, dependent=ON)
#Create Reference Points
print 'Defining Reference Points'
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
a = mdb.models['Model-1'].rootAssembly				# RPs for Middle to Left 1
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#12,13
RP5=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP6=RP.id
if Case==1:
	r1 = a.referencePoints
	refPoints1=(r1[RP6], )
	a.Set(referencePoints=refPoints1, name='RP_Z3-')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Middle to Left 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#14,15
RP7=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP8=RP.id
if Case==1:
	r1 = a.referencePoints
	refPoints1=(r1[RP7], )
	a.Set(referencePoints=refPoints1, name='RP_Z3+')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Right to Right
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#16,17
RP9=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP10=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Right to Right 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#18,19
RP11=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP12=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Left to Left 
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#RP13,21
RP13=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP14=RP.id
r1 = a.referencePoints
refPoints1=(r1[RP13], )
a.Set(referencePoints=refPoints1, name='RP_Z2-')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Left to Left  2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#22,23
RP15=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP16=RP.id
r1 = a.referencePoints
refPoints1=(r1[RP16], )
a.Set(referencePoints=refPoints1, name='RP_Z2+')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Rigid Body Constraints								#24,25,26
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
if Case==2:
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
if Case==2:
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

DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly													# RPs for Middle to Right 1
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#8,9
RPT1=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT2=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly													# RPs for Middle to Right 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))													#10,11
RPT3=RP.id															
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT4=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Middle to Left 1
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#12,13
RPT5=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT6=RP.id
if Case==1:
	r1 = a.referencePoints
	refPoints1=(r1[RPT6], )
	a.Set(referencePoints=refPoints1, name='RPT_Z3-')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Middle to Left 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#14,15
RPT7=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT8=RP.id
if Case==1:
	r1 = a.referencePoints
	refPoints1=(r1[RPT7], )
	a.Set(referencePoints=refPoints1, name='RPT_Z3+')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateT'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Right to Right
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#16,17
RPT9=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT10=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateT'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Right to Right 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#18,19
RPT11=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT12=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateT'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Left to Left 
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#RP13,21
RPT13=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT14=RP.id
r1 = a.referencePoints
refPoints1=(r1[RPT13], )
a.Set(referencePoints=refPoints1, name='RPT_Z2-')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateT'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Left to Left  2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#22,23
RPT15=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT16=RP.id
r1 = a.referencePoints
refPoints1=(r1[RPT16], )
a.Set(referencePoints=refPoints1, name='RPT_Z2+')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateT'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Rigid Body Constraints								#24,25,26
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT17=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateT'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT18=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT19=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2T'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly																	# RPs for Hinge Support at endswith								#27,28,29,30
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT20=RP.id
a = mdb.models['Model-1'].rootAssembly																	# RPs for Hinge Support at endswith								#27,28,29,30
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT26=RP.id
if Case==2:
	r1 = a.referencePoints
	refPoints1=(r1[RPT20], )
	a.Set(referencePoints=refPoints1, name='RPT_Z3+')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2T'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT21=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT27=RP.id
if Case==2:
	r1 = a.referencePoints
	refPoints1=(r1[RPT21], )
	a.Set(referencePoints=refPoints1, name='RPT_Z3-')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2T'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT22=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT28=RP.id
r1 = a.referencePoints
refPoints1=(r1[RPT22], )
a.Set(referencePoints=refPoints1, name='RPT_Z1-')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2T'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT23=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT29=RP.id
r1 = a.referencePoints
refPoints1=(r1[RPT23], )
a.Set(referencePoints=refPoints1, name='RPT_Z1+')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2T'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT24=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2T'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPT25=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP10].pointOn

DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly													# RPs for Middle to Right 1
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#8,9
RPL1=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL2=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly													# RPs for Middle to Right 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))													#10,11
RPL3=RP.id															
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL4=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Middle to Left 1
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#12,13
RPL5=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL6=RP.id
if Case==1:
	r1 = a.referencePoints
	refPoints1=(r1[RPL6], )
	a.Set(referencePoints=refPoints1, name='RPL_Z3-')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Middle to Left 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#14,15
RPL7=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL8=RP.id
if Case==1:
	r1 = a.referencePoints
	refPoints1=(r1[RPL7], )
	a.Set(referencePoints=refPoints1, name='RPL_Z3+')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateL'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Right to Right
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#16,17
RPL9=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL10=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateL'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Right to Right 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#18,19
RPL11=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL12=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateL'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Left to Left 
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#RP13,21
RPL13=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL14=RP.id
r1 = a.referencePoints
refPoints1=(r1[RPL13], )
a.Set(referencePoints=refPoints1, name='RPL_Z2-')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateL'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Left to Left  2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#22,23
RPL15=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL16=RP.id
r1 = a.referencePoints
refPoints1=(r1[RPL16], )
a.Set(referencePoints=refPoints1, name='RPL_Z2+')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateL'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Rigid Body Constraints								#24,25,26
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL17=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateL'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL18=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL19=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2L'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly																	# RPs for Hinge Support at endswith								#27,28,29,30
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL20=RP.id
a = mdb.models['Model-1'].rootAssembly																	# RPs for Hinge Support at endswith								#27,28,29,30
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL26=RP.id
if Case==2:
	r1 = a.referencePoints
	refPoints1=(r1[RPL20], )
	a.Set(referencePoints=refPoints1, name='RPL_Z3+')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2L'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL21=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL27=RP.id
if Case==2:
	r1 = a.referencePoints
	refPoints1=(r1[RPL21], )
	a.Set(referencePoints=refPoints1, name='RPL_Z3-')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2L'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL22=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL28=RP.id
r1 = a.referencePoints
refPoints1=(r1[RPL22], )
a.Set(referencePoints=refPoints1, name='RPL_Z1-')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2L'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL23=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL29=RP.id
r1 = a.referencePoints
refPoints1=(r1[RPL23], )
a.Set(referencePoints=refPoints1, name='RPL_Z1+')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2L'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL24=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2L'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RPL25=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP10].pointOn

## Create Elastomer
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP8].pointOn
s.Line(point1=(DVct[0], DVct2[2]), point2=(DVct[0], DVct2[2]-elast_thk))
s.Line(point1=(DVct[0], DVct2[2]-elast_thk), point2=(DVct2[0], DVct2[2]-elast_thk))
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP5].pointOn
s.Line(point1=(DVct2[0], DVct2[2]-elast_thk), point2=(DVct[0], DVct[2]-elast_thk))
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP8].pointOn
s.Line(point1=(DVct[0], DVct[2]-elast_thk), point2=(DVct2[0], DVct2[2]-elast_thk))
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP5].pointOn
s.Line(point1=(DVct2[0], DVct2[2]-elast_thk), point2=(DVct[0], DVct[2]-elast_thk))
DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP8].pointOn
s.Line(point1=(DVct[0], DVct[2]-elast_thk), point2=(DVct2[0], DVct2[2]-elast_thk))
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP5].pointOn
s.Line(point1=(DVct2[0], DVct2[2]-elast_thk), point2=(DVct[0], DVct[2]-elast_thk))
DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP8].pointOn
s.Line(point1=(DVct[0], DVct[2]-elast_thk), point2=(DVct2[0], DVct2[2]-elast_thk))
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP5].pointOn
s.Line(point1=(DVct2[0], DVct2[2]-elast_thk), point2=(DVct[0], DVct[2]-elast_thk))
DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP8].pointOn
s.Line(point1=(DVct[0], DVct[2]-elast_thk), point2=(DVct2[0], DVct2[2]-elast_thk))
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP5].pointOn
s.Line(point1=(DVct2[0], DVct2[2]-elast_thk), point2=(DVct[0], DVct[2]-elast_thk))
DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP9].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP8].pointOn
DVct3=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
s.Line(point1=(DVct[0], DVct[2]-elast_thk), point2=(-DVct3[0], DVct2[2]-elast_thk))

DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP8].pointOn
s.Line(point1=(DVct[0], DVct2[2]), point2=(DVct2[0], DVct2[2]))
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP5].pointOn
s.Line(point1=(DVct2[0], DVct2[2]), point2=(DVct[0], DVct[2]))
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP8].pointOn
s.Line(point1=(DVct[0], DVct[2]), point2=(DVct2[0], DVct2[2]))
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP5].pointOn
s.Line(point1=(DVct2[0], DVct2[2]), point2=(DVct[0], DVct[2]))
DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP8].pointOn
s.Line(point1=(DVct[0], DVct[2]), point2=(DVct2[0], DVct2[2]))
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP5].pointOn
s.Line(point1=(DVct2[0], DVct2[2]), point2=(DVct[0], DVct[2]))
DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP8].pointOn
s.Line(point1=(DVct[0], DVct[2]), point2=(DVct2[0], DVct2[2]))
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP5].pointOn
s.Line(point1=(DVct2[0], DVct2[2]), point2=(DVct[0], DVct[2]))
DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP8].pointOn
s.Line(point1=(DVct[0], DVct[2]), point2=(DVct2[0], DVct2[2]))
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP5].pointOn
s.Line(point1=(DVct2[0], DVct2[2]), point2=(DVct[0], DVct[2]))
DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP9].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP8].pointOn
DVct3=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
s.Line(point1=(DVct[0], DVct[2]), point2=(-DVct3[0], DVct2[2]))
s.Line(point1=(-DVct3[0], DVct2[2]), point2=(-DVct3[0], DVct2[2]-elast_thk))

p = mdb.models['Model-1'].Part(name='Elastomer', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Elastomer']
p.BaseSolidExtrude(sketch=s, depth=0.5*(V+2*LD))
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP8].pointOn
f1=p.faces
p.Mirror(mirrorPlane=f1.findAt(coordinates=(0,DVct[2]-0.5*elast_thk, 0)), 
    keepOriginal=ON)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']

a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Elastomer']
a.Instance(name='Elastomer', part=p, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Elastomer', ), axisPoint=(0,0, 0), 
	axisDirection=(1, 0, 0.0), angle=(90))	
a = mdb.models['Model-1'].rootAssembly
#a.translate(instanceList=('Elastomer', ), vector=(0,0.5*V,0))

#Defining the face partitions
print 'Partitioning part'
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
f = p.faces
pickedFaces = f.findAt(((DVct[0],DVct[2], DVct[1]), ))
v, e, d = p.vertices, p.edges, p.datums
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateL'].datums[DP5].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_PlateT'].datums[DP6].pointOn
# p.PartitionFaceByShortestPath(point1=v.findAt(coordinates=(DVct[0],DVct[2], DVct[1])), 
	# point2=v.findAt(coordinates=(DVct2[0],DVct2[2], DVct2[1])), faces=pickedFaces)
p.PartitionFaceByShortestPath(point1=(DVct[0],DVct[2], DVct[1]), 
	point2=(DVct2[0],DVct2[2], DVct2[1]), faces=pickedFaces)
	
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
f = p.faces
pickedFaces = f.findAt(((DVct[0],DVct[2], DVct[1]), ))
v, e, d = p.vertices, p.edges, p.datums
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP8].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP7].pointOn
# p.PartitionFaceByShortestPath(point1=v.findAt(coordinates=(DVct[0],DVct[2], DVct[1])), 
	# point2=v.findAt(coordinates=(DVct2[0],DVct2[2], DVct2[1])), faces=pickedFaces)
p.PartitionFaceByShortestPath(point1=(DVct[0],DVct[2], DVct[1]), 
	point2=(DVct2[0],DVct2[2], DVct2[1]), faces=pickedFaces)
	
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
f = p.faces
pickedFaces = f.findAt(((DVct[0],DVct[2], DVct[1]), ))
v, e, d = p.vertices, p.edges, p.datums
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP5].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP6].pointOn
# p.PartitionFaceByShortestPath(point1=v.findAt(coordinates=(DVct[0],DVct[2], DVct[1])), 
	# point2=v.findAt(coordinates=(DVct2[0],DVct2[2], DVct2[1])), faces=pickedFaces)
p.PartitionFaceByShortestPath(point1=(DVct[0],DVct[2], DVct[1]), 
	point2=(DVct2[0],DVct2[2], DVct2[1]), faces=pickedFaces)

DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
f = p.faces
pickedFaces = f.findAt(((DVct[0],DVct[2], DVct[1]), ))
v, e, d = p.vertices, p.edges, p.datums
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateL'].datums[DP8].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_PlateT'].datums[DP7].pointOn
# p.PartitionFaceByShortestPath(point1=v.findAt(coordinates=(DVct[0],DVct[2], DVct[1])), 
	# point2=v.findAt(coordinates=(DVct2[0],DVct2[2], DVct2[1])), faces=pickedFaces)
p.PartitionFaceByShortestPath(point1=(DVct[0],DVct[2], DVct[1]), 
	point2=(DVct2[0],DVct2[2], DVct2[1]), faces=pickedFaces)

DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP13].pointOn
DVct1=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP13].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
DVct3=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP13].pointOn
DVct4=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
f = p.faces
pickedFaces = f.findAt(((DVct[0],DVct[2], DVct[1]), ))
pickedFaces1 = f.findAt(((DVct1[0],DVct1[2], DVct1[1]), ))
pickedFaces2 = f.findAt(((DVct2[0],DVct2[2], DVct2[1]), ))
pickedFaces3 = f.findAt(((DVct3[0],DVct3[2], DVct3[1]), ))
pickedFaces4 = f.findAt(((DVct4[0],DVct4[2], DVct4[1]), ))
v, e, d = p.vertices, p.edges, p.datums
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2L'].datums[DP7].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2L'].datums[DP7].pointOn
# p.PartitionFaceByShortestPath(point1=v.findAt(coordinates=(DVct[0],DVct[2], DVct[1])), 
	# point2=v.findAt(coordinates=(DVct2[0],DVct2[2], DVct2[1])), faces=pickedFaces)
p.Set(faces=pickedFaces+pickedFaces1+pickedFaces2+pickedFaces3+pickedFaces4, name='Upsurf')
# p.PartitionFaceByShortestPath(point1=(DVct[0],DVct[2], DVct[1]), 
	# point2=(DVct2[0],DVct2[2], DVct2[1]), faces=pickedFaces+pickedFaces1+pickedFaces2+pickedFaces3+pickedFaces4)

DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP6].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
f, e1, d2 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f.findAt(coordinates=(DVct2[0],DVct2[2], 
    DVct2[1])), sketchUpEdge=e1.findAt(coordinates=(DVct[0], DVct[2], 
   0)), sketchPlaneSide=SIDE1, origin=(0.0, DVct2[2], 0.0))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=0.91, 
    gridSpacing=0.02, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Elastomer']
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2L'].datums[DP7].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2L'].datums[DP7].pointOn
s.Line(point1=(-0.5*Total, DVct[1]), point2=(0.5*Total, DVct2[1]))
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP8].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP8].pointOn
s.Line(point1=(-0.5*Total, DVct[1]), point2=(0.5*Total, DVct2[1]))
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP7].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP7].pointOn
s.Line(point1=(-0.5*Total, DVct[1]), point2=(0.5*Total, DVct2[1]))
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2T'].datums[DP8].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2T'].datums[DP8].pointOn
s.Line(point1=(-0.5*Total, DVct[1]), point2=(0.5*Total, DVct2[1]))
p = mdb.models['Model-1'].parts['Elastomer']
f = p.faces
f1, e, d1 = p.faces, p.edges, p.datums
DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP6].pointOn
p.PartitionFaceBySketchThruAll(sketchPlane=f1.findAt(coordinates=(DVct2[0],DVct2[2], 
    DVct2[1])), sketchUpEdge=e1.findAt(coordinates=(DVct[0], DVct[2], 0)),
	faces=pickedFaces+pickedFaces1+pickedFaces2+pickedFaces3+pickedFaces4, sketchPlaneSide=SIDE1, sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']	

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

# mdb.models['Model-1'].Material(name='Elastomer')
# mdb.models['Model-1'].materials['Elastomer'].Elastic(table=((3.0, 0.4999999), ))

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

# mdb.models['Model-1'].Material(name='Elastomer')
# mdb.models['Model-1'].materials['Elastomer'].Elastic(table=((3.0, 0.4999999), ))
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

p = mdb.models['Model-1'].parts['Plate']
f = p.faces
faces = f.getByBoundingBox(-H,-V,-G,H,V,G)
region = p.Set(faces=faces, name='Set-1')
p = mdb.models['Model-1'].parts['Plate']
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
DP=p.DatumCsysByThreePoints(point1=v1.findAt(coordinates=(ORTT, 0.0,-(0.5*V-D))), 
    name='Datum csys-1', coordSysType=CYLINDRICAL, origin=(0,0,-(0.5*V-D)), 
    point2=(ORTT,ORTT,-(0.5*V-D)))
DPCYS=DP.id
DP=p.DatumCsysByThreePoints(point1=v1.findAt(coordinates=(ORTT, 0.0,(0.5*V-D))), 
    name='Datum csys-2', coordSysType=CYLINDRICAL, origin=(0,0,(0.5*V-D)), 
    point2=(ORTT,ORTT,(0.5*V-D)))
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
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP5], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP7], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP8], v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP6], v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
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
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP20], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP21], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)	
# DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP3].pointOn
# a = mdb.models['Model-1'].rootAssembly
# r11 = a.referencePoints
# a.WirePolyLine(points=((r1[RP22], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP23], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)

a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
v1 = a.instances['Middle_PlateT'].vertices
v2 = a.instances['Right_PlateT'].vertices
v3 = a.instances['Left_PlateT'].vertices
v4 = a.instances['Right_Plate_2T'].vertices
v5 = a.instances['Left_Plate_2T'].vertices

DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT1], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateT'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT2], v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT3], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateT'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT4], v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT5], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT7], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateT'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT8], v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateT'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT6], v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateT'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[RPT16]), ), mergeType=IMPRINT, meshable=OFF)			
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateT'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[RPT14]), ), mergeType=IMPRINT, meshable=OFF)				
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateT'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[RPT10]), ), mergeType=IMPRINT, meshable=OFF)				
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateT'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[RPT12]), ), mergeType=IMPRINT, meshable=OFF)				
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2T'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT9], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)	
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2T'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT11], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)	
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2T'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT13], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)		
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2T'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT15], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)		
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2T'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT20], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2T'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT21], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)	
# DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2T'].datums[DP3].pointOn
# a = mdb.models['Model-1'].rootAssembly
# r11 = a.referencePoints
# a.WirePolyLine(points=((r1[RPT22], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2T'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT23], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)

a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
v1 = a.instances['Middle_PlateL'].vertices
v2 = a.instances['Right_PlateL'].vertices
v3 = a.instances['Left_PlateL'].vertices
v4 = a.instances['Right_Plate_2L'].vertices
v5 = a.instances['Left_Plate_2L'].vertices

DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL1], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateL'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL2], v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL3], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateL'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL4], v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL5], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL7], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateL'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL8], v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateL'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL6], v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateL'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[RPL16]), ), mergeType=IMPRINT, meshable=OFF)			
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateL'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[RPL14]), ), mergeType=IMPRINT, meshable=OFF)				
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateL'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[RPL10]), ), mergeType=IMPRINT, meshable=OFF)				
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateL'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v2.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[RPL12]), ), mergeType=IMPRINT, meshable=OFF)				
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2L'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL9], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)	
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2L'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL11], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)	
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2L'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL13], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)		
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2L'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL15], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)		
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2L'].datums[DP2].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL20], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2L'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL21], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)	
# DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2L'].datums[DP3].pointOn
# a = mdb.models['Model-1'].rootAssembly
# r11 = a.referencePoints
# a.WirePolyLine(points=((r1[RPL22], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2L'].datums[DP4].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL23], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)

a = mdb.models['Model-1'].rootAssembly
e1 = a.edges 
edges1 = e1.getByBoundingBox(-4*H,-4*V,-2*H,4*H,4*V,2*H)
a.Set(edges=edges1, name='Wire-Beams')

a = mdb.models['Model-1'].rootAssembly
for wires in range(68):
	a.suppressFeatures(('Wire-%i'% (wires+1), ))

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
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPT1], r11[RPT2]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP10].pointOn
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-1T')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPT3], r11[RPT4]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-2T')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPT5], r11[RPT6]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-3T')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPT7], r11[RPT8]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-4T')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateT'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPT9], r11[RPT10]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-5T')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateT'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPT11], r11[RPT12]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-6T')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateT'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPT13], r11[RPT14]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-7T')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateT'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPT15], r11[RPT16]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-8T')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2T'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPT21], r11[RPT27]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-9T')

DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2T'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPT20], r11[RPT26]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-10T')

DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2T'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPT22], r11[RPT28]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-11T')

DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2T'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPT23], r11[RPT29]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-12T')

a=mdb.models['Model-1'].rootAssembly
a.SetByBoolean(name='MRHT', sets=(a.sets['Wire-1T'], a.sets['Wire-2T'], ))

a=mdb.models['Model-1'].rootAssembly
a.SetByBoolean(name='MLHT', sets=(a.sets['Wire-3T'], a.sets['Wire-4T'], ))

a=mdb.models['Model-1'].rootAssembly
a.SetByBoolean(name='RRHT', sets=(a.sets['Wire-5T'], a.sets['Wire-6T'], ))

a=mdb.models['Model-1'].rootAssembly
a.SetByBoolean(name='LLHT', sets=(a.sets['Wire-7T'], a.sets['Wire-8T'], ))

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPL1], r11[RPL2]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP10].pointOn
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-1L')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPL3], r11[RPL4]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-2L')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPL5], r11[RPL6]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-3L')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPL7], r11[RPL8]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-4L')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateL'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPL9], r11[RPL10]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-5L')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateL'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPL11], r11[RPL12]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-6L')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateL'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPL13], r11[RPL14]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-7L')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateL'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPL15], r11[RPL16]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-8L')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2L'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPL21], r11[RPL27]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-9L')

DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2L'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPL20], r11[RPL26]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-10L')

DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2L'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPL22], r11[RPL28]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-11L')

DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2L'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RPL23], r11[RPL29]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
a.Set(edges=edges1, name='Wire-12L')

a=mdb.models['Model-1'].rootAssembly
a.SetByBoolean(name='MRHL', sets=(a.sets['Wire-1L'], a.sets['Wire-2L'], ))

a=mdb.models['Model-1'].rootAssembly
a.SetByBoolean(name='MLHL', sets=(a.sets['Wire-3L'], a.sets['Wire-4L'], ))

a=mdb.models['Model-1'].rootAssembly
a.SetByBoolean(name='RRHL', sets=(a.sets['Wire-5L'], a.sets['Wire-6L'], ))

a=mdb.models['Model-1'].rootAssembly
a.SetByBoolean(name='LLHL', sets=(a.sets['Wire-7L'], a.sets['Wire-8L'], ))

a = mdb.models['Model-1'].rootAssembly
for wires in range(68):
	a.resumeFeatures(('Wire-%i'% (wires+1), ))

#Create and Assign Connector Sections
print 'Create Connector Section'	
mdb.models['Model-1'].ConnectorSection(name='HINGE', assembledType=HINGE)

print 'Create Connector Section'	
mdb.models['Model-1'].ConnectorSection(name='BEAM', assembledType=BEAM)

print 'Define Datum Coordinate Systems'
for rowname in ['','T','L']:
	a = mdb.models['Model-1'].rootAssembly
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP10].pointOn
	a.DatumCsysByThreePoints(name='CSYS_MR%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for Middle to Right 1
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP9].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_MR2%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for Middle to R2
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP11].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_ML%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for Middle to Left
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP12].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_ML2%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for Middle to L2
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP10].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_RR%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for  RHS Support
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP9].pointOn	
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_RR2%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for RHS Support 2
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP11].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_LL%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP12].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_LL2%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2

	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP9].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_RS%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP10].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_RS2%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP11].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_LS%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP12].pointOn
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByThreePoints(name='CSYS_LS2%s'%rowname, coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2

# a = mdb.models['Model-1'].rootAssembly
# DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP10].pointOn
# a.DatumCsysByThreePoints(name='CSYS_MRT', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for Middle to Right 1
# DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP9].pointOn
# a = mdb.models['Model-1'].rootAssembly
# a.DatumCsysByThreePoints(name='CSYS_MR2T', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for Middle to R2
# DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP11].pointOn
# a = mdb.models['Model-1'].rootAssembly
# a.DatumCsysByThreePoints(name='CSYS_MLT', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for Middle to Left
# DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP12].pointOn
# a = mdb.models['Model-1'].rootAssembly
# a.DatumCsysByThreePoints(name='CSYS_ML2T', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for Middle to L2
# DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateT'].datums[DP10].pointOn
# a = mdb.models['Model-1'].rootAssembly
# a.DatumCsysByThreePoints(name='CSYS_RST', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for  RHS Support
# DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateT'].datums[DP9].pointOn	
# a = mdb.models['Model-1'].rootAssembly
# a.DatumCsysByThreePoints(name='CSYS_RS2T', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for RHS Support 2
# DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateT'].datums[DP11].pointOn
# a = mdb.models['Model-1'].rootAssembly
# a.DatumCsysByThreePoints(name='CSYS_LST', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support
# DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateT'].datums[DP12].pointOn
# a = mdb.models['Model-1'].rootAssembly
# a.DatumCsysByThreePoints(name='CSYS_LS2T', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2

# a = mdb.models['Model-1'].rootAssembly
# DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP10].pointOn
# a.DatumCsysByThreePoints(name='CSYS_MRL', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for Middle to Right 1
# DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP9].pointOn
# a = mdb.models['Model-1'].rootAssembly
# a.DatumCsysByThreePoints(name='CSYS_MR2L', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for Middle to R2
# DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP11].pointOn
# a = mdb.models['Model-1'].rootAssembly
# a.DatumCsysByThreePoints(name='CSYS_MLL', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for Middle to Left
# DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP12].pointOn
# a = mdb.models['Model-1'].rootAssembly
# a.DatumCsysByThreePoints(name='CSYS_ML2L', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for Middle to L2
# DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateL'].datums[DP10].pointOn
# a = mdb.models['Model-1'].rootAssembly
# a.DatumCsysByThreePoints(name='CSYS_RSL', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for  RHS Support
# DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateL'].datums[DP9].pointOn	
# a = mdb.models['Model-1'].rootAssembly
# a.DatumCsysByThreePoints(name='CSYS_RS2L', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for RHS Support 2
# DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateL'].datums[DP11].pointOn
# a = mdb.models['Model-1'].rootAssembly
# a.DatumCsysByThreePoints(name='CSYS_LSL', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support
# DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateL'].datums[DP12].pointOn
# a = mdb.models['Model-1'].rootAssembly
# a.DatumCsysByThreePoints(name='CSYS_LS2L', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2

print 'Assign Connector Sections'
for rowname in ['','T','L']:
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

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-1T']
# datum1 = a.datums[a.features['CSYS_MRT'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-2T']
# datum1 = a.datums[a.features['CSYS_MR2T'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-3T']
# datum1 = a.datums[a.features['CSYS_MLT'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-4T']
# datum1 = a.datums[a.features['CSYS_ML2T'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-5T']
# datum1 = a.datums[a.features['CSYS_RST'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-6T']
# datum1 = a.datums[a.features['CSYS_RS2T'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-7T']
# datum1 = a.datums[a.features['CSYS_LST'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-8T']
# datum1 = a.datums[a.features['CSYS_LS2T'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-1L']
# datum1 = a.datums[a.features['CSYS_MRL'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-2L']
# datum1 = a.datums[a.features['CSYS_MR2L'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-3L']
# datum1 = a.datums[a.features['CSYS_MLL'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-4L']
# datum1 = a.datums[a.features['CSYS_ML2L'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-5L']
# datum1 = a.datums[a.features['CSYS_RSL'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-6L']
# datum1 = a.datums[a.features['CSYS_RS2L'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-7L']
# datum1 = a.datums[a.features['CSYS_LSL'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-8L']
# datum1 = a.datums[a.features['CSYS_LS2L'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

#Mesh Parts
print 'Meshing the Part'
p = mdb.models['Model-1'].parts['Plate']
p.seedPart(size=seedsize, deviationFactor=0.1, minSizeFactor=0.1)		#G/6
p = mdb.models['Model-1'].parts['Plate']
p.generateMesh()

elemType1 = mesh.ElemType(elemCode=C3D8H, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['Elastomer']
p.seedPart(size=seedsize, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['Elastomer']
c = p.cells
cells = c[:]
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
    elemType3))
p = mdb.models['Model-1'].parts['Elastomer']
p.generateMesh()

a = mdb.models['Model-1'].rootAssembly
a.regenerate()
e6 = a.instances['Elastomer'].elements
elements6 = e6[:]
a.Set(elements=elements6, name='ALL_PART')

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
p = mdb.models['Model-1'].Part(name='Torque_Tube1', 
    objectToCopy=mdb.models['Model-1'].parts['Torque_Tube'])
p = mdb.models['Model-1'].parts['Torque_Tube1']
region = p.sets['TT']
del mdb.models['Model-1'].parts['Torque_Tube1'].sectionAssignments[0]
p.SectionAssignment(region=region, sectionName='TT_sect1', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)
	
p = mdb.models['Model-1'].parts['Torque_Tube1']
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
a.Instance(name='Torque_Tube_2', part=p, dependent=ON)
a.Instance(name='Torque_Tube_2T', part=p, dependent=ON)
a.Instance(name='Torque_Tube_2L', part=p, dependent=ON)
if Case==1:
	p = mdb.models['Model-1'].parts['Torque_Tube1']
else:
	p = mdb.models['Model-1'].parts['Torque_Tube']	
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
a.Instance(name='Torque_Tube_3', part=p, dependent=ON)
a.Instance(name='Torque_Tube_3T', part=p, dependent=ON)
a.Instance(name='Torque_Tube_3L', part=p, dependent=ON)

## Constraints##
for rowname in ['','T','L']:
	a = mdb.models['Model-1'].rootAssembly
	region1=a.sets['RP%s_Z1-'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	f1 = a.instances['Torque_Tube_1%s'%rowname].faces
	faces1 = f1.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 
		-(0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ))
	region2=a.Set(faces=faces1, name='s_Set-1-%s'%rowname)
	mdb.models['Model-1'].Coupling(name='TT1%s_KC-'%rowname, controlPoint=region1, 
		surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
		localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

	a = mdb.models['Model-1'].rootAssembly
	region1=a.sets['RP%s_Z1+'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	e1 = a.instances['Torque_Tube_1%s'%rowname].edges
	orientation = mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_1%s'%rowname].datums[DPCYS2]
	edges1 = e1.findAt(((-ORTT*cos(pi/4), ORTT*cos(pi/4), (0.5*V-D)), ), ((-ORTT*cos(pi/4), -ORTT*cos(pi/4), 
		(0.5*V-D)), ), ((ORTT*cos(pi/4), ORTT*cos(pi/4), (0.5*V-D)), ), ((ORTT*cos(pi/4), -ORTT*cos(pi/4), (0.5*V-D)), ))
	region2=a.Set(edges=edges1, name='s_Set-1+%s'%rowname)
	mdb.models['Model-1'].Coupling(name='TT1%s_KC+'%rowname, controlPoint=region1, 
		surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
		localCsys=orientation, u1=ON, u2=ON, u3=OFF, ur1=ON, ur2=ON, ur3=ON)

	a = mdb.models['Model-1'].rootAssembly
	region1=a.sets['RP%s_Z2-'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	e1 = a.instances['Torque_Tube_2%s'%rowname].edges
	orientation = mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_2%s'%rowname].datums[DPCYS]
	edges1 = e1.findAt(((-ORTT*cos(pi/4), ORTT*cos(pi/4), -(0.5*V-D)), ), ((-ORTT*cos(pi/4), -ORTT*cos(pi/4), 
		-(0.5*V-D)), ), ((ORTT*cos(pi/4), ORTT*cos(pi/4), -(0.5*V-D)), ), ((ORTT*cos(pi/4), -ORTT*cos(pi/4), -(0.5*V-D)), ))
	region2=a.Set(edges=edges1, name='s_Set-2-%s'%rowname)
	mdb.models['Model-1'].Coupling(name='TT2%s_KC-'%rowname, controlPoint=region1, 
		surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
		localCsys=orientation, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

	a = mdb.models['Model-1'].rootAssembly
	region1=a.sets['RP%s_Z2+'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	orientation = mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_2%s'%rowname].datums[DPCYS2]
	e1 = a.instances['Torque_Tube_2%s'%rowname].edges
	edges1 = e1.findAt(((-ORTT*cos(pi/4), ORTT*cos(pi/4), (0.5*V-D)), ), ((-ORTT*cos(pi/4), -ORTT*cos(pi/4), 
		(0.5*V-D)), ), ((ORTT*cos(pi/4), ORTT*cos(pi/4), (0.5*V-D)), ), ((ORTT*cos(pi/4), -ORTT*cos(pi/4), (0.5*V-D)), ))
	region2=a.Set(edges=edges1, name='s_Set-2+%s'%rowname)
	mdb.models['Model-1'].Coupling(name='TT2%s_KC+'%rowname, controlPoint=region1, 
		surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
		localCsys=orientation, u1=ON, u2=ON, u3=OFF, ur1=ON, ur2=ON, ur3=ON)

	a = mdb.models['Model-1'].rootAssembly
	region1=a.sets['RP%s_Z3-'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	orientation = mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_3%s'%rowname].datums[DPCYS]
	e1 = a.instances['Torque_Tube_3%s'%rowname].edges
	edges1 = e1.findAt(((-ORTT*cos(pi/4), ORTT*cos(pi/4), -(0.5*V-D)), ), ((-ORTT*cos(pi/4), -ORTT*cos(pi/4), 
		-(0.5*V-D)), ), ((ORTT*cos(pi/4), ORTT*cos(pi/4), -(0.5*V-D)), ), ((ORTT*cos(pi/4), -ORTT*cos(pi/4), -(0.5*V-D)), ))
	region2=a.Set(edges=edges1, name='s_Set-3-%s'%rowname)
	mdb.models['Model-1'].Coupling(name='TT3%s_KC-'%rowname, controlPoint=region1, 
		surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
		localCsys=orientation, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

	a = mdb.models['Model-1'].rootAssembly
	region1=a.sets['RP%s_Z3+'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	orientation = mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_3%s'%rowname].datums[DPCYS2]
	e1 = a.instances['Torque_Tube_3%s'%rowname].edges
	edges1 = e1.findAt(((-ORTT*cos(pi/4), ORTT*cos(pi/4), (0.5*V-D)), ), ((-ORTT*cos(pi/4), -ORTT*cos(pi/4), 
		(0.5*V-D)), ), ((ORTT*cos(pi/4), ORTT*cos(pi/4), (0.5*V-D)), ), ((ORTT*cos(pi/4), -ORTT*cos(pi/4), (0.5*V-D)), ))
	region2=a.Set(edges=edges1, name='s_Set-3+%s'%rowname)
	mdb.models['Model-1'].Coupling(name='TT3%s_KC+'%rowname, controlPoint=region1, 
		surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
		localCsys=orientation, u1=ON, u2=ON, u3=OFF, ur1=ON, ur2=ON, ur3=ON)

	a = mdb.models['Model-1'].rootAssembly
	a.rotate(instanceList=('Torque_Tube_1%s'%rowname,'Torque_Tube_2%s'%rowname,'Torque_Tube_3%s'%rowname, ), axisPoint=(-ORTT,0, 0), 
		axisDirection=(ORTT,0, 0), angle=(90))		
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP12].pointOn
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_1%s'%rowname].datums[DP_TT1].pointOn
	a.translate(instanceList=('Torque_Tube_1%s'%rowname, ), vector=(DVct[0]-DVct2[0], DVct[1]-DVct2[1], DVct[2]-DVct2[2]))
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP12].pointOn
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_2%s'%rowname].datums[DP_TT1].pointOn
	a.translate(instanceList=('Torque_Tube_2%s'%rowname, ), vector=(DVct[0]-DVct2[0], DVct[1]-DVct2[1], DVct[2]-DVct2[2]))
	DVct2=mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_3%s'%rowname].datums[DP_TT1].pointOn
	if Case==1:
		DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP12].pointOn
	else:
		a.rotate(instanceList=('Torque_Tube_3%s'%rowname, ), axisPoint=(0,0,-ORTT), 
		axisDirection=(0,0,ORTT), angle=(180))	
		DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP9].pointOn	
	a.translate(instanceList=('Torque_Tube_3%s'%rowname, ), vector=(DVct[0]-DVct2[0], DVct[1]-DVct2[1], DVct[2]-DVct2[2]))

#Create Rigid Body Constraints
print 'Defining RBC'
for rowname in ['','T','L']:
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Right_Plate%s'%rowname].sets['Set-1']
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refpt=eval('RP%s17'%rowname)
	refPoints1=(r1[refpt], )
	region1=regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].RigidBody(name='RB_R%s'%rowname, refPointRegion=region1, bodyRegion=region2)
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Left_Plate%s'%rowname].sets['Set-1']
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refpt=eval('RP%s18'%rowname)
	refPoints1=(r1[refpt], )
	region1=regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].RigidBody(name='RB_L%s'%rowname, refPointRegion=region1, bodyRegion=region2)
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Middle_Plate%s'%rowname].sets['Set-1']
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refpt=eval('RP%s19'%rowname)
	refPoints1=(r1[refpt], )
	region1=regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].RigidBody(name='RB_M%s'%rowname, refPointRegion=region1, bodyRegion=region2)
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Right_Plate_2%s'%rowname].sets['Set-1']
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refpt=eval('RP%s24'%rowname)
	refPoints1=(r1[refpt], )
	region1=regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].RigidBody(name='RB_R_2%s'%rowname, refPointRegion=region1, bodyRegion=region2)
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Left_Plate_2%s'%rowname].sets['Set-1']
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refpt=eval('RP%s25'%rowname)
	refPoints1=(r1[refpt], )
	region1=regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].RigidBody(name='RB_L_2%s'%rowname, refPointRegion=region1, bodyRegion=region2)

p = mdb.models['Model-1'].parts['Plate']
s = p.elements
side2Elements = s[:]
p.Surface(side2Elements=side2Elements, name='Surf-Bottom')

DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], DVct[1]), ))
p.Surface(side1Faces=side1Faces, name='Surf_L2')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], DVct[1]), ))
p.Surface(side1Faces=side1Faces, name='Surf_L')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], DVct[1]), ))
p.Surface(side1Faces=side1Faces, name='Surf_M')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], DVct[1]), ))
p.Surface(side1Faces=side1Faces, name='Surf_R')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], DVct[1]), ))
p.Surface(side1Faces=side1Faces, name='Surf_R2')

DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2T'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], -DVct[1]), ))
p.Surface(side1Faces=side1Faces, name='Surf_L2T')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateT'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], -DVct[1]), ))
p.Surface(side1Faces=side1Faces, name='Surf_LT')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], -DVct[1]), ))
p.Surface(side1Faces=side1Faces, name='Surf_MT')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateT'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], -DVct[1]), ))
p.Surface(side1Faces=side1Faces, name='Surf_RT')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2T'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], -DVct[1]), ))
p.Surface(side1Faces=side1Faces, name='Surf_R2T')

DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2L'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], -DVct[1]), ))
p.Surface(side1Faces=side1Faces, name='Surf_L2L')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateL'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], -DVct[1]), ))
p.Surface(side1Faces=side1Faces, name='Surf_LL')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], -DVct[1]), ))
p.Surface(side1Faces=side1Faces, name='Surf_ML')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateL'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], -DVct[1]), ))
p.Surface(side1Faces=side1Faces, name='Surf_RL')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2L'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], -DVct[1]), ))
p.Surface(side1Faces=side1Faces, name='Surf_R2L')

for rowname in ['','T','L']:
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%rowname].datums[DP13].pointOn
	a = mdb.models['Model-1'].rootAssembly
	s1 = a.instances['Elastomer'].faces
	side2Faces1 = s1.findAt(((DVct[0],DVct[1],DVct[2]), ))
	#region1=a.Surface(side2Faces=side2Faces1, name='s_Surf-1')
	#region1=a.Set(faces=side2Faces1, name='s_Surf-1')
	region1=a.instances['Elastomer'].surfaces['Surf_L2%s'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Left_Plate_2%s'%rowname].surfaces['Surf-Bottom']
	mdb.models['Model-1'].Tie(name='Tie_L2%s'%rowname, master=region2, slave=region1, 
		positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%rowname].datums[DP13].pointOn
	a = mdb.models['Model-1'].rootAssembly
	s1 = a.instances['Elastomer'].faces
	side2Faces1 = s1.findAt(((DVct[0],DVct[1],DVct[2]), ))
	#region1=a.Surface(side2Faces=side2Faces1, name='s_Surf-2')
	#region1=a.Set(faces=side2Faces1, name='s_Surf-2')
	region1=a.instances['Elastomer'].surfaces['Surf_L%s'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Left_Plate%s'%rowname].surfaces['Surf-Bottom']
	mdb.models['Model-1'].Tie(name='Tie_L%s'%rowname, master=region2, slave=region1, 
		positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)		
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%rowname].datums[DP13].pointOn
	a = mdb.models['Model-1'].rootAssembly
	s1 = a.instances['Elastomer'].faces
	side2Faces1 = s1.findAt(((DVct[0],DVct[1],DVct[2]), ))
	#region1=a.Surface(side2Faces=side2Faces1, name='s_Surf-3')
	#region1=a.Set(faces=side2Faces1, name='s_Surf-3')
	region1=a.instances['Elastomer'].surfaces['Surf_M%s'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Middle_Plate%s'%rowname].surfaces['Surf-Bottom']
	mdb.models['Model-1'].Tie(name='Tie_M%s'%rowname, master=region2, slave=region1, 
		positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)		
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%rowname].datums[DP13].pointOn
	a = mdb.models['Model-1'].rootAssembly
	s1 = a.instances['Elastomer'].faces
	side2Faces1 = s1.findAt(((DVct[0],DVct[1],DVct[2]), ))
	#region1=a.Surface(side2Faces=side2Faces1, name='s_Surf-4')
	#region1=a.Set(faces=side2Faces1, name='s_Surf-4')
	region1=a.instances['Elastomer'].surfaces['Surf_R%s'%rowname]
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Right_Plate%s'%rowname].surfaces['Surf-Bottom']
	mdb.models['Model-1'].Tie(name='Tie_R%s'%rowname, master=region2, slave=region1, 
		positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)		
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%rowname].datums[DP13].pointOn
	a = mdb.models['Model-1'].rootAssembly
	s1 = a.instances['Elastomer'].faces
	side2Faces1 = s1.findAt(((DVct[0],DVct[1],DVct[2]), ))
	#region1=a.Surface(side2Faces=side2Faces1, name='s_Surf-5')
	#region1=a.Set(faces=side2Faces1, name='s_Surf-5')
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
mdb.models['Model-1'].ImplicitDynamicsStep(name='RBM2', previous='RBM', 
    application=QUASI_STATIC, nohaf=OFF, amplitude=RAMP, alpha=DEFAULT, 
    initialConditions=OFF,timePeriod=10.0)
mdb.models['Model-1'].steps['RBM2'].setValues(nlgeom=ON)	
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'E', 'PE', 'LE', 'U', 'RF', 'CF','NT','COORD','THE','EE'))

#Define BCs				
print 'Defining all BCs'
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP8].pointOn
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Elastomer'].faces 
edges1 = e1.findAt(((DVct[0],DVct[1],DVct2[2]-0.5*elast_thk), ))
e2 = a.instances['Elastomer'].faces
edges2 = e1.findAt(((-DVct[0],DVct[1],DVct2[2]-0.5*elast_thk), ))
region = a.Set(faces=edges1+edges2, name='Set-43')
mdb.models['Model-1'].EncastreBC(name='Fix_Elast', createStepName='Initial', 
	region=region, localCsys=None)
	
for rowname in ['','T','L']:
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refPoints1=(r1[eval('RP%s27'%rowname)], r1[eval('RP%s26'%rowname)], r1[eval('RP%s28'%rowname)], r1[eval('RP%s29'%rowname)], )
	region = regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].DisplacementBC(name='Fixed hinge Support%s'%rowname, createStepName='Initial', 
		region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
		amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
		localCsys=None)
	if Case==1:
		a = mdb.models['Model-1'].rootAssembly
		r1 = a.referencePoints
		refPoints1=(r1[eval('RP%s22'%rowname)], )
		region = regionToolset.Region(referencePoints=refPoints1)
		mdb.models['Model-1'].DisplacementBC(name='Fixed TT1%s'%rowname, createStepName='Initial', 
			region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
			amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
			localCsys=None)
	else:
		a = mdb.models['Model-1'].rootAssembly
		r1 = a.referencePoints
		refPoints1=(r1[eval('RP%s21'%rowname)], r1[eval('RP%s22'%rowname)], )
		region = regionToolset.Region(referencePoints=refPoints1)
		mdb.models['Model-1'].DisplacementBC(name='Fixed hinge Support%s'%rowname, createStepName='Initial', 
			region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
			amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
			localCsys=None)	
	# a = mdb.models['Model-1'].rootAssembly
	# r1 = a.referencePoints
	# if Case==1:
		# refPoints1=(r1[eval('RP%s1'%rowname)],r1[eval('RP%s2'%rowname)],r1[eval('RP%s3'%rowname)],r1[eval('RP%s4'%rowname)],r1[eval('RP%s6'%rowname)],r1[eval('RP%s7'%rowname)],r1[eval('RP%s9'%rowname)],r1[eval('RP%s10'%rowname)],r1[eval('RP%s11'%rowname)],r1[eval('RP%s12'%rowname)], r1[eval('RP%s13'%rowname)],r1[eval('RP%s16'%rowname)],r1[eval('RP%s17'%rowname)],r1[eval('RP%s18'%rowname)],r1[eval('RP%s19'%rowname)],r1[eval('RP%s23'%rowname)],r1[eval('RP%s24'%rowname)],r1[eval('RP%s25'%rowname)], )
	# else:
		# refPoints1=(r1[eval('RP%s1'%rowname)],r1[eval('RP%s2'%rowname)],r1[eval('RP%s3'%rowname)],r1[eval('RP%s4'%rowname)],r1[eval('RP%s5'%rowname)],r1[eval('RP%s6'%rowname)],r1[eval('RP%s7'%rowname)],r1[eval('RP%s8'%rowname)],r1[eval('RP%s9'%rowname)],r1[eval('RP%s10'%rowname)],r1[eval('RP%s11'%rowname)],r1[eval('RP%s12'%rowname)], r1[eval('RP%s13'%rowname)],r1[eval('RP%s16'%rowname)],r1[eval('RP%s17'%rowname)],r1[eval('RP%s18'%rowname)],r1[eval('RP%s19'%rowname)],r1[eval('RP%s20'%rowname)],r1[eval('RP%s23'%rowname)],r1[eval('RP%s24'%rowname)],r1[eval('RP%s25'%rowname)], )
	# region = regionToolset.Region(referencePoints=refPoints1)
	# mdb.models['Model-1'].DisplacementBC(name='Fix Lateral%s'%rowname, createStepName='RBM', 
		# region=region, u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
		# amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
		# localCsys=None)

## BCs##
mdb.models['Model-1'].SmoothStepAmplitude(name='Amp-1', timeSpan=STEP, data=((
    0.0, 0.0), (2.0, 0.2), (4.0, 0.4), (6.0, 0.6), (8.0, 0.8), (10.0, 1.0)))
mdb.models['Model-1'].SmoothStepAmplitude(name='Amp-2', timeSpan=STEP, data=((
    0.0,T1), (2.0,0.2*(T12-T1)+T1), (4.0,0.4*(T12-T1)+T1), (6.0,0.6*(T12-T1)+T1), (8.0,0.8*(T12-T1)+T1), (10.0,T12)))	
mdb.models['Model-1'].SmoothStepAmplitude(name='Amp-3', timeSpan=STEP, data=((
    0.0,T2), (2.0,0.2*(T22-T2)+T2), (4.0,0.4*(T22-T2)+T2), (6.0,0.6*(T22-T2)+T2), (8.0,0.8*(T22-T2)+T2), (10.0,T22)))	
mdb.models['Model-1'].SmoothStepAmplitude(name='Amp-4', timeSpan=STEP, data=((
    0.0,T3), (2.0,0.2*(T32-T3)+T3), (4.0,0.4*(T32-T3)+T3), (6.0,0.6*(T32-T3)+T3), (8.0,0.8*(T32-T3)+T3), (10.0,T32)))		
# mdb.models['Model-1'].SmoothStepAmplitude(name='Amp-1', timeSpan=STEP, data=((
    # 0.0, 0.0), (1.0, 0.1), (2.0, 0.2), (3.0, 0.3), (4.0, 0.4), (5.0, 0.5), (6.0, 0.6), (7.0, 0.7), (8.0, 0.8), (9.0, 0.9), (10.0, 1.0)))
for rowname in ['','T','L']:	
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
	# mdb.models['Model-1'].predefinedFields['Predefined Field-1'].setValuesInStep(
		# stepName='RBM', magnitudes=(T1, ))
		
	datum = mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_2%s'%rowname].datums[DPCYS]
	mdb.models['Model-1'].ExpressionField(name='AnalyticalField-2%s'%rowname, 
		localCsys=datum, description='', expression='R/%f' % ORTT)
	a = mdb.models['Model-1'].rootAssembly
	region = a.instances['Torque_Tube_2%s'%rowname].sets['TT']
	mdb.models['Model-1'].Temperature(name='Predefined Field-2%s'%rowname, 
		createStepName='Initial', region=region, distributionType=UNIFORM, 
		crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(0.0, ))
	mdb.models['Model-1'].predefinedFields['Predefined Field-2%s'%rowname].setValues(
		distributionType=FIELD, field='AnalyticalField-2%s'%rowname)
	mdb.models['Model-1'].predefinedFields['Predefined Field-2%s'%rowname].setValuesInStep(
		stepName='RBM', amplitude='Amp-1',magnitudes=(eval('T2%s'%rowname), ))
	# mdb.models['Model-1'].predefinedFields['Predefined Field-2'].setValuesInStep(
		# stepName='RBM', magnitudes=(T2, ))
		
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
		
mdb.models['Model-1'].predefinedFields['Predefined Field-1'].setValuesInStep(
	stepName='RBM2', amplitude='Amp-2',magnitudes=(1))
mdb.models['Model-1'].predefinedFields['Predefined Field-2'].setValuesInStep(
	stepName='RBM2', amplitude='Amp-3', magnitudes=(1))		
mdb.models['Model-1'].predefinedFields['Predefined Field-3'].setValuesInStep(
	stepName='RBM2',  amplitude='Amp-4',magnitudes=(1))

## Pressure Load
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['Elastomer'].faces
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP8].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP13].pointOn
DVct3=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP5].pointOn
DVct4=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP8].pointOn
DVct5=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
DVct6=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP5].pointOn
DVct7=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP13].pointOn
DVct8=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP8].pointOn
DVct9=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP5].pointOn
side1Faces1 = s1.findAt(((0.5*(DVct[0]-0.5*Total),0,DVct[2]-elast_thk), ), ((DVct2[0],0,DVct2[2]-elast_thk), ), ((0.5*(DVct3[0]+DVct4[0]),0,0.5*(DVct3[2]+DVct4[2])-elast_thk), ), 
((DVct5[0],0,DVct5[2]-elast_thk), ), ((0.5*(DVct8[0]+DVct9[0]),0,0.5*(DVct8[2]+DVct9[2])-elast_thk), ), ((DVct7[0],0,DVct7[2]-elast_thk), ), ((0.5*(DVct6[0]+0.5*Total),0,DVct6[2]-elast_thk), ))
region = a.Surface(side1Faces=side1Faces1, name='Elastomer_Surf')
mdb.models['Model-1'].Pressure(name='Aero_Pressure', createStepName='RBM', 
    region=region, distributionType=UNIFORM, field='', magnitude=13750.0, 
    amplitude=UNSET)
	
#Define Sets
print 'Defining Sets'

a = mdb.models['Model-1'].rootAssembly
v5 = a.instances['Left_Plate_2'].vertices
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP8].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP5].pointOn
verts5 = v5.findAt(((DVct[0], DVct[1],DVct[2]), ), ((DVct2[0], DVct2[1],DVct2[2]), ))
v1 = a.instances['Left_Plate'].vertices
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP8].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP5].pointOn
verts1 = v1.findAt(((DVct[0], DVct[1],DVct[2]), ), ((DVct2[0], DVct2[1],DVct2[2]), ))
v3 = a.instances['Middle_Plate'].vertices
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP8].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP5].pointOn
verts3 = v3.findAt(((DVct[0], DVct[1],DVct[2]), ), ((DVct2[0], DVct2[1],DVct2[2]), ))
v4 = a.instances['Right_Plate'].vertices
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP8].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP5].pointOn
verts4 = v4.findAt(((DVct[0], DVct[1],DVct[2]), ), ((DVct2[0], DVct2[1],DVct2[2]), ))
a = mdb.models['Model-1'].rootAssembly
v6 = a.instances['Right_Plate_2'].vertices
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP8].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP5].pointOn
verts6 = v6.findAt(((DVct[0], DVct[1],DVct[2]), ), ((DVct2[0], DVct2[1],DVct2[2]), ))
a.Set(vertices=verts5+verts1+verts3+verts4+verts6, name='TIPNODE')


#####################################
### Creation/Execution of the Job ###
#####################################
print 'Creating/Running Job'

ModelName='Model-1'

mdb.Job(name=ModelName, model=ModelName, description='', type=ANALYSIS, 
	atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
	memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
	explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
	modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
	scratch='', multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)

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
tipDisp=[0 for x in range(21)]
if job.status==COMPLETED:
	tipDisp = getResults(ModelName)
else:
	tipDisp=[9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999]
DataFile = open('PostData.txt','a')
DataFile.write('%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f\n' % (tipDisp[2],tipDisp[12],tipDisp[3],tipDisp[13],tipDisp[0],tipDisp[10],tipDisp[1],tipDisp[11],tipDisp[8],tipDisp[18],tipDisp[9],tipDisp[19],tipDisp[6],tipDisp[16],tipDisp[7],tipDisp[17],tipDisp[4],tipDisp[14],tipDisp[5],tipDisp[15],tipDisp[20], ))
DataFile.close()
print 'DONE!!'