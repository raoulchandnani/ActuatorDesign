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
#
# Variables for Plates
H=0.55; # Horizontal plate length 
V=0.125; # Vertical plate length
G=0.040;   # Gap length
HNG=0.030; # Hinge height
D=0.015; #Hinge distance from edge
phiL=0.2; #initial angle
phiR=0.2;
theta_m=0.0;
theta_l=-0.4;
theta_l2=0.4;
elast_thk=0.005
Total = 1.000
seedsize=8
Case=2
# Variables for TT
ORTT=0.02
IRTT=0.01
theta=asin((sin(phiR)-sin(phiL))/3)
Hnet=(Total)/(3*cos(theta)+cos(phiL)+cos(phiR));
H=Hnet-G;
T1=340
T2=340
T3=340
##########################


### Write data file column headings

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
DataFile = open('PostData.txt','w')
DataFile.write('X1		Y1		X2		Y2		X3		Y3		X4		Y4		X5		Y5		X6		Y6		X7		Y7		X8		Y8		X9		Y9		X10		Y10		MaxVM\n')
DataFile.close()

###Sketch Geometry and Create Parts
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
a = mdb.models['Model-1'].rootAssembly				# RPs for Middle to Left 1
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#12,13
RP5=RP.id
if Case==2:
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	RP6=RP.id
else:
	r1 = a.referencePoints
	refPoints1=(r1[RP5], )
	a.Set(referencePoints=refPoints1, name='RP_Z3-')

DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Middle to Left 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#14,15
RP8=RP.id
if Case==2:
	a = mdb.models['Model-1'].rootAssembly
	RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
	RP7=RP.id
else:
	r1 = a.referencePoints
	refPoints1=(r1[RP8], )
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
# a = mdb.models['Model-1'].rootAssembly			# RPs for Left to Left 
# RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#RP13,21
# RP13=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP14=RP.id
r1 = a.referencePoints
refPoints1=(r1[RP14], )
a.Set(referencePoints=refPoints1, name='RP_Z2-')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Left to Left  2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#22,23
RP15=RP.id
# a = mdb.models['Model-1'].rootAssembly
# RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
# RP16=RP.id
r1 = a.referencePoints
refPoints1=(r1[RP15], )
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
if Case==2:
	r1 = a.referencePoints
	refPoints1=(r1[RP20], )
	a.Set(referencePoints=refPoints1, name='RP_Z3+')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP21=RP.id
if Case==2:
	r1 = a.referencePoints
	refPoints1=(r1[RP21], )
	a.Set(referencePoints=refPoints1, name='RP_Z3-')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP22=RP.id
r1 = a.referencePoints
refPoints1=(r1[RP22], )
a.Set(referencePoints=refPoints1, name='RP_Z1-')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP23=RP.id
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
p.BaseSolidExtrude(sketch=s, depth=V)
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
a.translate(instanceList=('Elastomer', ), vector=(0,0.5*V,0))
#Defining the face partitions
print 'Partitioning part'
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
f = p.faces
pickedFaces = f.findAt(((DVct[0],DVct[2], DVct[1]+0.5*V), ))
v, e, d = p.vertices, p.edges, p.datums
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP5].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP6].pointOn
p.PartitionFaceByShortestPath(point1=v.findAt(coordinates=(DVct[0],DVct[2], DVct[1]+0.5*V)), 
	point2=v.findAt(coordinates=(DVct2[0],DVct2[2], DVct2[1]+0.5*V)), faces=pickedFaces)
	
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
f = p.faces
pickedFaces = f.findAt(((DVct[0],DVct[2], DVct[1]+0.5*V), ))
v, e, d = p.vertices, p.edges, p.datums
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP8].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP7].pointOn
p.PartitionFaceByShortestPath(point1=v.findAt(coordinates=(DVct[0],DVct[2], DVct[1]+0.5*V)), 
	point2=v.findAt(coordinates=(DVct2[0],DVct2[2], DVct2[1]+0.5*V)), faces=pickedFaces)
	
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
f = p.faces
pickedFaces = f.findAt(((DVct[0],DVct[2], DVct[1]+0.5*V), ))
v, e, d = p.vertices, p.edges, p.datums
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP5].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP6].pointOn
p.PartitionFaceByShortestPath(point1=v.findAt(coordinates=(DVct[0],DVct[2], DVct[1]+0.5*V)), 
	point2=v.findAt(coordinates=(DVct2[0],DVct2[2], DVct2[1]+0.5*V)), faces=pickedFaces)

DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
f = p.faces
pickedFaces = f.findAt(((DVct[0],DVct[2], DVct[1]+0.5*V), ))
v, e, d = p.vertices, p.edges, p.datums
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP8].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP7].pointOn
p.PartitionFaceByShortestPath(point1=v.findAt(coordinates=(DVct[0],DVct[2], DVct[1]+0.5*V)), 
	point2=v.findAt(coordinates=(DVct2[0],DVct2[2], DVct2[1]+0.5*V)), faces=pickedFaces)


# Create Material
print 'Creating the Materials'
# Create Aluminium
mdb.models['Model-1'].Material(name='Aluminium')
mdb.models['Model-1'].materials['Aluminium'].Elastic(table=((70000.0, 0.3), ))

#Create Elastomer
mdb.models['Model-1'].Material(name='Elastomer')
mdb.models['Model-1'].materials['Elastomer'].Hyperelastic(
    materialType=ISOTROPIC, testData=OFF, type=MOONEY_RIVLIN, 
    volumetricResponse=VOLUMETRIC_DATA, table=((0.69,0.173, 0.0124), ))

# mdb.models['Model-1'].Material(name='Elastomer')
# mdb.models['Model-1'].materials['Elastomer'].Elastic(table=((3.0, 0.4999999), ))

## Create Material##
mdb.models['Model-1'].Material(name='SMA')
mdb.models['Model-1'].materials['SMA'].Density(table=((6500.0, ), ))
mdb.models['Model-1'].materials['SMA'].Depvar(n=31)
mdb.models['Model-1'].materials['SMA'].UserMaterial(mechanicalConstants=(6500.00
,6.62E+10,2.56E+10,0.33,300,0.00000000,0.00000000,1.50E+08,9017722.07,8342551.97
,470.00,470.00,285.000,310.000,315.000,330.000,0.048202422,1.29E-08,5.00E+07
,0.00E+00,0.00E+00,0.00E+00,2.50E+07,5.00E+07,0.9500,0.9500,0.9500,0.9500
,1.00E-06,1.00E-09,0.9999,6,0,0,1,0,0,0))
mdb.models['Model-1'].materials['SMA'].Expansion(table=((0.0, ), ))
mdb.models['Model-1'].materials['SMA'].Conductivity(table=((10.0, ), ))
mdb.models['Model-1'].materials['SMA'].SpecificHeat(table=((0.32, ), ))

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
f = p.faces
faces=f[:]
#faces = f.getByBoundingBox(-3*H,-3V,-G,H,V,G)
region = p.Set(faces=faces, name='Set-1')
p = mdb.models['Model-1'].parts['Elastomer']
p.SectionAssignment(region=region, sectionName='Stretchable', offset=0.0, 
	offsetType=TOP_SURFACE, offsetField='', 
	thicknessAssignment=FROM_SECTION)

## Create Section TT##
mdb.models['Model-1'].HomogeneousSolidSection(name='SMA_Section', material='SMA', 
    thickness=None)
p = mdb.models['Model-1'].parts['Torque_Tube']
v1, e = p.vertices, p.edges
DP=p.DatumCsysByThreePoints(point1=v1.findAt(coordinates=(ORTT, 0.0,-(0.5*V-D))), 
    name='Datum csys-1', coordSysType=CYLINDRICAL, origin=(0,0,-(0.5*V-D)), 
    point2=(ORTT,ORTT,-(0.5*V-D)))
DPCYS=DP.id
p = mdb.models['Model-1'].parts['Torque_Tube']
c = p.cells
cells = c.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), 0), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), 
    -0.5*(IRTT+ORTT)*cos(pi/4), 0), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 0), ), ((
    0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4),0), ))
region = regionToolset.Region(cells=cells)
orientation = mdb.models['Model-1'].parts['Torque_Tube'].datums[DPCYS]
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
p.SectionAssignment(region=region, sectionName='SMA_Section', offset=0.0, 
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
if Case==2:
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP4].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[RP7], v1.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP8], v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
if Case==2:
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP2].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[RP6], v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
# DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP4].pointOn
# a = mdb.models['Model-1'].rootAssembly
# r11 = a.referencePoints
# a.WirePolyLine(points=((v3.findAt(coordinates=(DVct[0],DVct[1], DVct[2])), r1[RP16]), ), mergeType=IMPRINT, meshable=OFF)			
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
# DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP2].pointOn
# a = mdb.models['Model-1'].rootAssembly
# r11 = a.referencePoints
# a.WirePolyLine(points=((r1[RP13], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)		
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP15], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)		
if Case==1:
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP2].pointOn
	a = mdb.models['Model-1'].rootAssembly
	r11 = a.referencePoints
	a.WirePolyLine(points=((r1[RP20], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP1].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP21], v4.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)	
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP22], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
# DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP4].pointOn
# a = mdb.models['Model-1'].rootAssembly
# r11 = a.referencePoints
# a.WirePolyLine(points=((r1[RP23], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)

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
if Case==2:
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
# DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP11].pointOn
# a = mdb.models['Model-1'].rootAssembly
# r11 = a.referencePoints
# a.WirePolyLine(points=((r11[RP13], r11[RP14]), ), mergeType=IMPRINT, meshable=OFF)
# a = mdb.models['Model-1'].rootAssembly
# e1 = a.edges
# edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
# a.Set(edges=edges1, name='Wire-7')
# DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP12].pointOn
# a = mdb.models['Model-1'].rootAssembly
# r11 = a.referencePoints
# a.WirePolyLine(points=((r11[RP15], r11[RP16]), ), mergeType=IMPRINT, meshable=OFF)
# a = mdb.models['Model-1'].rootAssembly
# e1 = a.edges
# edges1 = e1.findAt(((DVct[0],DVct[1], DVct[2]), ))
# a.Set(edges=edges1, name='Wire-8')

a=mdb.models['Model-1'].rootAssembly
a.SetByBoolean(name='MRH', sets=(a.sets['Wire-1'], a.sets['Wire-2'], ))
if Case==2:
	a=mdb.models['Model-1'].rootAssembly
	a.SetByBoolean(name='MLH', sets=(a.sets['Wire-3'], a.sets['Wire-4'], ))

a=mdb.models['Model-1'].rootAssembly
a.SetByBoolean(name='RRH', sets=(a.sets['Wire-5'], a.sets['Wire-6'], ))

# a=mdb.models['Model-1'].rootAssembly
# a.SetByBoolean(name='LLH', sets=(a.sets['Wire-7'], a.sets['Wire-8'], ))

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
a.DatumCsysByThreePoints(name='CSYS_RS', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for  RHS Support
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP9].pointOn	
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_RS2', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for RHS Support 2
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_LS', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP12].pointOn
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
if Case==2:
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
datum1 = a.datums[a.features['CSYS_RS'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-6']
datum1 = a.datums[a.features['CSYS_RS2'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-7']
# datum1 = a.datums[a.features['CSYS_LS'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['Wire-8']
# datum1 = a.datums[a.features['CSYS_LS2'].id]
# csa = a.SectionAssignment(sectionName='HINGE', region=region)
# a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-Beams']
csa = a.SectionAssignment(sectionName='BEAM', region=region)

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
## Mesh##
elemType1 = mesh.ElemType(elemCode=C3D8T, elemLibrary=STANDARD, 
    secondOrderAccuracy=OFF, distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=C3D6T, elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D4T, elemLibrary=STANDARD)
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

a = mdb.models['Model-1'].rootAssembly
a.regenerate()
e6 = a.instances['Elastomer'].elements
elements6 = e6[:]
a.Set(elements=elements6, name='ALL_PART')

#Assemble TTs
p = mdb.models['Model-1'].Part(name='Torque_Tube1', 
    objectToCopy=mdb.models['Model-1'].parts['Torque_Tube'])

p = mdb.models['Model-1'].parts['Torque_Tube1']
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
a.Instance(name='Torque_Tube_2', part=p, dependent=ON)

p = mdb.models['Model-1'].Part(name='Torque_Tube2', 
    objectToCopy=mdb.models['Model-1'].parts['Torque_Tube'])
p = mdb.models['Model-1'].parts['Torque_Tube2']
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
a.Instance(name='Torque_Tube_3', part=p, dependent=ON)


## Constraints##
a = mdb.models['Model-1'].rootAssembly
region1=a.sets['RP_Z1-']
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Torque_Tube_1'].faces
faces1 = f1.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 
    -(0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ))
region2=a.Set(faces=faces1, name='s_Set-1-')
mdb.models['Model-1'].Coupling(name='TT1_KC-', controlPoint=region1, 
    surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
    localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

a = mdb.models['Model-1'].rootAssembly
region1=a.sets['RP_Z1+']
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Torque_Tube_1'].faces
faces1 = f1.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), (0.5*V-D)), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 
    (0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), (0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), (0.5*V-D)), ))
region2=a.Set(faces=faces1, name='s_Set-1+')
mdb.models['Model-1'].Coupling(name='TT1_KC+', controlPoint=region1, 
    surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
    localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

a = mdb.models['Model-1'].rootAssembly
region1=a.sets['RP_Z2-']
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Torque_Tube_2'].faces
faces1 = f1.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 
    -(0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ))
region2=a.Set(faces=faces1, name='s_Set-2-')
mdb.models['Model-1'].Coupling(name='TT2_KC-', controlPoint=region1, 
    surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
    localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

a = mdb.models['Model-1'].rootAssembly
region1=a.sets['RP_Z2+']
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Torque_Tube_2'].faces
faces1 = f1.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), (0.5*V-D)), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 
    (0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), (0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), (0.5*V-D)), ))
region2=a.Set(faces=faces1, name='s_Set-2+')
mdb.models['Model-1'].Coupling(name='TT2_KC+', controlPoint=region1, 
    surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
    localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

a = mdb.models['Model-1'].rootAssembly
region1=a.sets['RP_Z3-']
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Torque_Tube_3'].faces
faces1 = f1.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 
    -(0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ))
region2=a.Set(faces=faces1, name='s_Set-3-')
mdb.models['Model-1'].Coupling(name='TT3_KC-', controlPoint=region1, 
    surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
    localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

a = mdb.models['Model-1'].rootAssembly
region1=a.sets['RP_Z3+']
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Torque_Tube_3'].faces
faces1 = f1.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), (0.5*V-D)), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 
    (0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), (0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), (0.5*V-D)), ))
region2=a.Set(faces=faces1, name='s_Set-3+')
mdb.models['Model-1'].Coupling(name='TT3_KC+', controlPoint=region1, 
    surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
    localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Torque_Tube_1','Torque_Tube_2','Torque_Tube_3', ), axisPoint=(-ORTT,0, 0), 
	axisDirection=(ORTT,0, 0), angle=(90))		
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_1'].datums[DP_TT1].pointOn
a.translate(instanceList=('Torque_Tube_1', ), vector=(DVct[0]-DVct2[0], DVct[1]-DVct2[1], DVct[2]-DVct2[2]))
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP12].pointOn
DVct2=mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_2'].datums[DP_TT1].pointOn
a.translate(instanceList=('Torque_Tube_2', ), vector=(DVct[0]-DVct2[0], DVct[1]-DVct2[1], DVct[2]-DVct2[2]))
DVct2=mdb.models['Model-1'].rootAssembly.instances['Torque_Tube_3'].datums[DP_TT1].pointOn
if Case==1:
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP12].pointOn
else:
	a.rotate(instanceList=('Torque_Tube_3', ), axisPoint=(-ORTT,0, 0), 
	axisDirection=(ORTT,0, 0), angle=(180))	
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP9].pointOn	

a.translate(instanceList=('Torque_Tube_3', ), vector=(DVct[0]-DVct2[0], DVct[1]-DVct2[1], DVct[2]-DVct2[2]))
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

p = mdb.models['Model-1'].parts['Plate']
s = p.elements
side2Elements = s[:]
p.Surface(side2Elements=side2Elements, name='Surf-Bottom')


DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], DVct[1]+0.5*V), ))
p.Surface(side1Faces=side1Faces, name='Surf_L2')
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], DVct[1]+0.5*V), ))
p.Surface(side1Faces=side1Faces, name='Surf_L')
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], DVct[1]+0.5*V), ))
p.Surface(side1Faces=side1Faces, name='Surf_M')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], DVct[1]+0.5*V), ))
p.Surface(side1Faces=side1Faces, name='Surf_R')
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP13].pointOn
p = mdb.models['Model-1'].parts['Elastomer']
s = p.faces
side1Faces = s.findAt(((DVct[0], DVct[2], DVct[1]+0.5*V), ))
p.Surface(side1Faces=side1Faces, name='Surf_R2')


DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['Elastomer'].faces
side2Faces1 = s1.findAt(((DVct[0],DVct[1],DVct[2]), ))
#region1=a.Surface(side2Faces=side2Faces1, name='s_Surf-1')
#region1=a.Set(faces=side2Faces1, name='s_Surf-1')
region1=a.instances['Elastomer'].surfaces['Surf_L2']
a = mdb.models['Model-1'].rootAssembly
region2=a.instances['Left_Plate_2'].surfaces['Surf-Bottom']
mdb.models['Model-1'].Tie(name='Tie_L2', master=region2, slave=region1, 
    positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['Elastomer'].faces
side2Faces1 = s1.findAt(((DVct[0],DVct[1],DVct[2]), ))
#region1=a.Surface(side2Faces=side2Faces1, name='s_Surf-2')
#region1=a.Set(faces=side2Faces1, name='s_Surf-2')
region1=a.instances['Elastomer'].surfaces['Surf_L']
a = mdb.models['Model-1'].rootAssembly
region2=a.instances['Left_Plate'].surfaces['Surf-Bottom']
mdb.models['Model-1'].Tie(name='Tie_L', master=region2, slave=region1, 
    positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
	
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['Elastomer'].faces
side2Faces1 = s1.findAt(((DVct[0],DVct[1],DVct[2]), ))
#region1=a.Surface(side2Faces=side2Faces1, name='s_Surf-3')
#region1=a.Set(faces=side2Faces1, name='s_Surf-3')
region1=a.instances['Elastomer'].surfaces['Surf_M']
a = mdb.models['Model-1'].rootAssembly
region2=a.instances['Middle_Plate'].surfaces['Surf-Bottom']
mdb.models['Model-1'].Tie(name='Tie_M', master=region2, slave=region1, 
    positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
	
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['Elastomer'].faces
side2Faces1 = s1.findAt(((DVct[0],DVct[1],DVct[2]), ))
#region1=a.Surface(side2Faces=side2Faces1, name='s_Surf-4')
#region1=a.Set(faces=side2Faces1, name='s_Surf-4')
region1=a.instances['Elastomer'].surfaces['Surf_R']
a = mdb.models['Model-1'].rootAssembly
region2=a.instances['Right_Plate'].surfaces['Surf-Bottom']
mdb.models['Model-1'].Tie(name='Tie_R', master=region2, slave=region1, 
    positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
	
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['Elastomer'].faces
side2Faces1 = s1.findAt(((DVct[0],DVct[1],DVct[2]), ))
#region1=a.Surface(side2Faces=side2Faces1, name='s_Surf-5')
#region1=a.Set(faces=side2Faces1, name='s_Surf-5')
region1=a.instances['Elastomer'].surfaces['Surf_R2']
a = mdb.models['Model-1'].rootAssembly
region2=a.instances['Right_Plate_2'].surfaces['Surf-Bottom']
mdb.models['Model-1'].Tie(name='Tie_R2', master=region2, slave=region1, 
    positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)


#Define Steps
print 'Defining the Steps'
# mdb.models['Model-1'].StaticStep(name='RBM', previous='Initial')
# mdb.models['Model-1'].steps['RBM'].setValues(nlgeom=ON)
# mdb.models['Model-1'].steps['RBM'].setValues(solutionTechnique=QUASI_NEWTON)

mdb.models['Model-1'].CoupledTempDisplacementStep(name='Temp-Disp', 
    previous='Initial', timePeriod=20.0, maxNumInc=1000000, initialInc=0.1, 
    minInc=1e-08, maxInc=0.1, deltmx=100.0)
mdb.models['Model-1'].steps['Temp-Disp'].setValues(nlgeom=ON)
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'E', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 
    'NT', 'HFL', 'RFL', 'HBF', 'MISES', 'SDV', 'TEMP'))
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
if Case==1:
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refPoints1=(r1[RP20], r1[RP21], r1[RP23], )
	region = regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].DisplacementBC(name='Fixed hinge Support', createStepName='Initial', 
		region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
		amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
		localCsys=None)
else:
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refPoints1=(r1[RP20], r1[RP23], )
	region = regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].DisplacementBC(name='Fixed hinge Support', createStepName='Initial', 
		region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
		amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
		localCsys=None)	

# a = mdb.models['Model-1'].rootAssembly
# r1 = a.referencePoints
# refPoints1=(r1[RP22], r1[RP23], )
# region = regionToolset.Region(referencePoints=refPoints1)
# mdb.models['Model-1'].DisplacementBC(name='Theta_L2', createStepName='RBMRBM', 
	# region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=theta_l2, ur3=UNSET, 
	# amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
	# localCsys=None)	
# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['LLH']
# mdb.models['Model-1'].ConnDisplacementBC(name='Theta_L', createStepName='RBM', 
	# region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=theta_l, ur2=UNSET, ur3=UNSET, 
	# amplitude=UNSET, distributionType=UNIFORM)	

# a = mdb.models['Model-1'].rootAssembly
# region=a.sets['MLH']
# mdb.models['Model-1'].ConnDisplacementBC(name='Theta_M', createStepName='RBM', 
	# region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=theta_m, ur2=UNSET, ur3=UNSET, 
	# amplitude=UNSET, distributionType=UNIFORM)

## BCs##

a = mdb.models['Model-1'].rootAssembly
region = a.instances['Torque_Tube_1'].sets['TT']
mdb.models['Model-1'].Temperature(name='Predefined Field-1', 
    createStepName='Initial', region=region, distributionType=UNIFORM, 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(280.0, ))

a = mdb.models['Model-1'].rootAssembly
region = a.instances['Torque_Tube_2'].sets['TT']
mdb.models['Model-1'].Temperature(name='Predefined Field-2', 
    createStepName='Initial', region=region, distributionType=UNIFORM, 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(280.0, ))

a = mdb.models['Model-1'].rootAssembly
region = a.instances['Torque_Tube_3'].sets['TT']
mdb.models['Model-1'].Temperature(name='Predefined Field-3', 
    createStepName='Initial', region=region, distributionType=UNIFORM, 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(280.0, ))
	
mdb.models['Model-1'].TabularAmplitude(name='Amp-1', timeSpan=STEP, 
    smooth=SOLVER_DEFAULT, data=((0.0, 280.0), (20.0, T1)))
mdb.models['Model-1'].TabularAmplitude(name='Amp-2', timeSpan=STEP, 
    smooth=SOLVER_DEFAULT, data=((0.0, 280.0), (20.0, T2)))
mdb.models['Model-1'].TabularAmplitude(name='Amp-3', timeSpan=STEP, 
    smooth=SOLVER_DEFAULT, data=((0.0, 280.0), (20.0, T3)))
	
a = mdb.models['Model-1'].rootAssembly
region = a.instances['Torque_Tube_1'].sets['TT']
mdb.models['Model-1'].TemperatureBC(name='Temp_change1', createStepName='Temp-Disp', 
    region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    magnitude=1.0, amplitude='Amp-1')
a = mdb.models['Model-1'].rootAssembly
region = a.instances['Torque_Tube_2'].sets['TT']
mdb.models['Model-1'].TemperatureBC(name='Temp_change2', createStepName='Temp-Disp', 
    region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    magnitude=1.0, amplitude='Amp-2')
a = mdb.models['Model-1'].rootAssembly
region = a.instances['Torque_Tube_3'].sets['TT']
mdb.models['Model-1'].TemperatureBC(name='Temp_change3', createStepName='Temp-Disp', 
    region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    magnitude=1.0, amplitude='Amp-3')
	
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
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, 
    userSubroutine='C:\\Temp\\Torque Tube Panel- Raoul\\Working Database\\UMAT_3D_Coupled_ML_IP_Original.for', 
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
    numGPUs=0)
job=mdb.jobs[ModelName]
# delete lock file, which for some reason tends to hang around, if it exists
if os.access('%s.lck'%ModelName,os.F_OK):
	os.remove('%s.lck'%ModelName)	
# Run the job, then process the results.
staho 
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