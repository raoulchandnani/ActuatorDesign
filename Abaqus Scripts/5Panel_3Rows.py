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
V=0.08; # Vertical plate length
G=0.030;   # Gap length
HNG=0.030; # Hinge height
D=0.015; #Hinge distance from edge
phiL=0.2; #initial angle
phiR=0.2;
theta_m=0.0;
theta_l=-0.4;
theta_l2=0.4;
theta_ml=0.0;
theta_ll=-0.3;
theta_l2l=0.3;
theta_mt=0.0;
theta_lt=-0.3;
theta_l2t=0.3;
elast_thk=0.020
Total = 1.000
seedsize=G/4
LD=0.1
theta=asin((sin(phiR)-sin(phiL))/3)
Hnet=(Total)/(3*cos(theta)+cos(phiL)+cos(phiR));
H=Hnet-G;
test=1
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

for name in ['Middle_Plate','Right_Plate','Right_Plate_2','Left_Plate','Left_Plate_2']:
	mdb.models['Model-1'].rootAssembly.features.changeKey(
	fromName='%s-lin-1-2'%name, toName='%sT'%name)
for name in ['Middle_Plate','Right_Plate','Right_Plate_2','Left_Plate','Left_Plate_2']:
	mdb.models['Model-1'].rootAssembly.features.changeKey(
	fromName='%s-lin-1-2-1'%name, toName='%sL'%name)

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
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Middle to Left 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#14,15
RP7=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP8=RP.id
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
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Left to Left  2
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))															#22,23
RP15=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP16=RP.id
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
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP21=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP22=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP23=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP24=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1], DVct[2]))
RP25=RP.id

DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly													# RPs for Middle to Right 1
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))															#8,9
RPT1=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT2=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly													# RPs for Middle to Right 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))													#10,11
RPT3=RP.id															
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT4=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Middle to Left 1
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))															#12,13
RPT5=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT6=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Middle to Left 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))															#14,15
RPT7=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT8=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Right to Right
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))															#16,17
RPT9=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT10=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Right to Right 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))															#18,19
RPT11=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT12=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Left to Left 
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))															#RP13,21
RPT13=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT14=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Left to Left  2
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))															#22,23
RPT15=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT16=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Rigid Body Constraints								#24,25,26
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT17=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT18=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT19=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly																	# RPs for Hinge Support at endswith								#27,28,29,30
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT20=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT21=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT22=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT23=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT24=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]+LD, DVct[2]))
RPT25=RP.id

DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly													# RPs for Middle to Right 1
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))															#8,9
RPL1=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL2=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly													# RPs for Middle to Right 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))													#10,11
RPL3=RP.id															
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL4=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Middle to Left 1
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))															#12,13
RPL5=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL6=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Middle to Left 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))															#14,15
RPL7=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL8=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Right to Right
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))															#16,17
RPL9=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL10=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Right to Right 2
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))															#18,19
RPL11=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL12=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Left to Left 
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))															#RP13,21
RPL13=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL14=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly			# RPs for Left to Left  2
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))															#22,23
RPL15=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL16=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly				# RPs for Rigid Body Constraints								#24,25,26
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL17=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL18=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL19=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly																	# RPs for Hinge Support at endswith								#27,28,29,30
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL20=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL21=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL22=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL23=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL24=RP.id
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP13].pointOn
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(DVct[0],DVct[1]-LD, DVct[2]))
RPL25=RP.id

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
mdb.models['Model-1'].materials['Aluminium'].Elastic(table=((70000.0, 0.3), ))

#Create Elastomer
mdb.models['Model-1'].Material(name='Elastomer')
mdb.models['Model-1'].materials['Elastomer'].Hyperelastic(
    materialType=ISOTROPIC, testData=OFF, type=MOONEY_RIVLIN, 
    volumetricResponse=VOLUMETRIC_DATA, table=((0.3339E06,-0.000337E06, 0.0015828E-06), ))

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
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP22], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
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
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2T'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPT22], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
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
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2L'].datums[DP3].pointOn
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RPL22], v5.findAt(coordinates=(DVct[0],DVct[1], DVct[2]))), ), mergeType=IMPRINT, meshable=OFF)
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

a = mdb.models['Model-1'].rootAssembly
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP10].pointOn
a.DatumCsysByThreePoints(name='CSYS_MRT', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for Middle to Right 1
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_MR2T', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for Middle to R2
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_MLT', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for Middle to Left
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateT'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_ML2T', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for Middle to L2
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateT'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_RST', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for  RHS Support
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateT'].datums[DP9].pointOn	
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_RS2T', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for RHS Support 2
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateT'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_LST', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateT'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_LS2T', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2

a = mdb.models['Model-1'].rootAssembly
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP10].pointOn
a.DatumCsysByThreePoints(name='CSYS_MRL', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for Middle to Right 1
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP9].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_MR2L', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for Middle to R2
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_MLL', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for Middle to Left
DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_PlateL'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_ML2L', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for Middle to L2
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateL'].datums[DP10].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_RSL', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for  RHS Support
DVct=mdb.models['Model-1'].rootAssembly.instances['Right_PlateL'].datums[DP9].pointOn	
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_RS2L', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]-V, DVct[2]), point2=(0.0, 0.0, DVct[2]))			#CSYS for RHS Support 2
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateL'].datums[DP11].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_LSL', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support
DVct=mdb.models['Model-1'].rootAssembly.instances['Left_PlateL'].datums[DP12].pointOn
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_LS2L', coordSysType=CARTESIAN, origin=(DVct[0],DVct[1], DVct[2]), point1=(DVct[0],DVct[1]+V, DVct[2]), point2=(0.0, 0.0, DVct[2]))		#CSYS for LHS Support 2

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
datum1 = a.datums[a.features['CSYS_RS'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-6']
datum1 = a.datums[a.features['CSYS_RS2'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-7']
datum1 = a.datums[a.features['CSYS_LS'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-8']
datum1 = a.datums[a.features['CSYS_LS2'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-Beams']
csa = a.SectionAssignment(sectionName='BEAM', region=region)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-1T']
datum1 = a.datums[a.features['CSYS_MRT'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-2T']
datum1 = a.datums[a.features['CSYS_MR2T'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-3T']
datum1 = a.datums[a.features['CSYS_MLT'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-4T']
datum1 = a.datums[a.features['CSYS_ML2T'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-5T']
datum1 = a.datums[a.features['CSYS_RST'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-6T']
datum1 = a.datums[a.features['CSYS_RS2T'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-7T']
datum1 = a.datums[a.features['CSYS_LST'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-8T']
datum1 = a.datums[a.features['CSYS_LS2T'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-1L']
datum1 = a.datums[a.features['CSYS_MRL'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-2L']
datum1 = a.datums[a.features['CSYS_MR2L'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-3L']
datum1 = a.datums[a.features['CSYS_MLL'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-4L']
datum1 = a.datums[a.features['CSYS_ML2L'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-5L']
datum1 = a.datums[a.features['CSYS_RSL'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-6L']
datum1 = a.datums[a.features['CSYS_RS2L'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-7L']
datum1 = a.datums[a.features['CSYS_LSL'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['Wire-8L']
datum1 = a.datums[a.features['CSYS_LS2L'].id]
csa = a.SectionAssignment(sectionName='HINGE', region=region)
a.ConnectorOrientation(region=csa.getSet(), localCsys1=datum1)

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

#Create Rigid Body Constraints
print 'Defining RBC'
for name in ['','T','L']:
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Right_Plate%s'%name].sets['Set-1']
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refpt=eval('RP%s17'%name)
	refPoints1=(r1[refpt], )
	region1=regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].RigidBody(name='RB_R%s'%name, refPointRegion=region1, bodyRegion=region2)
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Left_Plate%s'%name].sets['Set-1']
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refpt=eval('RP%s18'%name)
	refPoints1=(r1[refpt], )
	region1=regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].RigidBody(name='RB_L%s'%name, refPointRegion=region1, bodyRegion=region2)
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Middle_Plate%s'%name].sets['Set-1']
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refpt=eval('RP%s19'%name)
	refPoints1=(r1[refpt], )
	region1=regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].RigidBody(name='RB_M%s'%name, refPointRegion=region1, bodyRegion=region2)
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Right_Plate_2%s'%name].sets['Set-1']
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refpt=eval('RP%s24'%name)
	refPoints1=(r1[refpt], )
	region1=regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].RigidBody(name='RB_R_2%s'%name, refPointRegion=region1, bodyRegion=region2)
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Left_Plate_2%s'%name].sets['Set-1']
	a = mdb.models['Model-1'].rootAssembly
	r1 = a.referencePoints
	refpt=eval('RP%s25'%name)
	refPoints1=(r1[refpt], )
	region1=regionToolset.Region(referencePoints=refPoints1)
	mdb.models['Model-1'].RigidBody(name='RB_L_2%s'%name, refPointRegion=region1, bodyRegion=region2)

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

for name in ['','T','L']:
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate_2%s'%name].datums[DP13].pointOn
	a = mdb.models['Model-1'].rootAssembly
	s1 = a.instances['Elastomer'].faces
	side2Faces1 = s1.findAt(((DVct[0],DVct[1],DVct[2]), ))
	#region1=a.Surface(side2Faces=side2Faces1, name='s_Surf-1')
	#region1=a.Set(faces=side2Faces1, name='s_Surf-1')
	region1=a.instances['Elastomer'].surfaces['Surf_L2%s'%name]
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Left_Plate_2%s'%name].surfaces['Surf-Bottom']
	mdb.models['Model-1'].Tie(name='Tie_L2%s'%name, master=region2, slave=region1, 
		positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
	DVct=mdb.models['Model-1'].rootAssembly.instances['Left_Plate%s'%name].datums[DP13].pointOn
	a = mdb.models['Model-1'].rootAssembly
	s1 = a.instances['Elastomer'].faces
	side2Faces1 = s1.findAt(((DVct[0],DVct[1],DVct[2]), ))
	#region1=a.Surface(side2Faces=side2Faces1, name='s_Surf-2')
	#region1=a.Set(faces=side2Faces1, name='s_Surf-2')
	region1=a.instances['Elastomer'].surfaces['Surf_L%s'%name]
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Left_Plate%s'%name].surfaces['Surf-Bottom']
	mdb.models['Model-1'].Tie(name='Tie_L%s'%name, master=region2, slave=region1, 
		positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)		
	DVct=mdb.models['Model-1'].rootAssembly.instances['Middle_Plate%s'%name].datums[DP13].pointOn
	a = mdb.models['Model-1'].rootAssembly
	s1 = a.instances['Elastomer'].faces
	side2Faces1 = s1.findAt(((DVct[0],DVct[1],DVct[2]), ))
	#region1=a.Surface(side2Faces=side2Faces1, name='s_Surf-3')
	#region1=a.Set(faces=side2Faces1, name='s_Surf-3')
	region1=a.instances['Elastomer'].surfaces['Surf_M%s'%name]
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Middle_Plate%s'%name].surfaces['Surf-Bottom']
	mdb.models['Model-1'].Tie(name='Tie_M%s'%name, master=region2, slave=region1, 
		positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)		
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate%s'%name].datums[DP13].pointOn
	a = mdb.models['Model-1'].rootAssembly
	s1 = a.instances['Elastomer'].faces
	side2Faces1 = s1.findAt(((DVct[0],DVct[1],DVct[2]), ))
	#region1=a.Surface(side2Faces=side2Faces1, name='s_Surf-4')
	#region1=a.Set(faces=side2Faces1, name='s_Surf-4')
	region1=a.instances['Elastomer'].surfaces['Surf_R%s'%name]
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Right_Plate%s'%name].surfaces['Surf-Bottom']
	mdb.models['Model-1'].Tie(name='Tie_R%s'%name, master=region2, slave=region1, 
		positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)		
	DVct=mdb.models['Model-1'].rootAssembly.instances['Right_Plate_2%s'%name].datums[DP13].pointOn
	a = mdb.models['Model-1'].rootAssembly
	s1 = a.instances['Elastomer'].faces
	side2Faces1 = s1.findAt(((DVct[0],DVct[1],DVct[2]), ))
	#region1=a.Surface(side2Faces=side2Faces1, name='s_Surf-5')
	#region1=a.Set(faces=side2Faces1, name='s_Surf-5')
	region1=a.instances['Elastomer'].surfaces['Surf_R2%s'%name]
	a = mdb.models['Model-1'].rootAssembly
	region2=a.instances['Right_Plate_2%s'%name].surfaces['Surf-Bottom']
	mdb.models['Model-1'].Tie(name='Tie_R2%s'%name, master=region2, slave=region1, 
		positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)


#Define Steps
print 'Defining the Steps'
mdb.models['Model-1'].StaticStep(name='RBM', previous='Initial')
#mdb.models['Model-1'].ImplicitDynamicsStep(name='RBM', previous='Initial')
mdb.models['Model-1'].steps['RBM'].setValues(nlgeom=ON)
mdb.models['Model-1'].steps['RBM'].setValues(solutionTechnique=QUASI_NEWTON)
mdb.models['Model-1'].StaticStep(name='RBM2', previous='RBM')
#mdb.models['Model-1'].ImplicitDynamicsStep(name='RBM', previous='Initial')
mdb.models['Model-1'].steps['RBM2'].setValues(nlgeom=ON)
mdb.models['Model-1'].steps['RBM2'].setValues(solutionTechnique=QUASI_NEWTON)

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
mdb.models['Model-1'].EncastreBC(name='Fix_Elast', createStepName='RBM', 
	region=region, localCsys=None)

a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[RP20], r1[RP21], r1[RP22], r1[RP23],r1[RPL20], r1[RPL21], r1[RPL22], r1[RPL23],r1[RPT20], r1[RPT21], r1[RPT22], r1[RPT23], )
region = regionToolset.Region(referencePoints=refPoints1)
mdb.models['Model-1'].DisplacementBC(name='Fixed hinge Support', createStepName='Initial', 
	region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
	amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
	localCsys=None)

a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[RP22], r1[RP23], )
region = regionToolset.Region(referencePoints=refPoints1)
mdb.models['Model-1'].DisplacementBC(name='Theta_L2', createStepName='RBM', 
	region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=theta_l2l, ur3=UNSET, 
	amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
	localCsys=None)	
a = mdb.models['Model-1'].rootAssembly
region=a.sets['LLH']
mdb.models['Model-1'].ConnDisplacementBC(name='Theta_L', createStepName='RBM', 
	region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=theta_ll, ur2=UNSET, ur3=UNSET, 
	amplitude=UNSET, distributionType=UNIFORM)	

a = mdb.models['Model-1'].rootAssembly
region=a.sets['MLH']
mdb.models['Model-1'].ConnDisplacementBC(name='Theta_M', createStepName='RBM', 
	region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=theta_ml, ur2=UNSET, ur3=UNSET, 
	amplitude=UNSET, distributionType=UNIFORM)

a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[RPL22], r1[RPL23], )
region = regionToolset.Region(referencePoints=refPoints1)
mdb.models['Model-1'].DisplacementBC(name='Theta_L2L', createStepName='RBM', 
	region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=theta_l2l, ur3=UNSET, 
	amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
	localCsys=None)	
a = mdb.models['Model-1'].rootAssembly
region=a.sets['LLHL']
mdb.models['Model-1'].ConnDisplacementBC(name='Theta_LL', createStepName='RBM', 
	region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=theta_ll, ur2=UNSET, ur3=UNSET, 
	amplitude=UNSET, distributionType=UNIFORM)	

a = mdb.models['Model-1'].rootAssembly
region=a.sets['MLHL']
mdb.models['Model-1'].ConnDisplacementBC(name='Theta_ML', createStepName='RBM', 
	region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=theta_ml, ur2=UNSET, ur3=UNSET, 
	amplitude=UNSET, distributionType=UNIFORM)
	
a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[RPT22], r1[RPT23], )
region = regionToolset.Region(referencePoints=refPoints1)
mdb.models['Model-1'].DisplacementBC(name='Theta_L2T', createStepName='RBM', 
	region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=theta_l2t, ur3=UNSET, 
	amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
	localCsys=None)	
a = mdb.models['Model-1'].rootAssembly
region=a.sets['LLHT']
mdb.models['Model-1'].ConnDisplacementBC(name='Theta_LT', createStepName='RBM', 
	region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=theta_lt, ur2=UNSET, ur3=UNSET, 
	amplitude=UNSET, distributionType=UNIFORM)	

a = mdb.models['Model-1'].rootAssembly
region=a.sets['MLHT']
mdb.models['Model-1'].ConnDisplacementBC(name='Theta_MT', createStepName='RBM', 
	region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=theta_mt, ur2=UNSET, ur3=UNSET, 
	amplitude=UNSET, distributionType=UNIFORM)

mdb.models['Model-1'].boundaryConditions['Theta_L2'].setValuesInStep(
    stepName='RBM2', ur2=theta_l2)
mdb.models['Model-1'].boundaryConditions['Theta_L'].setValuesInStep(
    stepName='RBM2', ur1=theta_l)
	
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