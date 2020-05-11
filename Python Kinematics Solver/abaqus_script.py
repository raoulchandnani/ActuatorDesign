# read the files in the input file
fileobject = open('input.txt','rb')
DVs = []
for line in fileobject:
	DVs.append(float(line))
fileobject.close()
# Script an amazing Abaqus script
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
H=193; # Horizontal plate length 
V=125; # Vertical plate length
G=35;   # Gap length
HNG=35; # Hinge height
D=15; #Hinge distance from edge
phiL=0.35; #initial angle
phiR=0.35;
theta_r=0.10;
theta_l=0.10;
elast_thk=5
Total = 3000;
Mdb()
theta_l=DVs[0]
theta_r=DVs[1]
phiL=DVs[2]
phiR=DVs[3]
G=DVs[4]
H=(Total-G*(3+cos(phiL)+cos(phiR)))/(3+cos(phiL)+cos(phiR));
HNG=DVs[5]
D = 5

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

s.Line(point1=(0.0, 0.5*V), point2=(0.0, -0.5*V))
s.Line(point1=(-0.5*H, 0.0), point2=(0.5*H, 0.0))
s.rectangle(point1=(-(0.5*H-D), -(0.5*V-D)), point2=((0.5*H-D), (0.5*V-D)))
p = mdb.models['Model-1'].parts['Plate']
f = p.faces
pickedFaces = f.findAt(((0,0, 0.0), ))
e, d1 = p.edges, p.datums
p.PartitionFaceBySketch(sketchUpEdge=e.findAt(coordinates=(0.5*H,0, 0.0)), 
	faces=pickedFaces, sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']


###Sketch Geometry and Create Elastomer
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
	sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=((-0.5*G), (-0.5*V)), point2=(0.5*G,0.5*V))
p = mdb.models['Model-1'].Part(name='Elastomer1', dimensionality=THREE_D, 
	type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Elastomer1']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Elastomer1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']

###Sketch Geometry and Create Elastomer
LH=sqrt((0.25*G*G)+(HNG*HNG))
phi02=atan2((0.5*G),HNG)		
phi1=2*phi02-phiL
phi2=2*phi02-phiR
LD1=sqrt(2*LH*LH*(1-cos(phi1)))
LD2=sqrt(2*LH*LH*(1-cos(phi2)))
LHNG2=sqrt((0.5*G)*(0.5*G)+(HNG*HNG))
LD3=(G)*cos(phiL)-LHNG2*sin(phi02-phiL)
LD4=(G)*cos(phiR)-LHNG2*sin(phi02-phiR)

s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
	sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=((-0.5*LD1), (-0.5*V)), point2=(0.5*LD1,0.5*V))
p = mdb.models['Model-1'].Part(name='Elastomer2L', dimensionality=THREE_D, 
	type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Elastomer2L']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Elastomer2L']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']

s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
	sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=((-0.5*LD2), (-0.5*V)), point2=(0.5*LD2,0.5*V))
p = mdb.models['Model-1'].Part(name='Elastomer2R', dimensionality=THREE_D, 
	type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Elastomer2R']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Elastomer2R']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']


s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
	sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=((-0.5*LD3), (-0.5*V)), point2=(0.5*LD3,0.5*V))
p = mdb.models['Model-1'].Part(name='Elastomer3L', dimensionality=THREE_D, 
	type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Elastomer3L']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Elastomer3L']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']

s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
	sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=((-0.5*LD4), (-0.5*V)), point2=(0.5*LD4,0.5*V))
p = mdb.models['Model-1'].Part(name='Elastomer3R', dimensionality=THREE_D, 
	type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Elastomer3R']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Elastomer3R']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']

# Create Material
print 'Creating the Materials'
# Create Aluminium
mdb.models['Model-1'].Material(name='Aluminium')
mdb.models['Model-1'].materials['Aluminium'].Elastic(table=((70000.0, 0.3), ))

# Create Elastomer
mdb.models['Model-1'].Material(name='Elastomer')
mdb.models['Model-1'].materials['Elastomer'].Elastic(table=((3.0, 0.4999999), ))

#Create/Assign Section
print 'Creating the Sections'
mdb.models['Model-1'].HomogeneousShellSection(name='Solid', 
	preIntegrate=OFF, material='Aluminium', thicknessType=UNIFORM, 
	thickness=1.0, thicknessField='', idealization=NO_IDEALIZATION, 
	poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
	useDensity=OFF, integrationRule=SIMPSON, numIntPts=5)

#Create/Assign Section
print 'Creating the Sections'
mdb.models['Model-1'].HomogeneousShellSection(name='Stretchable', 
	preIntegrate=OFF, material='Elastomer', thicknessType=UNIFORM, 
	thickness=elast_thk, thicknessField='', idealization=NO_IDEALIZATION, 
	poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
	useDensity=OFF, integrationRule=SIMPSON, numIntPts=5)
	
print 'Assigning the Sections'

p = mdb.models['Model-1'].parts['Plate']
f = p.faces
faces = f.getByBoundingBox(-H,-V,-G,H,V,G)
region = p.Set(faces=faces, name='Set-1')
p = mdb.models['Model-1'].parts['Plate']
p.SectionAssignment(region=region, sectionName='Solid', offset=0.0, 
	offsetType=TOP_SURFACE, offsetField='', 
	thicknessAssignment=FROM_SECTION)

p = mdb.models['Model-1'].parts['Elastomer1']
f = p.faces
faces = f.getByBoundingBox(-H,-V,-G,H,V,G)
region = p.Set(faces=faces, name='Set-1')
p = mdb.models['Model-1'].parts['Elastomer1']
p.SectionAssignment(region=region, sectionName='Stretchable', offset=0.0, 
	offsetType=TOP_SURFACE, offsetField='', 
	thicknessAssignment=FROM_SECTION)
	
p = mdb.models['Model-1'].parts['Elastomer2L']
f = p.faces
faces = f.getByBoundingBox(-H,-V,-G,H,V,G)
region = p.Set(faces=faces, name='Set-1')
p = mdb.models['Model-1'].parts['Elastomer2L']
p.SectionAssignment(region=region, sectionName='Stretchable', offset=0.0, 
	offsetType=TOP_SURFACE, offsetField='', 
	thicknessAssignment=FROM_SECTION)

p = mdb.models['Model-1'].parts['Elastomer2R']
f = p.faces
faces = f.getByBoundingBox(-H,-V,-G,H,V,G)
region = p.Set(faces=faces, name='Set-1')
p = mdb.models['Model-1'].parts['Elastomer2R']
p.SectionAssignment(region=region, sectionName='Stretchable', offset=0.0, 
	offsetType=TOP_SURFACE, offsetField='', 
	thicknessAssignment=FROM_SECTION)

p = mdb.models['Model-1'].parts['Elastomer3L']
f = p.faces
faces = f.getByBoundingBox(-H,-V,-G,H,V,G)
region = p.Set(faces=faces, name='Set-1')
p = mdb.models['Model-1'].parts['Elastomer3L']
p.SectionAssignment(region=region, sectionName='Stretchable', offset=0.0, 
	offsetType=TOP_SURFACE, offsetField='', 
	thicknessAssignment=FROM_SECTION)

p = mdb.models['Model-1'].parts['Elastomer3R']
f = p.faces
faces = f.getByBoundingBox(-H,-V,-G,H,V,G)
region = p.Set(faces=faces, name='Set-1')
p = mdb.models['Model-1'].parts['Elastomer3R']
p.SectionAssignment(region=region, sectionName='Stretchable', offset=0.0, 
	offsetType=TOP_SURFACE, offsetField='', 
	thicknessAssignment=FROM_SECTION)
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

a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Right_Plate_2', ), axisPoint=((1.5*H+1.5*G),(V-D), HNG), 
	axisDirection=(0.0, (2*V), 0.0), angle=(phiR*180/pi))

a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Left_Plate_2', ), axisPoint=(-(1.5*H+1.5*G),-(V-D), HNG), 
	axisDirection=(0.0, -(2*V), 0.0), angle=(phiL*180/pi))
	
#Define Steps
print 'Defining the Steps'
mdb.models['Model-1'].StaticStep(name='RBM', previous='Initial')
mdb.models['Model-1'].steps['RBM'].setValues(nlgeom=ON)
mdb.models['Model-1'].steps['RBM'].setValues(solutionTechnique=FULL_NEWTON)

#Create Reference Points
print 'Defining Reference Points'

a = mdb.models['Model-1'].rootAssembly			# RPs for Middle to Right 1
RP=a.ReferencePoint(point=(0.5*(H+G),0.5*V-D,HNG))															#8,9
RP1=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(0.5*(H+G),0.5*V-D, HNG))
RP2=RP.id
a = mdb.models['Model-1'].rootAssembly			# RPs for Middle to Right 2
RP=a.ReferencePoint(point=(0.5*(H+G), -(0.5*V-D),HNG))													#10,11
RP3=RP.id															
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(0.5*(H+G), -(0.5*V-D), HNG))
RP4=RP.id

	
a = mdb.models['Model-1'].rootAssembly				# RPs for Middle to Left 1
RP=a.ReferencePoint(point=(-0.5*(H+G),0.5*V-D, HNG))															#12,13
RP5=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(-0.5*(H+G),0.5*V-D, HNG))
RP6=RP.id
a = mdb.models['Model-1'].rootAssembly				# RPs for Middle to Left 2
RP=a.ReferencePoint(point=(-0.5*(H+G),-(0.5*V-D), HNG))															#14,15
RP7=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(-0.5*(H+G),-(0.5*V-D), HNG))
RP8=RP.id


a = mdb.models['Model-1'].rootAssembly			# RPs for RHS Support
RP=a.ReferencePoint(point=(3*0.5*(H+G),(0.5*V-D),HNG))															#16,17
RP9=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(3*0.5*(H+G),(0.5*V-D), HNG))
RP10=RP.id
a = mdb.models['Model-1'].rootAssembly			# RPs for RHS Support 2
RP=a.ReferencePoint(point=(3*0.5*(H+G),-(0.5*V-D),HNG))															#18,19
RP11=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(3*0.5*(H+G),-(0.5*V-D), HNG))
RP12=RP.id


a = mdb.models['Model-1'].rootAssembly			# RPs for LHS Support
RP=a.ReferencePoint(point=(-3*0.5*(H+G),(0.5*V-D),HNG))															#RP13,21
RP13=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(-3*0.5*(H+G),(0.5*V-D), HNG))
RP14=RP.id
a = mdb.models['Model-1'].rootAssembly			# RPs for LHS Support 2
RP=a.ReferencePoint(point=(-3*0.5*(H+G),-(0.5*V-D),HNG))															#22,23
RP15=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(-3*0.5*(H+G),-(0.5*V-D), HNG))
RP16=RP.id

a = mdb.models['Model-1'].rootAssembly				# RPs for Rigid Body Constraints								#24,25,26
RP=a.ReferencePoint(point=((H+G), 0.0, 0.0))
RP17=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(-(H+G), 0.0, 0.0))
RP18=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(0.0, 0.0, 0.0))
RP19=RP.id
a = mdb.models['Model-1'].rootAssembly																	# RPs for Hinge Support at endswith								#27,28,29,30
RP=a.ReferencePoint(point=((1.5*H+1.5*G+(H+G)*cos(phiR)), (0.5*V-D), (-(H+G)*sin(phiR)+HNG)))
RP20=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=((1.5*H+1.5*G+(H+G)*cos(phiR)), -(0.5*V-D), (-(H+G)*sin(phiR)+HNG)))
RP21=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(-(1.5*H+1.5*G+(H+G)*cos(phiL)), (0.5*V-D), (-(H+G)*sin(phiL)+HNG)))
RP22=RP.id
a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(-(1.5*H+1.5*G+(H+G)*cos(phiL)), -(0.5*V-D), (-(H+G)*sin(phiL)+HNG)))
RP23=RP.id

a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=((1.5*H+1.5*G+0.5*(H+G)*cos(phiR)), 0, (-0.5*(H+G)*sin(phiR))))
RP24=RP.id

a = mdb.models['Model-1'].rootAssembly
RP=a.ReferencePoint(point=(-(1.5*H+1.5*G+0.5*(H+G)*cos(phiL)), 0, (-0.5*(H+G)*sin(phiL))))
RP25=RP.id

#Create Wires
print 'Defining Wires'

a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
v1 = a.instances['Middle_Plate'].vertices
v2 = a.instances['Right_Plate'].vertices
v3 = a.instances['Left_Plate'].vertices
v4 = a.instances['Right_Plate_2'].vertices
v5 = a.instances['Left_Plate_2'].vertices

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP1], v1.findAt(coordinates=((0.5*H-D),(0.5*V-D), 0.0))), ), mergeType=IMPRINT, meshable=OFF)

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP2], v2.findAt(coordinates=((0.5*H+D+G), (0.5*V-D), 0.0))), ), mergeType=IMPRINT, meshable=OFF)

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP3], v1.findAt(coordinates=((0.5*H-D), -(0.5*V-D), 0.0))), ), mergeType=IMPRINT, meshable=OFF)

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP4], v2.findAt(coordinates=((0.5*H+D+G), -(0.5*V-D), 0.0))), ), mergeType=IMPRINT, meshable=OFF)

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP5], v1.findAt(coordinates=(-(0.5*H-D), (0.5*V-D), 0.0))), ), mergeType=IMPRINT, meshable=OFF)

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP7], v1.findAt(coordinates=(-(0.5*H-D), -(0.5*V-D), 0.0))), ), mergeType=IMPRINT, meshable=OFF)

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP8], v3.findAt(coordinates=(-(0.5*H+D+G), -(0.5*V-D), 0.0))), ), mergeType=IMPRINT, meshable=OFF)

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP6], v3.findAt(coordinates=(-(0.5*H+D+G), (0.5*V-D), 0.0))), ), mergeType=IMPRINT, meshable=OFF)

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v3.findAt(coordinates=(-(1.5*H+G-D), -(0.5*V-D), 0.0)), r1[RP16]), ), mergeType=IMPRINT, meshable=OFF)			

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v3.findAt(coordinates=(-(1.5*H+G-D), (0.5*V-D), 0.0)), r1[RP14]), ), mergeType=IMPRINT, meshable=OFF)				

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v2.findAt(coordinates=((1.5*H+G-D), (0.5*V-D), 0.0)), r1[RP10]), ), mergeType=IMPRINT, meshable=OFF)				

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((v2.findAt(coordinates=((1.5*H+G-D), -(0.5*V-D), 0.0)), r1[RP12]), ), mergeType=IMPRINT, meshable=OFF)				

phi0=atan2((0.5*G+D),HNG)
phi1=phi0-phiL
phi2=phi0-phiR
LHNG=sqrt((0.5*G+D)*(0.5*G+D)+(HNG*HNG))
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP9], v4.findAt(coordinates=((1.5*H+1.5*G+LHNG*sin(phi2)), (0.5*V-D),(HNG-LHNG*cos(phi2))))), ), mergeType=IMPRINT, meshable=OFF)	
a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP11], v4.findAt(coordinates=((1.5*H+1.5*G+LHNG*sin(phi2)), -(0.5*V-D),(HNG-LHNG*cos(phi2))))), ), mergeType=IMPRINT, meshable=OFF)	


a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP13], v5.findAt(coordinates=(-(1.5*H+1.5*G+LHNG*sin(phi1)), (0.5*V-D),(HNG-LHNG*cos(phi1))))), ), mergeType=IMPRINT, meshable=OFF)		

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP15], v5.findAt(coordinates=(-(1.5*H+1.5*G+LHNG*sin(phi1)), -(0.5*V-D),(HNG-LHNG*cos(phi1))))), ), mergeType=IMPRINT, meshable=OFF)		

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP20], v4.findAt(coordinates=((1.5*H+1.5*G+LHNG*sin(phi2)+(H-2*D)*cos(phiR)), (0.5*V-D),(HNG-LHNG*cos(phi2)-(H-2*D)*sin(phiR))))), ), mergeType=IMPRINT, meshable=OFF)

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP21], v4.findAt(coordinates=((1.5*H+1.5*G+LHNG*sin(phi2)+(H-2*D)*cos(phiR)), -(0.5*V-D),(HNG-LHNG*cos(phi2)-(H-2*D)*sin(phiR))))), ), mergeType=IMPRINT, meshable=OFF)	

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP22], v5.findAt(coordinates=(-(1.5*H+1.5*G+LHNG*sin(phi1)+(H-2*D)*cos(phiL)), (0.5*V-D),(HNG-LHNG*cos(phi1)-(H-2*D)*sin(phiL))))), ), mergeType=IMPRINT, meshable=OFF)

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r1[RP23], v5.findAt(coordinates=(-(1.5*H+1.5*G+LHNG*sin(phi1)+(H-2*D)*cos(phiL)), -(0.5*V-D),(HNG-LHNG*cos(phi1)-(H-2*D)*sin(phiL))))), ), mergeType=IMPRINT, meshable=OFF)

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
e1 = a.edges
edges1 = e1.findAt(((0.5*(H+G),0.5*V-D,HNG), ))
a.Set(edges=edges1, name='Wire-1')

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP3], r11[RP4]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((0.5*(H+G), -(0.5*V-D),HNG), ))
a.Set(edges=edges1, name='Wire-2')

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP5], r11[RP6]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((-0.5*(H+G),0.5*V-D, HNG), ))
a.Set(edges=edges1, name='Wire-3')

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP7], r11[RP8]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((-0.5*(H+G),-(0.5*V-D), HNG), ))
a.Set(edges=edges1, name='Wire-4')

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP9], r11[RP10]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((3*0.5*(H+G),(0.5*V-D),HNG), ))
a.Set(edges=edges1, name='Wire-5')

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP11], r11[RP12]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((3*0.5*(H+G),-(0.5*V-D),HNG), ))
a.Set(edges=edges1, name='Wire-6')

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP13], r11[RP14]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((-3*0.5*(H+G),(0.5*V-D),HNG), ))
a.Set(edges=edges1, name='Wire-7')

a = mdb.models['Model-1'].rootAssembly
r11 = a.referencePoints
a.WirePolyLine(points=((r11[RP15], r11[RP16]), ), mergeType=IMPRINT, meshable=OFF)
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((-3*0.5*(H+G),-(0.5*V-D),HNG), ))
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
a.resumeFeatures(('Wire-1', 'Wire-2', 'Wire-3', 'Wire-4', 'Wire-5', 'Wire-6', 
	'Wire-7', 'Wire-8', 'Wire-9', 'Wire-10','Wire-11','Wire-12','Wire-13','Wire-14','Wire-15','Wire-16','Wire-17','Wire-18','Wire-19','Wire-20', ))

#Create and Assign Connector Sections
print 'Create Connector Section'	
mdb.models['Model-1'].ConnectorSection(name='HINGE', assembledType=HINGE)

print 'Create Connector Section'	
mdb.models['Model-1'].ConnectorSection(name='BEAM', assembledType=BEAM)

print 'Define Datum Coordinate Systems'
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_MR', coordSysType=CARTESIAN, origin=(0.5*(H+G),0.5*V-D,HNG), point1=(0.5*(H+G),-0.5*V,HNG), point2=(0.0, 0.0, HNG))			#CSYS for Middle to Right 1
	
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_MR2', coordSysType=CARTESIAN, origin=(0.5*(H+G), -(0.5*V-D),HNG), point1=(0.5*(H+G),-0.5*V,HNG), point2=(0.0, 0.0, HNG))			#CSYS for Middle to R2

a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_ML', coordSysType=CARTESIAN, origin=(-0.5*(H+G),0.5*V-D, HNG), point1=(-0.5*(H+G),0.5*V, HNG), point2=(0.0, 0.0, HNG))		#CSYS for Middle to Left

a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_ML2', coordSysType=CARTESIAN, origin=(-0.5*(H+G),-(0.5*V-D), HNG), point1=(-0.5*(H+G),HNG, HNG), point2=(0.0, 0.0, HNG))		#CSYS for Middle to L2

a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_RS', coordSysType=CARTESIAN, origin=(3*0.5*(H+G),(0.5*V-D),HNG), point1=(3*0.5*(H+G),HNG, HNG), point2=(0.0, 0.0,HNG))			#CSYS for  RHS Support
	
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_RS2', coordSysType=CARTESIAN, origin=(3*0.5*(H+G),-(0.5*V-D),HNG), point1=(3*0.5*(H+G),HNG,HNG), point2=(0.0, 0.0,HNG))			#CSYS for RHS Support 2

a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_LS', coordSysType=CARTESIAN, origin=(-3*0.5*(H+G),(0.5*V-D), HNG), point1=(-3*0.5*(H+G),HNG,HNG), point2=(0.0, 0.0, HNG))		#CSYS for LHS Support

a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByThreePoints(name='CSYS_LS2', coordSysType=CARTESIAN, origin=(-3*0.5*(H+G),-(0.5*V-D),HNG), point1=(-3*0.5*(H+G),HNG,HNG), point2=(0.0, 0.0,HNG))		#CSYS for LHS Support 2


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

#Add Elastomer

a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Elastomer1']
a.Instance(name='Elastomer_R1', part=p, dependent=ON)

a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Elastomer1']
a.Instance(name='Elastomer_L1', part=p, dependent=ON)

a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Elastomer2R']
a.Instance(name='Elastomer_R2', part=p, dependent=ON)

a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Elastomer2L']
a.Instance(name='Elastomer_L2', part=p, dependent=ON)

a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Elastomer3R']
a.Instance(name='Elastomer_R3', part=p, dependent=ON)

a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Elastomer3L']
a.Instance(name='Elastomer_L3', part=p, dependent=ON)

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Middle_Plate'].edges
edges1 = e1.findAt((((0.5*H), (0.25*V), 0.0), ), (((0.5*H), -(0.25*V), 0.0), ))
region1=a.Surface(side1Edges=edges1, name='MR')

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Elastomer_R1'].edges
edges1 = e1.findAt((((-0.5*G), 0, 0.0), ))
region2=a.Surface(side1Edges=edges1, name='ER1L')
mdb.models['Model-1'].Tie(name='ER1_1', master=region1, slave=region2, 
	positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Right_Plate'].edges
edges1 = e1.findAt((((0.5*H+G), (0.25*V), 0.0), ), (((0.5*H+G), -(0.25*V), 0.0), ))
region1=a.Surface(side1Edges=edges1, name='RL')

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Elastomer_R1'].edges
edges1 = e1.findAt((((0.5*G), 0, 0.0), ))
region2=a.Surface(side1Edges=edges1, name='ER1R')
mdb.models['Model-1'].Tie(name='ER1_2', master=region1, slave=region2, 
	positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
	
a1 = mdb.models['Model-1'].rootAssembly
a1.translate(instanceList=('Elastomer_R1', ), vector=((0.5*H+0.5*G), 0.0, 0.0))

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Middle_Plate'].edges
edges1 = e1.findAt(((-(0.5*H), (0.25*V), 0.0), ), ((-(0.5*H), -(0.25*V), 0.0), ))
region1=a.Surface(side1Edges=edges1, name='ML')

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Elastomer_L1'].edges
edges1 = e1.findAt((((0.5*G), 0, 0.0), ))
region2=a.Surface(side1Edges=edges1, name='EL1R')
mdb.models['Model-1'].Tie(name='EL1_1', master=region1, slave=region2, 
	positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Left_Plate'].edges
edges1 = e1.findAt(((-(0.5*H+G), (0.25*V), 0.0), ), ((-(0.5*H+G), -(0.25*V), 0.0), ))
region1=a.Surface(side1Edges=edges1, name='LR')

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Elastomer_L1'].edges
edges1 = e1.findAt(((-(0.5*G), 0, 0.0), ))
region2=a.Surface(side1Edges=edges1, name='EL1L')
mdb.models['Model-1'].Tie(name='EL1_2', master=region1, slave=region2, 
	positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

a1 = mdb.models['Model-1'].rootAssembly
a1.translate(instanceList=('Elastomer_L1', ), vector=(-(0.5*H+0.5*G), 0.0, 0.0))

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Right_Plate'].edges
edges1 = e1.findAt((((1.5*H+G), (0.25*V), 0.0), ), (((1.5*H+G), -(0.25*V), 0.0), ))
region1=a.Surface(side1Edges=edges1, name='RR')

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Elastomer_R2'].edges
edges1 = e1.findAt((((-0.5*LD2), 0, 0.0), ))
region2=a.Surface(side1Edges=edges1, name='ER2L')
mdb.models['Model-1'].Tie(name='ER2_1', master=region1, slave=region2, 
	positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Right_Plate_2'].edges
edges1 = e1.findAt((((1.5*H+1.5*G+LHNG*sin(phi2)-D*cos(phiR)), (0.25*V),(HNG-LHNG*cos(phi2)+D*sin(phiR))), ), (((1.5*H+1.5*G+LHNG*sin(phi2)-D*cos(phiR)), -(0.25*V),(HNG-LHNG*cos(phi2)+D*sin(phiR))), ))
region1=a.Surface(side1Edges=edges1, name='R2L')

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Elastomer_R2'].edges
edges1 = e1.findAt((((0.5*LD2), 0, 0.0), ))
region2=a.Surface(side1Edges=edges1, name='ER2R')
mdb.models['Model-1'].Tie(name='ER2_2', master=region1, slave=region2, 
	positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

a1 = mdb.models['Model-1'].rootAssembly
a1.translate(instanceList=('Elastomer_R2', ), vector=((1.5*H+G+0.5*LD2), 0.0, 0.0))	

a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Elastomer_R2', ), axisPoint=((1.5*H+G),(-V),0), 
	axisDirection=(0.0, (2*V), 0.0), angle=(phiR*90/pi))	

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Left_Plate'].edges
edges1 = e1.findAt(((-(1.5*H+G), (0.25*V), 0.0), ), ((-(1.5*H+G), -(0.25*V), 0.0), ))
region1=a.Surface(side1Edges=edges1, name='LL')

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Elastomer_L2'].edges
edges1 = e1.findAt((((0.5*LD1), 0, 0.0), ))
region2=a.Surface(side1Edges=edges1, name='EL2R')
mdb.models['Model-1'].Tie(name='EL2_1', master=region1, slave=region2, 
	positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Left_Plate_2'].edges
edges1 = e1.findAt(((-(1.5*H+1.5*G+LHNG*sin(phi1)-D*cos(phiL)), (0.25*V),(HNG-LHNG*cos(phi1)+D*sin(phiL))), ), ((-(1.5*H+1.5*G+LHNG*sin(phi1)-D*cos(phiL)), -(0.25*V),(HNG-LHNG*cos(phi1)+D*sin(phiL))), ))
region1=a.Surface(side1Edges=edges1, name='L2R')

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Elastomer_L2'].edges
edges1 = e1.findAt(((-(0.5*LD1), 0, 0.0), ))
region2=a.Surface(side1Edges=edges1, name='EL2L')
mdb.models['Model-1'].Tie(name='EL2_2', master=region1, slave=region2, 
	positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)	
	
a1 = mdb.models['Model-1'].rootAssembly
a1.translate(instanceList=('Elastomer_L2', ), vector=(-(1.5*H+G+0.5*LD1), 0.0, 0.0))		
	
a = mdb.models['Model-1'].rootAssembly
a.rotate(instanceList=('Elastomer_L2', ), axisPoint=(-(1.5*H+G),(V),0), 
	axisDirection=(0.0, (-2*V), 0.0), angle=(phiL*90/pi))

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Left_Plate_2'].edges
edges1 = e1.findAt(((-(1.5*H+1.5*G+LHNG*sin(phi1)+(H-D)*cos(phiL)), (0.25*V),(HNG-LHNG*cos(phi1)-(H-D)*sin(phiL))), ), ((-(1.5*H+1.5*G+LHNG*sin(phi1)+(H-D)*cos(phiL)), -(0.25*V),(HNG-LHNG*cos(phi1)-(H-D)*sin(phiL))), ))
region1=a.Surface(side1Edges=edges1, name='L2L')

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Elastomer_L3'].edges
edges1 = e1.findAt((((0.5*LD3), 0, 0.0), ))
region2=a.Surface(side1Edges=edges1, name='EL3R')
mdb.models['Model-1'].Tie(name='EL3', master=region1, slave=region2, 
	positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)	

a1 = mdb.models['Model-1'].rootAssembly
a1.translate(instanceList=('Elastomer_L3', ), vector=(-(0.5*LD3+1.5*H+1.5*G+LHNG2*sin(phi02-phiL)+(H)*cos(phiL)), 0,(HNG-LHNG*cos(phi1)-(H-D)*sin(phiL))))
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Right_Plate_2'].edges
edges1 = e1.findAt((((1.5*H+1.5*G+LHNG*sin(phi2)+(H-D)*cos(phiR)), (0.25*V),(HNG-LHNG*cos(phi2)-(H-D)*sin(phiR))), ), (((1.5*H+1.5*G+LHNG*sin(phi2)+(H-D)*cos(phiR)), -(0.25*V),(HNG-LHNG*cos(phi2)-(H-D)*sin(phiR))), ))
region1=a.Surface(side1Edges=edges1, name='R2R')
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Elastomer_R3'].edges
edges1 = e1.findAt(((-(0.5*LD4), 0, 0.0), ))
region2=a.Surface(side1Edges=edges1, name='ER3L')
mdb.models['Model-1'].Tie(name='ER3', master=region1, slave=region2, 
	positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)	
a1 = mdb.models['Model-1'].rootAssembly
a1.translate(instanceList=('Elastomer_R3', ), vector=((0.5*LD4+1.5*H+1.5*G+LHNG2*sin(phi02-phiR)+(H)*cos(phiR)), 0,(HNG-LHNG*cos(phi2)-(H-D)*sin(phiR))))

#Create Rigid Body Constraints
print 'Defining RBC'
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

#Define BCs				
print 'Defining all BCs'

a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[RP20], r1[RP21], r1[RP22], r1[RP23], )
region = regionToolset.Region(referencePoints=refPoints1)
mdb.models['Model-1'].DisplacementBC(name='Fixed hinge Support', createStepName='Initial', 
	region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
	amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
	localCsys=None)

a = mdb.models['Model-1'].rootAssembly
region=a.sets['MRH']
mdb.models['Model-1'].ConnDisplacementBC(name='Theta_R', createStepName='RBM', 
	region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=theta_r, ur2=UNSET, ur3=UNSET, 
	amplitude=UNSET, distributionType=UNIFORM)


a = mdb.models['Model-1'].rootAssembly
region=a.sets['MLH']
mdb.models['Model-1'].ConnDisplacementBC(name='Theta_L', createStepName='RBM', 
	region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=theta_l, ur2=UNSET, ur3=UNSET, 
	amplitude=UNSET, distributionType=UNIFORM)

a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Elastomer_L3'].edges 
edges1 = e1.findAt(((-(LD3+1.5*H+1.5*G+LHNG2*sin(phi02-phiL)+(H)*cos(phiL)), 0,(HNG-LHNG*cos(phi1)-(H-D)*sin(phiL))), ))
e2 = a.instances['Elastomer_R3'].edges
edges2 = e2.findAt((((LD4+1.5*H+1.5*G+LHNG2*sin(phi02-phiR)+(H)*cos(phiR)), 0,(HNG-LHNG*cos(phi2)-(H-D)*sin(phiR))), ))
region = a.Set(edges=edges1+edges2, name='Set-43')
mdb.models['Model-1'].EncastreBC(name='Fix_Elast', createStepName='RBM', 
	region=region, localCsys=None)


#Define Sets
print 'Defining Sets'

a = mdb.models['Model-1'].rootAssembly
v5 = a.instances['Left_Plate_2'].vertices
verts5 = v5.findAt(((-(1.5*H+1.5*G+LHNG*sin(phi1)+(H-D)*cos(phiL)), 0.0,(HNG-LHNG*cos(phi1)-(H-D)*sin(phiL))), ), ((-(1.5*H+1.5*G+LHNG*sin(phi1)-D*cos(phiL)), 0,(HNG-LHNG*cos(phi1)+D*sin(phiL))), ))
v1 = a.instances['Left_Plate'].vertices
verts1 = v1.findAt(((-(1.5*H+G),0,0), ), ((-(0.5*H+G),0,0), ))
v3 = a.instances['Middle_Plate'].vertices
verts3 = v3.findAt(((-(0.5*H),0,0), ), ((0.5*H,0,00), ))
v4 = a.instances['Right_Plate'].vertices
verts4 = v4.findAt(((+(0.5*H+G),0,0), ), (((1.5*H+G),0,0), ))
a = mdb.models['Model-1'].rootAssembly
v6 = a.instances['Right_Plate_2'].vertices
verts6 = v6.findAt((((1.5*H+1.5*G+LHNG*sin(phi2)+(H-D)*cos(phiR)), 0.0,(HNG-LHNG*cos(phi2)-(H-D)*sin(phiR))), ), (((1.5*H+1.5*G+LHNG*sin(phi2)-D*cos(phiR)), 0,(HNG-LHNG*cos(phi2)+D*sin(phiR))), ))
a.Set(vertices=verts5+verts1+verts3+verts4+verts6, name='TIPNODE')


#Mesh Parts
print 'Meshing the Part'
p = mdb.models['Model-1'].parts['Plate']
p.seedPart(size=V/8, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['Plate']
p.generateMesh()

elemType1 = mesh.ElemType(elemCode=S8R, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=STRI65, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['Elastomer1']
f = p.faces
faces = f.findAt(((0, 0, 0.0), ))
pickedRegions =(faces, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
p.seedPart(size=0.11*G, deviationFactor=0.1, minSizeFactor=0.1)
p.generateMesh()

p = mdb.models['Model-1'].parts['Elastomer2R']
f = p.faces
faces = f.findAt(((0, 0, 0.0), ))
pickedRegions =(faces, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
e = p.edges
pickedEdges = e.findAt(((0, 0.5*V, 0.0), ), ((0, -(0.5*V), 0.0), ))
p.seedEdgeByNumber(edges=pickedEdges, number=12, constraint=FINER)
p.generateMesh()

p = mdb.models['Model-1'].parts['Elastomer2L']
f = p.faces
faces = f.findAt(((0, 0, 0.0), ))
pickedRegions =(faces, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
e = p.edges
pickedEdges = e.findAt(((0, 0.5*V, 0.0), ), ((0, -(0.5*V), 0.0), ))
p.seedEdgeByNumber(edges=pickedEdges, number=12, constraint=FINER)
p.generateMesh()

p = mdb.models['Model-1'].parts['Elastomer3L']
f = p.faces
faces = f.findAt(((0, 0, 0.0), ))
pickedRegions =(faces, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
p.seedPart(size=0.11*G, deviationFactor=0.1, minSizeFactor=0.1)
p.generateMesh()

p = mdb.models['Model-1'].parts['Elastomer3R']
f = p.faces
faces = f.findAt(((0, 0, 0.0), ))
pickedRegions =(faces, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
p.seedPart(size=0.11*G, deviationFactor=0.1, minSizeFactor=0.1)
p.generateMesh()

a = mdb.models['Model-1'].rootAssembly
a.regenerate()
e6 = a.instances['Elastomer_R1'].elements
elements6 = e6[:]
e7 = a.instances['Elastomer_L1'].elements
elements7 = e7[:]
e8 = a.instances['Elastomer_R2'].elements
elements8 = e8[:]
e9 = a.instances['Elastomer_L2'].elements
elements9 = e9[:]
e10 = a.instances['Elastomer_R3'].elements
elements10 = e10[:]
e11 = a.instances['Elastomer_L3'].elements
elements11 = e11[:]
a.Set(elements=elements6+elements7+elements8+elements9+elements10+elements11, name='ALL_PART')

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
job.submit()
job.waitForCompletion()
print 'Completed job'
tipDisp=[0 for lenn in range(27)]
tipDisp = getResults(ModelName)
phi02=atan2((0.5*G),HNG)
phi1=phi02-phiL
phi2=phi02-phiR
LHNG2=sqrt((0.5*G)*(0.5*G)+(HNG*HNG))
Ax=(0.5*H)
Az=0
Bx=(0.5*H+G)
Bz=0
Cx=(1.5*H+G)
Cz=0
Dx1=(1.5*H+1.5*G+LHNG2*sin(phi1))
Dz1=(HNG-LHNG2*cos(phi1))
Dx2=(1.5*H+1.5*G+LHNG2*sin(phi2))
Dz2=(HNG-LHNG2*cos(phi2))
Ex1=(1.5*H+1.5*G+LHNG2*sin(phi1)+(H)*cos(phiL))
Ez1=(HNG-LHNG2*cos(phi1)-(H)*sin(phiL))
Ex2=(1.5*H+1.5*G+LHNG2*sin(phi2)+(H)*cos(phiR))
Ez2=(HNG-LHNG2*cos(phi2)-(H)*sin(phiR))
A = ((-Ex1+tipDisp[2]),(Ez1+tipDisp[12]),(-Dx1+tipDisp[3]),(Dz1+tipDisp[13]),(-Cx+tipDisp[0]),(Cz+tipDisp[10]),(-Bx+tipDisp[1]),(Bz+tipDisp[11]),(-Ax+tipDisp[8]),(Az+tipDisp[18]),(Ax+tipDisp[9]),(Az+tipDisp[19]),(Bx+tipDisp[6]),(Bz+tipDisp[16]),(Cx+tipDisp[7]),(Cz+tipDisp[17]),(Dx2+tipDisp[4]),(Dz2+tipDisp[14]),(Ex2+tipDisp[5]),(Ez2+tipDisp[15]))
x=[0 for ind in range(11)]
x[0]=-0.5*Total
x[1:10]=A[0:19:2]
x[11]=0.5*Total
z=[0 for ind in range(11)]
z[0]=Ez1
z[1:10]=A[1:20:2]
z[11]=Ez2
fileobject = open('temp.txt','wb')
for xloc in x:
	fileobject.write('%.4f\n' % xloc)
for zloc in z:
	fileobject.write('%.4f\n' % zloc)
fileobject.close()

#tipDisp =[9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999,9999]
if os.access('%s.lck'%ModelName,os.F_OK):
	os.remove('%s.lck'%ModelName)
DataFile = open('PostData.txt','a')
#DataFile.write('\n')
DataFile.write('%10f,%10f,%10f,%10f,%10f,%10f,' % (DVs[0],DVs[1],DVs[2],DVs[3],DVs[4],DVs[5], ))
DataFile.write('%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f,' % ((-Ex1+tipDisp[2]),(Ez1+tipDisp[12]),(-Dx1+tipDisp[3]),(Dz1+tipDisp[13]),(-Cx+tipDisp[0]),(Cz+tipDisp[10]),(-Bx+tipDisp[1]),(Bz+tipDisp[11]),(-Ax+tipDisp[8]),(Az+tipDisp[18]),(Ax+tipDisp[9]),(Az+tipDisp[19]),(Bx+tipDisp[6]),(Bz+tipDisp[16]),(Cx+tipDisp[7]),(Cz+tipDisp[17]),(Dx2+tipDisp[4]),(Dz2+tipDisp[14]),(Ex2+tipDisp[5]),(Ez2+tipDisp[15]),tipDisp[20], ))
DataFile.close()
Mdb()
print 'DONE!!'
outputs = tipDisp
# Output data
# This creates the input txt file for the Abaqus script
# fileobject = open('output.txt','wb')
# for output in outputs:
	# fileobject.write('%.4f\n' % output)
# fileobject.close()