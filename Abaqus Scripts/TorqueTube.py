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
H=550; # Horizontal plate length 
V=125; # Vertical plate length
G=40;   # Gap length
HNG=30; # Hinge height
D=15; #Hinge distance from edge
phiL=0.2; #initial angle
phiR=0.2;
theta_m=0.0;
theta_l=-0.4;
theta_l2=0.4;
elast_thk=5
Total = 1000
seedsize=8
theta=asin((sin(phiR)-sin(phiL))/3)
Hnet=(Total)/(3*cos(theta)+cos(phiL)+cos(phiR));
H=Hnet-G;
##########################


### Write data file column headings

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
# Variables
ORTT=12.5
IRTT=6.25
V=125
D=15
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
## Create Section##
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
## Reference points##
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Torque_Tube']
a.Instance(name='Torque_Tube-1', part=p, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
e11 = a.instances['Torque_Tube-1'].edges
RP=a.ReferencePoint(point=(0,0,-(0.5*V-D)))
RPTT1=RP.id
a = mdb.models['Model-1'].rootAssembly
e21 = a.instances['Torque_Tube-1'].edges
RP=a.ReferencePoint(point=(0,0,(0.5*V-D)))
RPTT2=RP.id

a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[RPTT2], )
a.Set(referencePoints=refPoints1, name='RP_Z+')

a = mdb.models['Model-1'].rootAssembly
r1 = a.referencePoints
refPoints1=(r1[RPTT1], )
a.Set(referencePoints=refPoints1, name='RP_Z-')
## Constraints##
a = mdb.models['Model-1'].rootAssembly
region1=a.sets['RP_Z-']
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Torque_Tube-1'].faces
faces1 = f1.findAt(((-0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ), ((-0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), 
    -(0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), 0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ), ((0.5*(IRTT+ORTT)*cos(pi/4), -0.5*(IRTT+ORTT)*cos(pi/4), -(0.5*V-D)), ))
region2=a.Set(faces=faces1, name='s_Set-3')
mdb.models['Model-1'].Coupling(name='Kinematic_Coupling', controlPoint=region1, 
    surface=region2, influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC, 
    localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

## BCs##
a = mdb.models['Model-1'].rootAssembly
region = a.sets['RP_Z-']
datum = mdb.models['Model-1'].rootAssembly.datums[1]
mdb.models['Model-1'].DisplacementBC(name='Fix_TT', createStepName='Initial', 
    region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
    amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=datum)
a = mdb.models['Model-1'].rootAssembly
region = a.instances['Torque_Tube-1'].sets['TT']
mdb.models['Model-1'].Temperature(name='Predefined Field-1', 
    createStepName='Initial', region=region, distributionType=UNIFORM, 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(280.0, ))

mdb.models['Model-1'].CoupledTempDisplacementStep(name='Temp-Disp', 
    previous='Initial', timePeriod=20.0, maxNumInc=1000000, initialInc=0.1, 
    minInc=1e-08, maxInc=0.1, deltmx=100.0)

mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'E', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP', 
    'NT', 'HFL', 'RFL', 'HBF', 'MISES', 'SDV', 'TEMP'))
mdb.models['Model-1'].TabularAmplitude(name='Amp-1', timeSpan=STEP, 
    smooth=SOLVER_DEFAULT, data=((0.0, 280.0), (20.0, 400.0)))

a = mdb.models['Model-1'].rootAssembly
region = a.instances['Torque_Tube-1'].sets['TT']
mdb.models['Model-1'].TemperatureBC(name='Temp_change', createStepName='Temp-Disp', 
    region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    magnitude=1.0, amplitude='Amp-1')
## Job##
mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, 
    userSubroutine='C:\\Temp\\Torque Tube Panel- Raoul\\Working Database\\UMAT_3D_Coupled_ML_IP_Original.for', 
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
    numGPUs=0)

