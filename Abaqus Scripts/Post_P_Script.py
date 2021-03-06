"""
...
"""

from abaqus import *
from abaqusConstants import *
import visualization
from viewerModules import *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getResults(ModelName):
	disppTip=[0 for x in range(21)]
	"""
	This ODB reading script does the following:
	-Retrieves the displacement at TIPNODE
	-Scans for max. Mises stress in a part (if set exists)
	"""

	# Open the output database.
	print 'Made it in'
	odbName = ModelName+'.odb'
	print odbName
	odb = visualization.openOdb(odbName)
	lastFrame = odb.steps['RBM'].frames[-1]
	print 'odb open'

	# Selecting the node(s) to be queried
	pTip = odb.rootAssembly.nodeSets['TIPNODE']
		
	# Retrieve Y-displacements at the splines/connectors
	print 'Retrieving ALL final displacements at ALL points'
	dispField = lastFrame.fieldOutputs['U']

	print 'Retrieving ALL displacements at TIPNODE'
	dFieldpTip = dispField.getSubset(region=pTip)

	print 'Retrieving only U2 at TIPNODE'
	#Note, U1=data[0], U2=data[1], U3=data[2]

	for i in range(10):
		disppTip[i] = dFieldpTip.values[i].data[0]
		disppTip[i+10] = dFieldpTip.values[i].data[2]			
	#The following is left for use in later probems/projects
	print 'Scanning the PART for maximum VM STRESS'
	elsetName='ALL_PART'
	elset = elemset = None
	region = "over the entire model"
	assembly = odb.rootAssembly

	#Check to see if the element set exists
	#in the assembly

	if elsetName:
	   try:
		   elemset = assembly.elementSets[elsetName]
		   region = " in the element set : " + elsetName;
	   except KeyError:
		   print 'An assembly level elset named %s does' \
				  'not exist in the output database %s' \
				  % (elsetName, odbName)
		   odb.close()
		   exit(0)
		   
	""" Initialize maximum values """
	maxMises = -0.1
	maxVMElem = 0
	maxStep = "_None_"
	maxFrame = -1
	Stress = 'S'
	isStressPresent = 0
	for step in odb.steps.values():
	   print 'Processing Step:', step.name
	   for frame in step.frames:
		   allFields = frame.fieldOutputs
		   if (allFields.has_key(Stress)):
			   isStressPresent = 1
			   stressSet = allFields[Stress]
			   if elemset:
				   stressSet = stressSet.getSubset(
					   region=elemset)      
			   for stressValue in stressSet.values:                
				   if (stressValue.mises > maxMises):
					   maxMises = stressValue.mises
					   maxVMElem = stressValue.elementLabel
					   maxStep = step.name
					   maxFrame = frame.incrementNumber
	if(isStressPresent):
	   print 'Maximum von Mises stress %s is %f in element %d'%(
		   region, maxMises, maxVMElem)
	   print 'Location: frame # %d  step:  %s '%(maxFrame,maxStep)
	else:
	   print 'Stress output is not available in' \
			 'the output database : %s\n' %(odb.name)     
	disppTip[20]=maxMises
	odb.close()
	return disppTip

def getResults2(ModelName):
	disppTip=[[0 for y in range(3)] for x in range(10000)]
	disppTip2=[[0 for y in range(3)] for x in range(10000)]
	"""
	This ODB reading script does the following:
	-Retrieves the displacement at TIPNODE
	-Scans for max. Mises stress in a part (if set exists)
	"""

	# Open the output database.
	print 'Made it in'
	odbName = ModelName+'.odb'
	print odbName
	odb = visualization.openOdb(odbName)
	lastFrame = odb.steps['RBM'].frames[-1]
	firstFrame = odb.steps['RBM'].frames[0]
	print 'odb open'

	# Selecting the node(s) to be queried
	pTip = odb.rootAssembly.nodeSets['TIPNODE']
	dtm = odb.rootAssembly.datumCsyses['Datum_Aero_Data']
	session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(transformationType=USER_SPECIFIED, datumCsys=dtm)	
	# Retrieve Y-displacements at the splines/connectors
	print 'Retrieving ALL final displacements at ALL points'
	dispField = firstFrame.fieldOutputs['COORD']
	dispField2 = lastFrame.fieldOutputs['COORD']
	print 'Retrieving ALL displacements at TIPNODE'
	dFieldpTip = dispField.getSubset(region=pTip)
	dFieldpTip2 = dispField2.getSubset(region=pTip)
	print 'Retrieving only U2 at TIPNODE'
	#Note, U1=data[0], U2=data[1], U3=data[2]
	nnodes=len(dFieldpTip.values)
	for i in range(nnodes):
		disppTip[i][0] = dFieldpTip.values[i].data[0]
		disppTip[i][1] = dFieldpTip.values[i].data[1]
		disppTip[i][2] = dFieldpTip.values[i].data[2]			
		disppTip2[i][0] = dFieldpTip2.values[i].data[0]
		disppTip2[i][1] = dFieldpTip2.values[i].data[1]
		disppTip2[i][2] = dFieldpTip2.values[i].data[2]
		#The following is left for use in later probems/projects
	odb.close()
	disppTip=np.array(disppTip)
	disppTip2=np.array(disppTip2)
	disppTip=disppTip[~np.all(disppTip == 0, axis=1)]
	disppTip2=disppTip2[~np.all(disppTip2 == 0, axis=1)]
	np.savetxt("initial.txt", disppTip, delimiter=",")
	np.savetxt("final.txt", disppTip2, delimiter=",")
	return disppTip,disppTip2
	
def find_strain(JobName, StepName, FOTYPE, SetName = None, object = 'assembly', 

                object_name = None, strain_metric = 'mises'):

    """Finds Strain for a Set if SetName defined.Currently only works for mises.

        Otherwise finds maxMises for the whole model
        - object = 'assembly' or 'instance' (in that case specify name of instance.

        - strain_metric = 'mises' to get maximum equivalent stress and 'principal'"""

    odbName = JobName + '.odb'

    odb = visualization.openOdb(odbName)

    lastFrame = odb.steps[StepName].frames[-1]

    deformation = lastFrame.fieldOutputs[FOTYPE]

    # If not Set defined, get the value for the whole thing

    if SetName == None:

        raise NotImplementedError

    else:

        elsetName = SetName

        elset = elemset = None

        assembly = odb.rootAssembly

        #Check to see if the element set exists

        #in the assembly

        if elsetName:

           try:

               if object == 'assembly':

                   elemset = assembly.elementSets[elsetName]

                   region = " in the element set : " + elsetName;

               elif object == 'instance':

                   elemset = assembly.instances[object_name.upper()].elementSets[elsetName.upper()]                    

           except KeyError:

               print assembly.elementSets

               print 'An assembly level elset named %s does' \
                      'not exist in the output database %s' \
                      % (elsetName, odbName)

               exit(0)

        """ Initialize maximum values """

        maxStrain = -0.1

        maxStrainElem = 0

        maxStep = "_None_"

        maxFrame = -1

        Stress = 'S'

        isStressPresent = 0

        for step in odb.steps.values():

           if StepName == None or step.name == StepName:

               if elemset:

                   strainSet = deformation.getSubset(region=elemset)

               for strainValue in strainSet.values:

                   if strain_metric == 'mises':

                      strain_i = strainValue.mises

                   elif strain_metric == 'principal':

                      strain_i = abs(strainValue.maxPrincipal)
                      strain_act=strainValue.maxPrincipal					  

                   elif strain_metric == 'Shear':

                      strain_i = abs(strainValue.data[5])
                      strain_act=strainValue.data[5]					
                   if (strain_i > maxStrain):
                       
					   maxStrain = strain_i
					   retStrain = strain_act

					   maxStrainElem = strainValue.elementLabel

					   maxStep = step.name

					   maxFrame = lastFrame.incrementNumber

        odb.close()

        return retStrain	



	"""
	This ODB reading script does the following:
	-Retrieves the displacement at TIPNODE
	-Scans for max. Mises stress in a part (if set exists)
	"""
def find_stress(ModelName,elsetName):
	# Open the output database.
	print 'Made it in'
	odbName = ModelName+'.odb'
	print odbName
	odb = visualization.openOdb(odbName)
	lastFrame = odb.steps['RBM'].frames[-1]
	#elsetName='ALL_PART'
	elset = elemset = None
	region = "over the entire model"
	assembly = odb.rootAssembly

	#Check to see if the element set exists
	#in the assembly

	if elsetName:
	   try:
		   elemset = assembly.elementSets[elsetName]
		   region = " in the element set : " + elsetName;
	   except KeyError:
		   print 'An assembly level elset named %s does' \
				  'not exist in the output database %s' \
				  % (elsetName, odbName)
		   odb.close()
		   exit(0)
		   
	""" Initialize maximum values """
	maxMises = -0.1
	maxVMElem = 0
	maxStep = "_None_"
	maxFrame = -1
	Stress = 'S'
	isStressPresent = 0
	for step in odb.steps.values():
	   print 'Processing Step:', step.name
	   for frame in step.frames:
		   allFields = frame.fieldOutputs
		   if (allFields.has_key(Stress)):
			   isStressPresent = 1
			   stressSet = allFields[Stress]
			   if elemset:
				   stressSet = stressSet.getSubset(
					   region=elemset)      
			   for stressValue in stressSet.values:                
				   if (stressValue.mises > maxMises):
					   maxMises = stressValue.mises
					   maxVMElem = stressValue.elementLabel
					   maxStep = step.name
					   maxFrame = frame.incrementNumber
	if(isStressPresent):
	   print 'Maximum von Mises stress %s is %f in element %d'%(
		   region, maxMises, maxVMElem)
	   print 'Location: frame # %d  step:  %s '%(maxFrame,maxStep)
	else:
	   print 'Stress output is not available in' \
			 'the output database : %s\n' %(odb.name)  		
        odb.close()
        return maxMises	
			 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def createXYPlot(vpOrigin, vpName, plotName, data):

	##NOTE: I have never used this but it might be interesting in some problems
    
    """
    Display curves of theoretical and computed results in
    a new viewport.
    """
    print 'In plotter'
    from visualization import  USER_DEFINED
    
    vp = session.Viewport(name=vpName, origin=vpOrigin, 
        width=150, height=100)
    print 'Viewport created'
    xyPlot = session.XYPlot(plotName)
    chart = xyPlot.charts.values()[0]
    curveList = []
    for elemName, xyValues in sorted(data.items()):
        xyData = session.XYData(elemName, xyValues)
        curve = session.Curve(xyData)
        curveList.append(curve)

    print 'Curves created'
    chart.setValues(curvesToPlot=curveList)
    chart.axes1[0].axisData.setValues(useSystemTitle=False,title='Arc Height')
    chart.axes2[0].axisData.setValues(useSystemTitle=False,title=plotName)
    vp.setValues(displayedObject=xyPlot)
    print 'Plot displayed'
    return

	
   
