############################################
###creat by XJTU-Zhou-group                #
###2021/5                                  #
###creat in Abaqus Version 2020            #
############################################
# -*- coding: mbcs -*-
from abaqus import *
from abaqusConstants import *
import __main__

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from odbAccess import*
import xyPlot
import regionToolset
import displayGroupOdbToolset as dgo
import displayGroupMdbToolset as dgm
import connectorBehavior

#y0, x1 fixed（influence of the width of slot）
x0=0
y0=5

#point0 y1=y0, x1 fixed（influence of the length of slot）
x1=45
y1=y0

#point2 
x2=55.55
y2=7.52

#point3 
x3=74.30
y3=5.62

#point4 y4=0 fixed
x4=60.30
y4=0


##creat part 
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
    0.0, 0.0), point1=(19.0, 0.0))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-1'].BaseShellExtrude(depth=220, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

mdb.models['Model-1'].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=XYPLANE)
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=YZPLANE)
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=XZPLANE)
mdb.models['Model-1'].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=19.0, principalPlane=YZPLANE)

##creat slot
mdb.models['Model-1'].parts['Part-1'].DatumCsysByThreePoints(coordSysType=
    CARTESIAN, line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0), name=
    'Datum csys-1', origin=(0.0, 0.0, 0.0))

mdb.models['Model-1'].ConstrainedSketch(gridSpacing=22.19, name='__profile__', 
    sheetSize=887.9, transform=
    mdb.models['Model-1'].parts['Part-1'].MakeSketchTransform(
    sketchPlane=mdb.models['Model-1'].parts['Part-1'].datums[5], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=mdb.models['Model-1'].parts['Part-1'].datums[6].axis2, 
    sketchOrientation=RIGHT, origin=(19.0, 0.0, 110.0)))

mdb.models['Model-1'].parts['Part-1'].projectReferencesOntoSketch(filter=
    COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(0.0, 
    0.0), point2=(1.0, 0.0))
mdb.models['Model-1'].sketches['__profile__'].HorizontalConstraint(
    addUndoState=False, entity=
    mdb.models['Model-1'].sketches['__profile__'].geometry[2])
mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(0.0, 
    0.0), point2=(0.0, 1.0))
mdb.models['Model-1'].sketches['__profile__'].VerticalConstraint(addUndoState=
    False, entity=mdb.models['Model-1'].sketches['__profile__'].geometry[3])

mdb.models['Model-1'].sketches['__profile__'].Line(point1=(x0, y0), point2=(45-12, y1))#g4
mdb.models['Model-1'].sketches['__profile__'].Spline(points=((45-12, y1), (x2-12, y2), (x3-12, y3), (x4-12, y4),(x3-12,-y3),(x2-12,-y2),(45-12,-y1)))#g5
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(x0, -y0), point2=(45-12, -y1))#g6

mdb.models['Model-1'].sketches['__profile__'].copyMirror(mirrorLine=
    mdb.models['Model-1'].sketches['__profile__'].geometry[3], objectList=(
    mdb.models['Model-1'].sketches['__profile__'].geometry[4], 
    mdb.models['Model-1'].sketches['__profile__'].geometry[5],
    mdb.models['Model-1'].sketches['__profile__'].geometry[6]))
mdb.models['Model-1'].parts['Part-1'].CutExtrude(flipExtrudeDirection=OFF, 
    sketch=mdb.models['Model-1'].sketches['__profile__'], sketchOrientation=
    RIGHT, sketchPlane=mdb.models['Model-1'].parts['Part-1'].datums[5], 
    sketchPlaneSide=SIDE1, sketchUpEdge=
    mdb.models['Model-1'].parts['Part-1'].datums[6].axis2)
del mdb.models['Model-1'].sketches['__profile__']

mdb.models['Model-1'].parts['Part-1'].PartitionFaceByDatumPlane(datumPlane=
    mdb.models['Model-1'].parts['Part-1'].datums[3], faces=
    mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(('[#1 ]', 
    ), ))
mdb.models['Model-1'].parts['Part-1'].PartitionFaceByDatumPlane(datumPlane=
    mdb.models['Model-1'].parts['Part-1'].datums[4], faces=
    mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(('[#3 ]', 
    ), ))
    
##setting ABD matrix

mdb.models['Model-1'].GeneralStiffnessSection(applyThermalStress=0, density=3.18e-10, name='Section-1', poissonDefinition=DEFAULT, referenceTemperature=None, stiffnessMatrix=(7714.0, 6380.0, 7714.0, 0.0, 0.0, 5962.0, 0.0, 0.0, 0.0, 23.6, 0.0, 0.0, 0.0, 19.1, 23.6, 0.0, 0.0, 0.0, 0.0, 0.0, 19.9), useDensity=ON)
mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=Region(faces=mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(mask=('[#f ]', ), )), sectionName='Section-1', thicknessAssignment=FROM_SECTION)

##creat local coordinate system & assign the direction of the material

mdb.models['Model-1'].parts['Part-1'].DatumCsysByThreePoints(coordSysType=
    CYLINDRICAL, name='Datum csys-2', origin=(0.0, 0.0, 0.0), point1=(0.0, 0.0, 
    1.0), point2=(1.0, 0.0, 0.0))

mdb.models['Model-1'].parts['Part-1'].MaterialOrientation(
    additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
    , axis=AXIS_3, fieldName='', localCsys=
    mdb.models['Model-1'].parts['Part-1'].datums[11], orientationType=SYSTEM, 
    region=Region(
    faces=mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(
    mask=('[#f ]', ), )))  

##assembling instances

mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1', 
    part=mdb.models['Model-1'].parts['Part-1'])


##dynamic steps

mdb.models['Model-1'].ExplicitDynamicsStep(improvedDtMethod=ON, 
    linearBulkViscosity=0.06, massScaling=((SEMI_AUTOMATIC, MODEL, 
    THROUGHOUT_STEP, 0.0, 1e-06, BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ), name=
    'folding', previous='Initial', quadBulkViscosity=1.2)
mdb.models['Model-1'].ExplicitDynamicsStep(improvedDtMethod=ON, 
    linearBulkViscosity=0.06, massScaling=((SEMI_AUTOMATIC, MODEL, 
    THROUGHOUT_STEP, 0.0, 1e-06, BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ), name=
    'deployment', previous='folding', quadBulkViscosity=1.2)
mdb.models['Model-1'].steps['folding'].setValues(timePeriod=3.0)
mdb.models['Model-1'].steps['deployment'].setValues(timePeriod=3.0)


mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(filter=
    ANTIALIASING, timeInterval=0.01, variables=('CF', 'RF', 'RM', 'RT', 'SE', 
    'SF', 'U', 'UR', 'UT'))
mdb.models['Model-1'].historyOutputRequests['H-Output-1'].setValues(filter=
    ANTIALIASING, timeInterval=0.01, variables=('ALLAE', 'ALLCD', 'ALLCW', 
    'ALLDC', 'ALLDMD', 'ALLFD', 'ALLIE', 'ALLKE', 'ALLMW', 'ALLPD', 'ALLPW', 
    'ALLSE', 'ALLVD', 'ALLWK'))
    
##creat interaction  
 
mdb.models['Model-1'].ContactProperty('IntProp-1')
mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
    formulation=FRICTIONLESS)
mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
    allowSeparation=ON, constraintEnforcementMethod=DEFAULT, 
    pressureOverclosure=HARD)
mdb.models['Model-1'].ContactExp(createStepName='Initial', name='Int-1')
mdb.models['Model-1'].interactions['Int-1'].includedPairs.setValuesInStep(
    stepName='Initial', useAllstar=ON)
mdb.models['Model-1'].interactions['Int-1'].contactPropertyAssignments.appendInStep(
    assignments=((GLOBAL, SELF, 'IntProp-1'), ), stepName='Initial')
    
##creat reference points 

mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(0, 0, 220))
mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(0, 0.0, 0))
mdb.models['Model-1'].rootAssembly.ReferencePoint(point=(0, 0.0, 110))
mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='RP-1', toName='RP-A')
mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='RP-2', toName='RP-B')
mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='RP-3', toName='RP-C')
mdb.models['Model-1'].rootAssembly.Set(name='RP-A', referencePoints=(mdb.models['Model-1'].rootAssembly.referencePoints[4], ))
mdb.models['Model-1'].rootAssembly.Set(name='RP-B', referencePoints=(mdb.models['Model-1'].rootAssembly.referencePoints[5], ))
mdb.models['Model-1'].rootAssembly.Set(name='RP-C', referencePoints=(mdb.models['Model-1'].rootAssembly.referencePoints[6], ))

##creat coupling controling & equation controling of reference points

mdb.models['Model-1'].Coupling(controlPoint=Region(referencePoints=(
    mdb.models['Model-1'].rootAssembly.referencePoints[4], )), couplingType=
    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
    'Constraint-1', surface=Region(
    side1Edges=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.getSequenceFromMask(
    mask=('[#20200408 ]', ), )), u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
mdb.models['Model-1'].Coupling(controlPoint=Region(referencePoints=(
    mdb.models['Model-1'].rootAssembly.referencePoints[5], )), couplingType=
    KINEMATIC, influenceRadius=WHOLE_SURFACE, localCsys=None, name=
    'Constraint-2', surface=Region(
    side1Edges=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.getSequenceFromMask(
    mask=('[#10800802 ]', ), )), u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
mdb.models['Model-1'].Equation(name='Constraint-3', terms=((1.0, 'RP-A', 4), (
    -1.0, 'RP-B', 4), (-1.0, 'RP-C', 4)))
    
##visous pressure damping 

mdb.models['Model-1'].Pressure(amplitude=UNSET, createStepName='folding', 
    distributionType=VISCOUS, field='', magnitude=1.4623e-06, name='Load-1', 
    region=Region(
    side1Faces=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.getSequenceFromMask(
    mask=('[[#f ] ]', ), )))

    
##creat boundary condition

mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-1', 
    region=mdb.models['Model-1'].rootAssembly.sets['RP-A'], u1=SET, u2=SET, u3=
    SET, ur1=UNSET, ur2=SET, ur3=SET)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-2', 
    region=mdb.models['Model-1'].rootAssembly.sets['RP-B'], u1=SET, u2=SET, u3=
    UNSET, ur1=UNSET, ur2=SET, ur3=SET)
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-3', 
    region=mdb.models['Model-1'].rootAssembly.sets['RP-C'], u1=SET, u2=SET, u3=
    SET, ur1=UNSET, ur2=SET, ur3=SET)
    
##smoothstep amplitude curve

mdb.models['Model-1'].SmoothStepAmplitude(data=((0.0, 0.0), (3.0, 1.0), (6.0, 0.0)), name='Amp-1', timeSpan=TOTAL)
mdb.models['Model-1'].boundaryConditions['BC-3'].setValuesInStep(amplitude='Amp-1', stepName='folding', ur1=2.996)

##mesh

mdb.models['Model-1'].parts['Part-1'].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=3)
mdb.models['Model-1'].parts['Part-1'].setMeshControls(algorithm=ADVANCING_FRONT, elemShape=QUAD, regions=mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(('[[#f ] ]', ), ))
mdb.models['Model-1'].parts['Part-1'].setElementType(elemTypes=(ElemType(
    elemCode=S4R, elemLibrary=STANDARD, secondOrderAccuracy=ON, hourglassControl=STIFFNESS), ElemType(
    elemCode=S3, elemLibrary=STANDARD)), regions=(
    mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(('[[#f ] ]', 
    ), ), ))
mdb.models['Model-1'].parts['Part-1'].generateMesh()
mdb.models['Model-1'].rootAssembly.regenerate()

##creat job

mdb.Job(activateLoadBalancing=False, atTime=None, contactPrint=OFF, 
    description='', echoPrint=OFF, explicitPrecision=DOUBLE_PLUS_PACK, 
    historyPrint=OFF, memory=90, memoryUnits=PERCENTAGE, model='Model-1', 
    modelPrint=OFF, multiprocessingMode=DEFAULT, name='Job-1', 
    nodalOutputPrecision=FULL, numCpus=14, numDomains=14, 
    parallelizationMethodExplicit=DOMAIN, queue=None, resultsFormat=ODB, 
    scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
mdb.jobs['Job-1'].waitForCompletion()
mdb.saveAs(pathName='Job-1.cae')

##post-process

odb = session.openOdb(name='Job-1.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=odb)

##export transverse curvature

odb = session.openOdb(name='Job-1.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
s1f300_SK_ANTIALIASING = session.odbs['Job-1.odb'].steps['folding'].frames[300].fieldOutputs['SK_ANTIALIASING']
tmpField = abs(s1f300_SK_ANTIALIASING.getScalarField(componentLabel="SK2"))
currentOdb = session.odbs['Job-1.odb']
scratchOdb = session.ScratchOdb(odb=currentOdb)
sessionStep = scratchOdb.Step(name='Session Step', description='Step for Viewer non-persistent fields', domain=TIME, timePeriod=1.0)
sessionFrame = sessionStep.Frame(frameId=0, frameValue=0.0, description='Session Frame')
sessionField = sessionFrame.FieldOutput(name='transverse curve', description='abs(s1f300_SK_ANTIALIASING.getScalarField(componentLabel="SK2"))', field=tmpField)
odbname = session.viewports['Viewport: 1'].odbDisplay.name
frame1 = session.scratchOdbs[odbname].steps['Session Step'].frames[0]
session.viewports['Viewport: 1'].odbDisplay.setFrame(frame=frame1)
odb = session.scratchOdbs['Job-1.odb']
session.fieldReportOptions.setValues(printXYData=OFF, printTotal=OFF)
session.writeFieldReport(fileName='1_transverse curve.txt', append=ON, sortItem='单元标签', odb=odb, step=0, frame=0, outputPosition=INTEGRATION_POINT, variable=(('transverse curve', INTEGRATION_POINT), ), stepFrame=SPECIFY)
session.odbs['Job-1.odb'].close()

##export longitudinal curvature

odb = session.openOdb(name='Job-1.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
s1f300_SK_ANTIALIASING = session.odbs['Job-1.odb'].steps['folding'].frames[300].fieldOutputs['SK_ANTIALIASING']
tmpField = abs(s1f300_SK_ANTIALIASING.getScalarField(componentLabel="SK3"))
currentOdb = session.odbs['Job-1.odb']
scratchOdb = session.ScratchOdb(odb=currentOdb)
sessionStep = scratchOdb.Step(name='Session Step', description='Step for Viewer non-persistent fields', domain=TIME, timePeriod=1.0)
sessionFrame = sessionStep.Frame(frameId=0, frameValue=0.0, description='Session Frame')
sessionField = sessionFrame.FieldOutput(name='longitudinal curve', description='abs(s1f300_SK_ANTIALIASING.getScalarField(componentLabel="SK3"))', field=tmpField)
odbname = session.viewports['Viewport: 1'].odbDisplay.name
frame1 = session.scratchOdbs[odbname].steps['Session Step'].frames[0]
session.viewports['Viewport: 1'].odbDisplay.setFrame(frame=frame1)
odb = session.scratchOdbs['Job-1.odb']
session.fieldReportOptions.setValues(printXYData=OFF, printTotal=OFF)
session.writeFieldReport(fileName='2_longitudinal curve.txt', append=ON, sortItem='单元标签', odb=odb, step=0, frame=0, outputPosition=INTEGRATION_POINT, variable=(('longitudinal curve', INTEGRATION_POINT), ), stepFrame=SPECIFY)
session.odbs['Job-1.odb'].close()


##export the maximum deployment moment

odb = session.openOdb(name='Job-1.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
session.mdbData.summary()
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
odbName=session.viewports[session.currentViewportName].odbDisplay.name
session.odbData[odbName].setValues(activeFrames=(('deployment', ('0:-1', )), ))
odb = session.odbs['Job-1.odb']
xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=((
'RM_ANTIALIASING', NODAL, ((COMPONENT, 'RM1'), )), ), nodeSets=("RP-C", ))
xyp = session.XYPlot('XYPlot-1')
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
curveList = session.curveSet(xyData=xyList)
chart.setValues(curvesToPlot=curveList)
session.charts[chartName].autoColor(lines=True, symbols=True)
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
x0 = chart.curves['_RM_ANTIALIASING:RM1 PI: ASSEMBLY N: 3']
session.xyReportOptions.setValues(xyData=OFF, minMax=ON)
session.writeXYReport(fileName='3_deployment_moment.txt', xyData=(x0, ))
del session.xyDataObjects['_temp_1']
session.odbs['Job-1.odb'].close()


##export the maximum stored strain energy


session.mdbData.summary()
xyp = session.xyPlots['XYPlot-1']
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
o1 = session.openOdb(name='Job-1.odb', readOnly=False)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
odb = session.odbs['Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, outputVariableName='Strain energy: ALLSE for Whole Model', steps=('folding', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
x0 = chart.curves['_temp_1']
session.xyReportOptions.setValues(xyData=OFF, minMax=ON)
session.writeXYReport(fileName='4_strain energy.txt', xyData=(x0, ))
del session.xyDataObjects['_temp_1']
session.odbs['Job-1.odb'].close()


##The maximum kinetic energy/The maximum strain energy


odb = session.openOdb(name='Job-1.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
odb = session.odbs['Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
outputVariableName='Kinetic energy: ALLKE for Whole Model', steps=(
'folding', 'deployment', ), suppressQuery=True, 
__linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
session.xyDataObjects.changeKey(fromName='_temp_1', toName='K')
del session.xyDataObjects['_temp_1']
odb = session.odbs['Job-1.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
outputVariableName='Internal energy: ALLIE for Whole Model', steps=(
'folding', 'deployment', ), suppressQuery=True, 
__linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
session.xyDataObjects.changeKey(fromName='_temp_1', toName='I')
del session.xyDataObjects['_temp_1']
xy1 = session.xyDataObjects['K']
xy2 = session.xyDataObjects['I']
xy3 = avg((xy1))/avg((xy2))
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
c1 = session.Curve(xyData=xy3)
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
x0 = chart.curves['_temp_3']
session.xyReportOptions.setValues(xyData=OFF, minMax=ON)
session.writeXYReport(fileName='5_K_I.txt', xyData=(x0, ))
session.odbs['Job-1.odb'].close()
del session.xyDataObjects['_temp_1']
del session.xyDataObjects['_temp_2']
del session.xyDataObjects['_temp_3']


##the maximum bending stiffness


odb = session.openOdb(name='Job-1.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
session.linkedViewportCommands.setValues(_highlightLinkedViewports=False)
odbName=session.viewports[session.currentViewportName].odbDisplay.name
session.odbData[odbName].setValues(activeFrames=(('deployment', ('0:-1', )), ))
odb = session.odbs['Job-1.odb']
xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RM_ANTIALIASING', NODAL, ((COMPONENT, 'RM1'), )), ), nodeSets=("RP-C", ))
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
curveList = session.curveSet(xyData=xyList)
chart.setValues(curvesToPlot=curveList)
session.charts[chartName].autoColor(lines=True, symbols=True)
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
xy1 = session.xyDataObjects['_RM_ANTIALIASING:RM1 PI: ASSEMBLY N: 3']
xy2 = differentiate(xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
c1 = session.Curve(xyData=xy2)
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
x0 = chart.curves['_temp_1']
session.xyReportOptions.setValues(xyData=OFF, minMax=ON)
session.writeXYReport(fileName='6_bending stiffness.txt', xyData=(x0, ))
session.odbs['Job-1.odb'].close()
del session.xyDataObjects['_temp_1']
del session.xyDataObjects['_temp_2']
del session.xyDataObjects['_temp_3']

##mass of the CTSH


openMdb(pathName='Job-1.cae')
a=mdb.models['Model-1'].rootAssembly
b=a.getArea(a.instances['Part-1-1'].faces)
c=3.18E-10*b*1000
file_handle=open('8_mass.txt',mode='a+')
print >>file_handle, c
file_handle.close()