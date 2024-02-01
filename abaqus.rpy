# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2018 replay file
# Internal Version: 2017_11_07-17.21.41 127140
# Run by vk20927 on Thu Aug 31 00:36:46 2023
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=183.06640625, 
    height=70.4861145019531)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
openMdb('Elastic.cae')
#: The model database "C:\Temp\forXp\Verification\CDM\nonlinearShear\Elastic.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Damage'].parts['indentor']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
mdb.models.changeKey(fromName='Damage', toName='Elastic')
p = mdb.models['Elastic'].parts['indentor']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
mdb.save()
#: The model database has been saved to "C:\Temp\forXp\Verification\CDM\nonlinearShear\Elastic.cae".
