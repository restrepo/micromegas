
# coding: utf-8

# # How To: Look up the efficinecy of a particular result, for a particular set of masses

# In[1]:

#Set up the path to SModelS installation folder if running on a different folder
import sys,os
sys.path.append(os.path.join(os.getenv("HOME"),"smodels/"))


# In[2]:

#Import those parts of smodels that are needed for this exercise
from smodels.tools.physicsUnits import GeV
from smodels.experiment.databaseObj import Database


# In[3]:

## Load the database:
databasePath = os.path.join(os.getenv("HOME"),"smodels-database/")
db = Database(databasePath)


# ## Look up efficiency for an Upper Limit-type result:

# In[4]:

#Select desired result:
resultID = ["CMS-PAS-SUS-13-016"]
txname = ["T1tttt"]
expResult = db.getExpResults(analysisIDs=resultID,txnames=txname,dataTypes='upperLimit')[0]
print 'selected result:',expResult


# In[5]:

#Define the desired mass vector (must be consistent with the txname/simplified model):
massesA = [[500*GeV, 150*GeV],[500*GeV, 150*GeV]]
massesB = [[5000*GeV, 150*GeV],[5000*GeV, 150*GeV]]
#For UL-type results, the efficiency is 1, if the mass is inside the grid or zero if it is outside:
print 'efficiency for mass\n',massesA,' is: ',expResult.getEfficiencyFor(mass=massesA,txname="T1tttt")
print 'efficiency for mass\n',massesB,' is: ',expResult.getEfficiencyFor(mass=massesB,txname="T1tttt")


# ## Look up efficiency for an Efficiency Map-type result:

# In[6]:

#Select desired result:
resultID = ["CMS-PAS-SUS-13-016"]
txname = ["T1tttt"]
expResult = db.getExpResults(analysisIDs=resultID,txnames=txname,dataTypes='efficiencyMap')[0]
print 'selected result:',expResult


# In[12]:

#For an efficiency map result one needs to specify the desired signal region (dataset) and mass
masses = [[500*GeV, 150*GeV],[500*GeV, 150*GeV]]
datasetID = 'sr0'
print 'efficinecy for mass\n ',masses,'\n in dataset',datasetID,' is: ',expResult.getEfficiencyFor(mass=masses,txname="T1tttt",dataset=datasetID)


# In[ ]:



