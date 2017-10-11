
# coding: utf-8

# # How To: Load the database, selecting only a few results.

# In[1]:

#Set up the path to SModelS installation folder if running on a different folder
import sys,os
sys.path.append(os.path.join(os.getenv("HOME"),"smodels/"))


# In[2]:

from smodels.experiment.databaseObj import Database
from smodels.tools.physicsUnits import GeV


# In[3]:

## Load the database:
dbPath = os.path.join(os.getenv("HOME"),"smodels-database/")
database = Database(dbPath)


# ## How to select results from one publication (or conference note)

# In[6]:

#Select only the CMS SUS-12-028 conference note
expID=["CMS-SUS-12-028"]


# In[7]:

#Loads the selected analyses
#(The INFO tells you that superseded analyses are not loaded, see below)
results = database.getExpResults(analysisIDs=expID)


# In[9]:

#Print all the results selected:
for exp in results:
    print exp
#Print the txnames constrained by the result in bracket notation:
exp = results[0]
for tx in exp.getTxNames():
    print tx,'=',tx.constraint


# In[10]:

#Print ALL info fields available:
exp = results[0]
print exp.getAttributes()


# In[11]:

#Print values for some of the info fields (always returned as a list):
print 'sqrts=',exp.getValuesFor('sqrts')
print 'lumi=',exp.getValuesFor('lumi')
print 'dataType=',exp.getValuesFor('dataType')
print 'txnames=',exp.getValuesFor('txName')


# In[12]:

#To obtain the upper limit for a given mass vector and a given simplified model (txname)
#Note that the number of masses in the mass vector must be consitent with the txname.
#For the T1 txname, for instance:
massesT1 = [[300*GeV,100*GeV],[300*GeV,100*GeV]]
print 'xsection upper limit = ',exp.getUpperLimitFor(mass=massesT1,txname='T1')


# In[13]:

#For the T2 analysis:
massesT2 = [[300*GeV,50*GeV],[300*GeV,50*GeV]]
print 'xsection upper limit = ',exp.getUpperLimitFor(mass=massesT2,txname='T2')


# In[14]:

#If you try with the wrong mass format, an error will be printed:
masses = [[300*GeV],[300*GeV,50*GeV]]
print 'xsection upper limit = ',exp.getUpperLimitFor(mass=masses,txname='T2')


# ## How to load results for one simplified model (txname)

# In[15]:

#It is also possible to load all the results for a single simplified (using the Txname convention)
Txnames = ["T1"]
T1results = database.getExpResults(txnames=Txnames)


# In[16]:

#Print all the results constraining the required Txname:
for exp in T1results:
    print exp.globalInfo.id #(or print exp.getValuesFor('id'))


# ## How to load all experimental results, including the superseded publications

# In[17]:

#By default only non-supersed analyses are loaded:
results = database.getExpResults()
print 'Number of non-superseded results = ',len(results)


# In[18]:

#To load all results (including the superseded ones), set useSuperseded=True
allResults = database.getExpResults(useSuperseded=True)
print 'Including superseded results =',len(allResults)


# ## How to selected upper-limit and efficiency map results:

# In[19]:

#Get only upper-limit results:
ULresults =  database.getExpResults(dataTypes=['upperLimit'])
for exp in ULresults:
    print exp.globalInfo.id,exp.getValuesFor('dataType')


# In[20]:

#Get only efficiency map results:
EMresults =  database.getExpResults(dataTypes=['efficiencyMap'])
for exp in EMresults:
    print exp.globalInfo.id,exp.getValuesFor('dataType')


# In[ ]:



