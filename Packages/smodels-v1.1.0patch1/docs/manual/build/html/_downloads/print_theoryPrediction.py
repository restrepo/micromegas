
# coding: utf-8

# # How To: Print out the theory predictions

# In[1]:

#Set up the path to SModelS installation folder if running on a different folder
import sys,os
sys.path.append(os.path.join(os.getenv("HOME"),"smodels/"))


# In[2]:

#Import those parts of smodels that are needed for this exercise
#(We will assume the input is a SLHA file. For LHE files, use the lheDecomposer instead)
from smodels.theory import slhaDecomposer
from smodels.installation import installDirectory
from smodels.tools.physicsUnits import fb, GeV
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObj import Database


# In[3]:

## Define the database path
databasepath = os.path.join(os.getenv("HOME"),"smodels-database/")
#and load the database:
database = Database(databasepath)
#Get list of desired experimental results (by default all results):
expResults = database.getExpResults()


# In[4]:

#Define the SLHA input file name
filename="%s/inputFiles/slha/gluino_squarks.slha" % installDirectory()


# In[5]:

#Perform the decomposition:
topList = slhaDecomposer.decompose(filename, sigcut = 0.03 * fb, doCompress=True, doInvisible=True,minmassgap = 5* GeV)


# In[6]:

#Compute the theory prediction for each selected experimental result using the output from the decomposition:
allThPredictions = [theoryPredictionsFor(exp, topList) for exp in expResults]


# In[7]:

#Print information about each theory prediction for each result:
for thPreds in allThPredictions:
    if not thPreds: continue #skip results with no predictions
    for theoryPred in thPreds:
        print "\nExperimental Result: ",theoryPred.expResult.globalInfo.id,"(%s-type)" %theoryPred.dataset.dataInfo.dataType #Result ID
        print "Theory prediction xsec = ",theoryPred.xsection.value #Signal xsection*efficiency*BR
        print "Conditions violation (if any) = ",theoryPred.conditions #Condition values (for UL-type results)
        print "Dataset:",theoryPred.dataset.dataInfo.dataId  #Corresponding signal region (for EM-type results)
        print "Txnames = ",[str(tx) for tx in theoryPred.txnames] #List of simplified models (txnames) contributing to the signal xsec        


# In[ ]:



