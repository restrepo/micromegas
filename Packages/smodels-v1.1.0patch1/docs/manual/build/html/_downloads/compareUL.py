
# coding: utf-8

# # How To: Compare theory predictions with experimental limits

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

#Define the SLHA input file name
filename="%s/inputFiles/slha/gluino_squarks.slha" % installDirectory()


# In[4]:

#Load the database, do the decomposition and compute theory predictions:
#(Look at the theory predictions HowTo to learn how to compute theory predictions)
databasepath = os.path.join(os.getenv("HOME"),"smodels-database/")
database = Database(databasepath)
expResults = database.getExpResults()
topList = slhaDecomposer.decompose(filename, sigcut = 0.03 * fb, doCompress=True, doInvisible=True,minmassgap = 5* GeV)
allThPredictions = [theoryPredictionsFor(exp, topList) for exp in expResults]


# In[5]:

#Print the value of each theory prediction for each experimental
#result and the corresponding upper limit (see the obtain experimental upper limits HowTo to learn how
#to compute the upper limits).
#Also print the expected upper limit, if available
for thPreds in allThPredictions:
    if not thPreds: continue #skip results with no predictions
    for theoryPred in thPreds:
        expID = theoryPred.expResult.globalInfo.id
        dataType = theoryPred.dataset.dataInfo.dataType
        print "\nExperimental Result: %s (%s-type)" %(expID,dataType) #Result ID and type
        print "Theory prediction xsec = ",theoryPred.xsection.value #Signal xsection*efficiency*BR
        print "Upper limit = ",theoryPred.expResult.getUpperLimitFor(mass = theoryPred.mass, 
                                                              txname=theoryPred.txnames[0], 
                                                              dataID = theoryPred.dataset.dataInfo.dataId)
        print "Expected Upper limit = ",theoryPred.expResult.getUpperLimitFor(mass = theoryPred.mass, 
                                                              txname=theoryPred.txnames[0], 
                                                              dataID = theoryPred.dataset.dataInfo.dataId,
                                                               expected= True)        


# In[ ]:



