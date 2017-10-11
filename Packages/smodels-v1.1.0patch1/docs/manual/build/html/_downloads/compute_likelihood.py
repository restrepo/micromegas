
# coding: utf-8

# # How To: Compute likelihood and chi2

# In[1]:

#Set up the path to SModelS installation folder if running on a different folder
import sys,os
sys.path.append(os.path.join(os.getenv("HOME"),"smodels/"))


# In[2]:

from smodels.tools.statistics import likelihood,chi2


# ## How to compute the likelihood and chi2 from the number of events

# In[5]:

#If the number of observed events, the number of expected background events,
#its error and the number of signal events and its error are known, the likelihood
#for the signal (assuming a truncated gaussian distribution for the background and signal uncertainties)
#can be computed as:
nSignal = 10.3
deltaSignal = 1.
nObserved = 5
nBG = 4.2
deltaBG = 0.71
print 'likelihood=',likelihood(nSignal, nObserved, nBG, deltaBG, deltaSignal)
print 'chi2=',chi2(nSignal, nObserved, nBG, deltaBG, deltaSignal)


# ## How to compute the likelihood and chi2 from a theory prediction

# In[6]:

#In most cases one wants to compute the likelihood and chi2 for a given theory prediction computed by SModelS.
#Below we generate theory predictions and compute the likelihood and chi2 values for them
#First we import those parts of smodels that are needed for this exercise
#(We will assume the input is a SLHA file. For LHE files, use the lheDecomposer instead)
from smodels.theory import slhaDecomposer
from smodels.installation import installDirectory
from smodels.tools.physicsUnits import fb, GeV
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObj import Database


# In[7]:

#Define the SLHA input file name
filename="%s/inputFiles/slha/gluino_squarks.slha" % installDirectory()


# In[8]:

#Load the database, do the decomposition and compute theory predictions:
#(Look at the theory predictions HowTo to learn how to compute theory predictions)
databasepath = os.path.join(os.getenv("HOME"),"smodels-database/")
database = Database(databasepath)
expResults = database.getExpResults()
topList = slhaDecomposer.decompose(filename, sigcut = 0.03 * fb, doCompress=True, doInvisible=True,minmassgap = 5* GeV)
allThPredictions = [theoryPredictionsFor(exp, topList) for exp in expResults]


# In[13]:

#For each theory prediction, compute the corresponding likelihood and chi2 values
#(This is only possible for efficiency map-type results):
for i,thPreds in enumerate(allThPredictions):
    if not thPreds: continue #skip results with no predictions
    expID = expResults[i].globalInfo.id
    dataType = expResults[i].getValuesFor('dataType')[0]    
    for theoryPred in thPreds:
        #Compute the likelihood and chi2:
        theoryPred.computeStatistics()
        print "\nExperimental Result: %s (%s-type)" %(expID,dataType) #Result ID and type
        print "Theory prediction xsec = ",theoryPred.xsection.value #Signal xsection*efficiency*BR
        if dataType == 'efficiencyMap':
            print 'likelihood =',theoryPred.likelihood,', chi2 =',theoryPred.chi2
        else:
            print "(likelihood not available)"


# In[ ]:



