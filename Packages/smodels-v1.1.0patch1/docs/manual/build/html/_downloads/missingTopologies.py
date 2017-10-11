
# coding: utf-8

# # How To: Find missing topologies that are not covered by the database

# In[1]:

# Set up the path to SModelS installation folder if running on a different folder
import sys,os
sys.path.append(os.path.join(os.getenv("HOME"),"smodels/"))


# In[2]:

# Import those parts of smodels that are needed for this exercise
from smodels.tools.physicsUnits import TeV, GeV, fb
from smodels.installation import installDirectory
from smodels.theory import slhaDecomposer
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.experiment.databaseObj import Database
from smodels.tools import coverage


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
#(Computing theory predictions will tag the elements which have been tested)
allThPredictions = [theoryPredictionsFor(exp, topList) for exp in expResults]


# In[5]:

# Create missing Topologies object
uncovered = coverage.Uncovered(topList)


# In[6]:

#Print basic information about coverage:
print "Total missing topology cross section (fb): %10.3E" %(uncovered.getMissingXsec())
print "Total cross section where we are outside the mass grid (fb): %10.3E" %(uncovered.getOutOfGridXsec())
print "Total cross section in long cascade decays (fb): %10.3E" %(uncovered.getLongCascadeXsec())
print "Total cross section in decays with asymmetric branches (fb): %10.3E" %(uncovered.getAsymmetricXsec())        


# In[7]:

# Get list of topologies which are not tested by any result:
missingTopos = uncovered.missingTopos
#Print information about the first few missing topologies and
#the elements contributing to the topology: 
for top in missingTopos.topos[:5]:
    print '\nmissing topology:',top.topo
    print 'Contributing elements:'
    for el in sorted(top.contributingElements, key=lambda el: el.weight):
        print el,', xsection:',el.weight


# In[10]:

# Get list of topologies which are outside the upper limit and efficiency map grids:
outsideGrid = uncovered.outsideGrid
#Print information about the first few missing topologies and
#the elements contributing to the topology: 
for top in outsideGrid.topos[:5]:
    print '\noutside the grid topology:',top.topo
    print 'Contributing elements:'
    for el in sorted(top.contributingElements, key=lambda el: el.weight):
        print el,'mass=',el.getMasses()

