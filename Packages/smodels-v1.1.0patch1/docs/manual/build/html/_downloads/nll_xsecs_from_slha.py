
# coding: utf-8

# # How To: Compute NLL cross sections for a given SLHA file

# In[1]:

#Set up the path to SModelS installation folder if running on a different folder
import sys,os
sys.path.append(os.path.join(os.getenv("HOME"),"smodels/"))


# In[2]:

#Import those parts of smodels that are needed for this exercise
from smodels.tools import xsecComputer
from smodels.tools.physicsUnits import TeV, fb
from smodels.installation import installDirectory
from smodels.tools.xsecComputer import LO, NLL


# In[3]:

#Define the SLHA file name
filename="%s/inputFiles/slha/gluino_squarks.slha" % installDirectory()


# In[6]:

#Lets compute the NLL cross-sections for 8 TeV. The xsecComputer will first use Pythia to compute
#the LO cross-sections and then NLLfast to compute the k-factors.
#The output will contain only the processes contained in NLLfast (gluino and squark production)
#For the Pythia step we have to define the number of MC events (1k)
#(To see how to also include LO cross-sections for all processes check the LO HowTo)
NLL = 2
xsecsNLL=xsecComputer.computeXSec(sqrts = 8*TeV, maxOrder=NLL, nevts=1000, slhafile=filename )


# In[7]:

# the output is a XSectionList ...
type(xsecsNLL)


# In[8]:

#Each entry in the list contains the cross-section value:
print(xsecsNLL[0].value)
#The PDGs of the particles produced:
print(xsecsNLL[0].pid)
#And some additional info
print("label =",xsecsNLL[0].info.label,"Sqrts =",xsecsNLL[0].info.sqrts, "QCD order =",xsecsNLL[0].info.order)


# In[9]:

#Finally, lets write the cross-sections back to the file 
#(will write only the cross-sections not overlapping the existing ones):
xsecComputer.addXSecToFile(xsecsNLL,filename)


# In[ ]:



