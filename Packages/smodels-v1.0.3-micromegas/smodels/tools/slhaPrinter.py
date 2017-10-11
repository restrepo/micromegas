#!/usr/bin/env python

"""
    .. module:: slhaPrinter
    :synopsis: Prints result, missing topologies in slha format
    .. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>
"""

from smodels.tools.physicsUnits import TeV, fb, GeV
from smodels.installation import version

def formatSLHAInput(maxcond, minmassgap, sigmacut, databaseversion, docompress):
    try:
        smodelsversion = version()
        if not smodelsversion.startswith("v"): smodelsversion = "v" + smodelsversion
    except:
	    smodelsversion = "unknown"

    output = "BLOCK SModelS_Settings\n"
    output += " 0 %-20s #SModelS version\n" %(smodelsversion)
    output += " 1 %-20s #database version\n" %(databaseversion)
    output += " 2 %-20s #maximum condition violation\n" % (maxcond)
    output += " 3 %-20s #compression (0 off, 1 on)\n" % (docompress)
    output += " 3 %-20s #minimum mass gap for mass compression [GeV]\n" % (minmassgap/GeV)
    output += " 4 %-20s #sigmacut [fb]\n\n" % (sigmacut/fb)
    return output

def formatSLHAResults(resultList):
    output = "BLOCK SModelS_Exclusion\n"
    if not resultList or resultList.isEmpty():
        excluded = -1
    else:
        firstResult = resultList.getBestResult()
        r = resultList.getR(firstResult)
        if r > 1: excluded = 1
        else: excluded = 0
    output += " 0 0 %-20s #output status (-1 not tested, 0 not excluded, 1 excluded)\n\n" % (excluded)
    if excluded == 0: rList = [firstResult]
    elif excluded == 1: rList = resultList.outputarray
    else: rList = []
    cter = 1
    for res in rList:
       if resultList.getR(res) <1 and not excluded == 0: break
       output += " %d 0 %-20s #txname (upper limit)\n" % (cter, res.analysis.label.split(":")[1] )
       output += " %d 1 %-20.3E #r value\n" % (cter, resultList.getR(res))
       output += " %d 2 %-20.2f #condition violation\n" % (cter, res.getmaxCondition())
       output += " %d 3 %-20s #analysis\n" % (cter, res.analysis.label.split(":")[0])
       output += "\n"
       cter += 1
    return output

def formatSLHAMissing(missingtopos):
    if not missingtopos: return ""
    if not missingtopos.topos: return ""
    output = "BLOCK SModelS_Missing_Topos #sqrts[TeV] weight[fb] description\n"
    cter = 0
    for t in sorted(missingtopos.topos, key=lambda x: x.value, reverse=True):
        output += " %d %d %10.3E %s\n" % (cter, missingtopos.sqrts/TeV, t.value, str(t.topo))
        cter += 1
        if cter > 9: break
    return output

def writeSLHA(resultList, maxcond, minmassgap, sigmacut, missingtopos, databaseversion, docompress, outfile="summary.slha"):
    f = open(outfile, "w")
    f.write(formatSLHAInput(maxcond, minmassgap, sigmacut, databaseversion, docompress))
    f.write(formatSLHAResults(resultList))
    f.write(formatSLHAMissing(missingtopos))
    f.close()
    return None

