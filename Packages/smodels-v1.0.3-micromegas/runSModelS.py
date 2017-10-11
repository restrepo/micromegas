#!/usr/bin/env python

"""
.. module:: runSModelS
   :synopsis: Main code for running SModelS.
   
"""

import sys, os, logging
from smodels.installation import installDirectory
sys.path.append(installDirectory()+"/Unum-4.1.3-py2.7.egg/")
import argparse
from ConfigParser import SafeConfigParser
from smodels.tools.physicsUnits import GeV, fb
from smodels.tools import slhaPrinter
log = logging.getLogger(__name__)


def main(inputFile, parameterFile, outputFile, slhaOutputFile, particlePath):
    """
    Provides a command line interface to basic SModelS functionalities.
    
    :param inputFile: input file name (either a SLHA or LHE file)
    :param parameterFile: File containing the input parameters (default = /etc/parameters_default.ini)
    :param outputFile: Output file to write a summary of results
    :param slhaOutputFile: Output file to write SLHA type summary of results
    :param particlePath: Path to directory where particles.py is stored
    
    """

    """
    Read and check input file
    =========================
    """
    parser = SafeConfigParser()
    parser.read(parameterFile)

    """ Minimum value of cross-section for an element to be considered eligible for decomposition.
        Too small sigmacut leads to too large decomposition time. """
    sigmacut = parser.getfloat("parameters", "sigmacut") * fb

    """ Minimum value for considering two states non-degenerate (only used for mass compression) """
    minmassgap = parser.getfloat("parameters", "minmassgap") * GeV

    if os.path.exists(outputFile):
        log.warning("Removing old output file in " + outputFile)
    outfile = open(outputFile, 'w')
    outfile.close()

    databaseVersion = "unknown" # set default database version that is printed in case of errors

    """ Set doCompress flag, only used for slha type output """
    if parser.getboolean("options", "doCompress") or parser.getboolean("options", "doInvisible"): docompress = 1
    else: docompress = 0

    """
    check if particles.py exists in specified path, and add to sys.path
    """
    if not os.path.isfile(os.path.join(particlePath,"particles.py")):
        log.error("particle.py not found in %s" %particlePath )
        return slhaPrinter.writeSLHA(None, parser.getfloat("parameters", "maxcond"), minmassgap, sigmacut, None, databaseVersion, docompress, slhaOutputFile)
    else:
        sys.path.insert(1, particlePath)
        from smodels.tools import ioObjects, missingTopologies
        from smodels.experiment import smsHelpers, smsAnalysisFactory
        from smodels.theory import slhaDecomposer, lheDecomposer
        from smodels.theory.theoryPrediction import theoryPredictionFor


    inputType = parser.get("options", "inputType").lower()
    if inputType != 'slha' and inputType != 'lhe':
        log.error("Unknown input type (must be SLHA or LHE): %s" % inputType)
        return slhaPrinter.writeSLHA(None, parser.getfloat("parameters", "maxcond"), minmassgap, sigmacut, None, databaseVersion, docompress, slhaOutputFile)

    """ Check input file for errors """
    inputStatus = ioObjects.FileStatus()
    if parser.getboolean("options", "checkInput"):
        inputStatus.checkFile(inputType, inputFile, sigmacut)

    """ Check database location """
    try:
        smsHelpers.base = parser.get("path", "databasePath")
        if smsHelpers.base == "./smodels-database" or smsHelpers.base == "./smodels-database/": smsHelpers.base = installDirectory()+"/smodels-database/"
        databaseVersion = smsHelpers.databaseVersion()
    except:
        log.error("Database not found in %s" % os.path.realpath(smsHelpers.base))
        return slhaPrinter.writeSLHA(None, parser.getfloat("parameters", "maxcond"), minmassgap, sigmacut, None, databaseVersion, docompress, slhaOutputFile)

    """ Initialize output status and exit if there were errors in the input """
    outputStatus = ioObjects.OutputStatus(inputStatus.status, inputFile, dict(parser.items("parameters")), databaseVersion, outputFile)
    if outputStatus.status < 0:
        return slhaPrinter.writeSLHA(None, parser.getfloat("parameters", "maxcond"), minmassgap, sigmacut, None, databaseVersion, docompress, slhaOutputFile)

    """
    Decompose input file
    ====================
    """
    try:
        """ Decompose input SLHA file, store the output elements in smstoplist """
        if inputType == 'slha':
            smstoplist = slhaDecomposer.decompose(inputFile, sigmacut, doCompress=parser.getboolean("options", "doCompress"),
                         doInvisible=parser.getboolean("options", "doInvisible"), minmassgap=minmassgap)
        else:
            smstoplist = lheDecomposer.decompose(inputFile, doCompress=parser.getboolean("options", "doCompress"),
                         doInvisible=parser.getboolean("options", "doInvisible"), minmassgap=minmassgap)
    except:
        """ Update status to fail, print error message and exit """
        outputStatus.updateStatus(-1)
        return slhaPrinter.writeSLHA(None, parser.getfloat("parameters", "maxcond"), minmassgap, sigmacut, None, databaseVersion, docompress, slhaOutputFile)

    """ Print Decomposition output.
        If no topologies with sigma > sigmacut are found, update status, write output file, stop running """
    if not smstoplist:
        outputStatus.updateStatus(-3)
        return slhaPrinter.writeSLHA(None, parser.getfloat("parameters", "maxcond"), minmassgap, sigmacut, None, databaseVersion, docompress, slhaOutputFile)

    outLevel = 0
    if parser.getboolean("stdout", "printDecomp"):
        outLevel = 1
        outLevel += parser.getboolean("stdout", "addElmentInfo")
    smstoplist.printout(outputLevel=outLevel)


    """
    Load analysis database
    ======================
    """
    
    """ In case that a list of analyses or txnames are given, retrieve list """
    analyses = parser.get("database", "analyses")
    if "," in analyses:
        analyses = analyses.split(",")
    txnames = parser.get("database", "txnames")
    if "," in txnames:
        txnames = txnames.split(",")
    
    """ Load analyses """
    listofanalyses = smsAnalysisFactory.load(analyses, txnames)

    """ Print list of analyses loaded """
    if parser.getboolean("stdout", "printAnalyses"):
        outLevel = 1
        outLevel += parser.getboolean("stdout", "addAnaInfo")
        print("=======================\n == List of Analyses   ====\n ================")
        for analysis in listofanalyses:
            analysis.printout(outputLevel=outLevel)


    """
    Compute theory predictions and anlalyses constraints
    ====================================================
    """

    """ Define result list that collects all theoryPrediction objects.
        Variables set to define printing options. """
    results = ioObjects.ResultList(bestresultonly=not parser.getboolean("file", "expandedSummary"),
                                   describeTopo=parser.getboolean("file", "addConstraintInfo"))

    """ Get theory prediction for each analysis and print basic output """
    for analysis in listofanalyses:
        theorypredictions = theoryPredictionFor(analysis, smstoplist)
        if not theorypredictions:
            continue
        if parser.getboolean("stdout", "printResults"):
            print "================================================================================"
            theorypredictions.printout()
        print "................................................................................"

        """ Create a list of results, to determine the best result """
        for theoryprediction in theorypredictions:
            results.addResult(theoryprediction, maxcond=parser.getfloat("parameters", "maxcond"))

    """ If there is no best result, this means that there are no matching experimental results for the point """
    if results.isEmpty():
        """ no experimental constraints found """
        outputStatus.updateStatus(0)
    else:
        outputStatus.updateStatus(1)

    """ Write output file """
    outputStatus.printout("file", outputFile)
    """ Add experimental constraints if found """
    if outputStatus.status == 1:
        results.printout("file", outputFile)

    sqrts = max([xsec.info.sqrts for xsec in smstoplist.getTotalWeight()])
    if parser.getboolean("options", "findMissingTopos"):
        """ Look for missing topologies, add them to the output file """
        missingtopos = missingTopologies.MissingTopoList(sqrts)
        missingtopos.findMissingTopos(smstoplist, listofanalyses, minmassgap, parser.getboolean("options", "doCompress"),
                         doInvisible=parser.getboolean("options", "doInvisible"))
        missingtopos.printout("file", outputFile)
    slhaPrinter.writeSLHA(results, parser.getfloat("parameters", "maxcond"), minmassgap, sigmacut, missingtopos, databaseVersion, docompress, slhaOutputFile)


if __name__ == "__main__":
    """ Set default input and output files """
    parameterFile = "%s/etc/parameters_default.ini" % installDirectory()
    outputFile = "summary.txt"
    slhaOutputFile = "summary.slha"

    """ Get the name of input slha file and parameter file """
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-f', '--filename', help='name of SLHA or LHE input file, necessary input', required=True)
    argparser.add_argument('-p', '--parameterFile',
                            help='name of parameter file, optional argument, if not set, use all parameters from etc/parameters_default.ini',
                            default=parameterFile)
    argparser.add_argument('-o', '--outputFile', help='name of output file, optional argument, default is: ' + outputFile,
                           default=outputFile)
    argparser.add_argument('-s', '--slhaOutputFile', help='name of SLHA type output file', default=slhaOutputFile)
    argparser.add_argument('-particles', '--particlePath', help='path to directory where particles.py is located', default=installDirectory()+"/smodels")
    args = argparser.parse_args()

    main(args.filename, args.parameterFile, args.outputFile, args.slhaOutputFile, args.particlePath)
