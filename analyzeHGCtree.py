#!/usr/bin/env python2

import glob, sys, os, array
# add -b to sys argument to force run in batch mode
#sys.argv.append( '-b' )
import ROOT as rt

def getTreeFromFiles(opts):

    chain = rt.TChain(opts.treename)
    for fname in [opts.infile]:
        if opts.verbose > 0: print("#Adding file %s" %fname)

        chain.Add(fname)

    if opts.verbose > 0:
        print("#Added %i files to chain %s" %(chain.GetNtrees(),opts.treename))
        print("#Chain contains %i entries" % chain.GetEntries())

    return chain

def anaTree(tree, opts):
    "Analyze entry"

    if opts.maxEntries == -1: opts.maxEntries = tree.GetEntries()
    if opts.verbose > 0: print("#Going to analyze %i entries" %opts.maxEntries)

    for ientry, entry in enumerate(tree):
        if opts.verbose > 1: print ("#Analyzing entry number %i" %ientry)
        if ientry > opts.maxEntries-1: break

        ## Do analysis

    return 1

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()

    parser.usage = '%prog [options]'
    parser.description="""
    Analyse HGCAL flat tuple
    """

    # General options
    parser.add_option("-b","--batch", dest="batch",default=False, action="store_true", help="Batch mode")
    parser.add_option("-f","--file", dest="infile",default="root://eosuser.cern.ch//eos/user/b/beaudett/local/hgcal/singleElePt15-extendedNtp.root",help="Input File name(s)")
    parser.add_option("-t","--tname", dest="treename",default="ana/hgc",help="Tree name")
    parser.add_option("-n","--maxEve","--maxEvents", "--maxEntries",dest="maxEntries", default=-1,  type="int",    help="Maximum entries to analyze")
    parser.add_option("-v","--verbose",  dest="verbose",  default=1,  type="int",    help="Verbosity level (0 = quiet, 1 = verbose, 2+ = more)")

    ## Options examples (NOT USED)
    # bool
    #parser.add_option("--mc","--mcData", dest="mcData",default=False, action="store_true", help="Use pseudo-data from MC")
    # string
    #parser.add_option("--cuts", dest="cutFile",default="trig_base.txt",help="Baseline cuts file")
    # int/floats
    #parser.add_option("-v","--verbose",  dest="verbose",  default=1,  type="int",    help="Verbosity level (0 = quiet, 1 = verbose, 2+ = more)")

    # Read options and args
    (options,args) = parser.parse_args()

    ###############################################################
    ### Start running
    ###############################################################

    tree = getTreeFromFiles(options)
    anaTree(tree,options)
