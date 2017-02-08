#!/usr/bin/env python2

import glob, sys, os
import array as ar
import numpy as np
import ROOT as rt
from clusterTools import *

def getTreeFromFiles(opts):

    chain = rt.TChain(opts.treename)
    for fname in [opts.infile]:
        if opts.verbose > 0: print("#Adding file %s" %fname)

        chain.Add(fname)

    if opts.verbose > 0:
        print("#Added %i files to chain %s" %(chain.GetNtrees(),opts.treename))
        print("#Chain contains %i entries" % chain.GetEntries())

    return chain


def makeHist(values,hname = "hist"):

    nbins = len(values)/100
    xmin = min(values)
    xmax = max(values)
    hist = rt.TH1F(hname,hname,nbins,xmin,xmax)

    for val in values: hist.Fill(val)

    return hist

def plotHists(hists,cname="hists",ctitle="hists"):

    canv = rt.TCanvas(cname,cname,600,600)
    canv.DivideSquare(len(hists),0.01,0.01)

    plotopt = "hist"
    for i,hist in enumerate(hists):
        canv.cd(i+1)
        hist.Draw(plotopt)

    return canv

def anaTree(tree, opts):
    "Analyze entry"

    if opts.maxEntries == -1: opts.maxEntries = tree.GetEntries()
    if opts.verbose > 0: print("#Going to analyze %i entries" %opts.maxEntries)

    # save output to root file
    ofile = rt.TFile(opts.plotdir+"/clust_plots.root","recreate")

    for ientry, entry in enumerate(tree):
        if ientry > opts.maxEntries-1: break
        if opts.verbose > 1:
            print ("#Analyzing entry number %i" %ientry)

        #if ientry != 44: continue

        ##############
        ## Do analysis
        ##############
        minE = 0.1
        minLayer = 4
        maxLayer = 15

        for layer in range(minLayer,maxLayer):

            for step in range(1,11):
                dcut = 0+step*1.0

                # create output dir
                #ofdir = "Event_%i" %ientry + "/Layer_%i" %layer + "/dcut_%0.1f" %dcut
                ofdir = "Event_%i" %ientry + "/Layer_%i" %layer + "/"
                ofile.mkdir(ofdir)
                ofile.cd(ofdir)

                canv = calcDensity(entry.rechits, dcut, minE, layer)
                #canv.SetName(canv.GetName()+"_eve%i"%ientry+"_%0.1f"%dcut )
                canv.SetName(canv.GetName()+"_ev%i_ly%i_dc%0.1f"% (ientry, layer,dcut) )

                if not opts.batch:
                    canv.Draw()
                    q = raw_input("exit")

                canv.Write()
    ofile.Close()

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
    parser.add_option("-n","--maxEve","--maxEvents", "--maxEntries",dest="maxEntries", default=-1,  type="int",    help="Maximum entries to analyze")
    parser.add_option("-v","--verbose",  dest="verbose",  default=1,  type="int",    help="Verbosity level (0 = quiet, 1 = verbose, 2+ = more)")

    # File options
    parser.add_option("-f","--file", dest="infile",default="root://eosuser.cern.ch//eos/user/b/beaudett/local/hgcal/singleElePt15-extendedNtp.root",help="Input File name(s)")
    parser.add_option("-t","--tname", dest="treename",default="ana/hgc",help="Tree name")

    # Plot options
    parser.add_option("-p","--pdir","--plotdir", dest="plotdir",default="plots/",help="Plot directory")

    # Read options and args
    (options,args) = parser.parse_args()

    # evaluate options
    if not os.path.exists(options.plotdir): os.makedirs(options.plotdir)

    ###############################################################
    ### Start running
    ###############################################################

    print("Starting")
    tree = getTreeFromFiles(options)
    anaTree(tree,options)

    print("Done")
