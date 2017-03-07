#!/usr/bin/env python2

import glob, sys, os
import array as ar
import numpy as np
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


def makeHist(values,hname = "hist"):
    '''
    nbins = len(values)/100
    xmin = min(values)
    xmax = max(values)
    hist = rt.TH1F(hname,hname,nbins,xmin,xmax)
    '''
    hist = rt.TH1F(hname,hname,400,-200,200)
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

    dxs = []
    dys = []
    dzs = []

    if opts.maxEntries == -1: opts.maxEntries = tree.GetEntries()
    if opts.verbose > 0: print("#Going to analyze %i entries" %opts.maxEntries)

    for ientry, entry in enumerate(tree):
        if ientry > opts.maxEntries-1: break
        if opts.verbose > 1:
            print ("#Analyzing entry number %i" %ientry)

        ##############
        ## Do analysis
        ##############

        ## Access particles
        for particle in entry.particles:

            # require initial particle
            if particle.gen < 1: continue
            # require particle to reach EE
            if not particle.reachedEE: continue

            ## Get PosZ
            # convert vector to array
            posx =  ar.array('d',particle.posx)
            posy =  ar.array('d',particle.posy)
            posz =  ar.array('d',particle.posz)
            #for pz in posz: print(pz)
            if len(posx) == 0:
                if opts.verbose > 0: print("Particle with no posx filled")
                continue
            if len(posy) == 0:
                if opts.verbose > 0: print("Particle with no posy filled")
                continue
            if len(posz) == 0:
                if opts.verbose > 0: print("Particle with no posz filled")
                continue

            ## Access cluster2d
            for clust2d in entry.cluster2d:
                # filter cluster2d by energy
                if clust2d.energy < 0.1: continue
                # ignore cluster2d after EE
                if clust2d.layer > 28: continue
                # filter cluster2d from correct side
                if clust2d.z * particle.eta < 0: continue

                # get corresponding posz
                if len(posz) < clust2d.layer:
                    if opts.verbose > 0: print("Problem! No posz for layer %i" %clust2d.layer )
                    continue

                # save values
                dx = posx[clust2d.layer-1]-clust2d.x; dxs.append(dx)
                dy = posy[clust2d.layer-1]-clust2d.y; dys.append(dy)
                dz = posz[clust2d.layer-1]-clust2d.z; dzs.append(dz)

                #dxlayers[clust2d.layer-1].append(dx)

    hists = []
    hists.append(makeHist(dxs,"dx"))
    hists.append(makeHist(dys,"dy"))
    hists.append(makeHist(dzs,"dz"))

    canv = plotHists(hists,"dxyz","delta(Particle pos, clust2d pos)")

    if not opts.batch:
        canv.Draw()
        q = raw_input("Exit")

    canv.SaveAs(opts.plotdir+canv.GetName()+".pdf")

    # save output to root file
    ofile = rt.TFile(opts.plotdir+"/clust2d_plots.root","recreate")
    for obj in hists+[canv]:
        obj.Write()
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
