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

def addParticleFootprints(canv, layer, particles):

    grPartXY = rt.TGraph(); rt.SetOwnership(grPartXY,0)
    grPartXY.SetMarkerStyle(34)
    grPartXY.SetMarkerColor(rt.kBlack)
    grPartXY.SetMarkerSize(2)

    indx = 0
    #print(particles)
    for particle in particles:
        # require initial particle
        if particle.gen < 1: continue
        # require particle to reach EE
        if not particle.reachedEE: continue
        # select hemisphere
        if particle.eta > 0: continue

        ## Get PosZ
        posx =  particle.posx[layer-1]
        posy =  particle.posy[layer-1]

        #print("Particle position: %f, %f" %(posx,posy))
        grPartXY.SetPoint(indx,posx,posy)

        indx+=1

        #print grPartXY.GetN(), indx

    canv.cd(1)
    grPartXY.Draw("p same")
    rt.gPad.Update()

    return 1

def addClusterFootprints(canv, layer, clusters):

    grClustXY = rt.TGraph(); rt.SetOwnership(grClustXY,0)
    grClustXY.SetMarkerStyle(25)
    grClustXY.SetMarkerColor(rt.kBlack)
    grClustXY.SetMarkerSize(2)

    indx = 0
    #print(particles)
    for cluster in clusters:
        if cluster.z > 0: continue
        if cluster.layer != layer: continue

        grClustXY.SetPoint(indx,cluster.x,cluster.y)
        indx+=1

    canv.cd(1)
    grClustXY.Draw("p same")
    rt.gPad.Update()

    return 1

def anaTree(tree, opts):
    "Analyze entry"

    if opts.maxEntries == -1: opts.maxEntries = tree.GetEntries()
    if opts.verbose > 0: print("#Going to analyze %i entries" %opts.maxEntries)

    # save output to root file
    ofile = rt.TFile(opts.plotdir+"/clust_plots.root","recreate")

    # accumulate rechits from different events
    accumhits = np.array([])
    particles = []
    clusters = []

    for ientry, entry in enumerate(tree):
        if ientry > opts.maxEntries-1: break
        if opts.verbose > 1:
            print ("#Analyzing entry number %i" %ientry)

        #if ientry != 44: continue

        ##############
        ## Do analysis
        ##############
        minE = 0.01
        minLayer = 10
        maxLayer = 11

        if False:
            # accumulate hits
            accumhits = np.append(accumhits, [hitpoint(rechit) for rechit in entry.rechits_raw if rechit.energy > minE])
            particles += [hitpart(particle) for particle in entry.particles if particle.gen > 0]
            clusters += [hitpoint(cluster) for cluster in entry.cluster2d if cluster.energy > minE]
        else:
            # don't accum hits
            accumhits = np.array([hitpoint(rechit) for rechit in entry.rechits_raw if rechit.energy > minE])
            particles = [hitpart(particle) for particle in entry.particles if particle.gen > 0]
            clusters = [hitpoint(cluster) for cluster in entry.cluster2d if cluster.energy > minE]

        if opts.verbose > 0:
            print("Accumulated %i hits, %i 2d clusters, %i particles" % (len(accumhits),len(clusters),len(particles)))

        for layer in range(minLayer,maxLayer):

            for step in range(2,3):
                dcut = 0+step*1.0

                # create output dir
                ofdir = "Event_%i" %ientry + "/Layer_%i" %layer + "/"
                ofile.mkdir(ofdir)
                ofile.cd(ofdir)

                #canv = calcDensity(entry.rechits, dcut, minE, layer)
                canv = calcDensity(accumhits, dcut, minE, layer)
                canv.SetName(canv.GetName()+"_ev%i_ly%i_dc%0.1f"% (ientry, layer,dcut) )

                addParticleFootprints(canv, layer, particles)
                addClusterFootprints(canv, layer, clusters)

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
