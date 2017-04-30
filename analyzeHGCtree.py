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

def createHistogramDictionnary():
    histDict = {}
    histDict['fbrem'] = rt.TH1F("fbrem","fbrem",110,-0.05,1.05)
    histDict['dx'] =  rt.TH1F("dx","dx",100,-0.25,0.25)
    histDict['dy'] =  rt.TH1F("dy","dy",100,-0.25,0.25)
    histDict['dz'] =  rt.TH1F("dz","dz",100,-0.25,0.25)
    histDict['detapos'] = rt.TH1F("detapos","detapos",100,-0.1,0.1)
    histDict['detaneg'] = rt.TH1F("detaneg","detaneg",100,-0.1,0.1)
    histDict['dphipos'] = rt.TH1F("dphipos","dphipos",100,-0.1,0.1)
    histDict['dphineg'] = rt.TH1F("dphineg","dphineg",100,-0.1,0.1)
    histDict['dphiposfbrem'] = rt.TH2F("dphiposfbrem","dphiposfbrem",25,0,1,100,-0.1,0.1)
    histDict['dphinegfbrem'] = rt.TH2F("dphinegfbrem","dphinegfbrem",25,0,1,100,-0.1,0.1)
    return histDict

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
    histDict = createHistogramDictionnary()
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
            histDict["fbrem"].Fill(particle.fbrem)

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

            ## Access rechits
            for rechit in entry.rechits_raw:
                # filter rechits by energy
                if rechit.energy < 0.1: continue
                # ignore rechits after EE
                if rechit.layer > 28: continue
                # filter rechits from correct side
                if rechit.z * particle.eta < 0: continue

                # get corresponding posz
                if len(posz) < rechit.layer:
                    if opts.verbose > 0: print("Problem! No posz for layer %i" %rechit.layer )
                    continue

                # save values
                dx = posx[rechit.layer-1]-rechit.x; histDict["dx"].Fill(dx)
                dy = posy[rechit.layer-1]-rechit.y; histDict["dy"].Fill(dy)
                dz = posz[rechit.layer-1]-rechit.z; histDict["dz"].Fill(dz)
                posextrap = rt.TVector3(posx[rechit.layer-1],posy[rechit.layer-1],posz[rechit.layer-1])
                posrechit = rt.TVector3(rechit.x,rechit.y,rechit.z)
                deta=posextrap.Eta()-posrechit.Eta()
                dphi=posextrap.Phi()-posrechit.Phi()
                if particle.pid>0:
                    histDict['detapos'].Fill(deta)
                    histDict['dphipos'].Fill(dphi)
                    histDict['dphiposfbrem'].Fill(particle.fbrem,dphi)
                else:
                    histDict['detaneg'].Fill(deta)
                    histDict['dphineg'].Fill(dphi)
                    histDict['dphinegfbrem'].Fill(particle.fbrem,dphi)
            

    if not opts.batch:
        canv.Draw()
        q = raw_input("Exit")


    # save output to root file
    ofile = rt.TFile(opts.plotdir+"/plots.root","recreate")
    for obj in histDict:
        histDict[obj].Write()
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
