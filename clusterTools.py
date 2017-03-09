#!/usr/bin/env python2
import math
import numpy as np
import ROOT as rt
rt.gStyle.SetPalette(rt.kRainBow)

class hitpart:

    def __init__(self, genpart):
        self.gen = genpart.gen
        self.eta = genpart.eta
        self.posx = np.array(genpart.posx)
        self.posy = np.array(genpart.posy)
        self.posz = np.array(genpart.posz)
        self.reachedEE = genpart.reachedEE

class hitpoint:
    ## simplified rechit with extra parameters

    def __init__(self, rechit):
        self.x = rechit.x
        self.y = rechit.y
        self.z = rechit.z
        self.layer = rechit.layer
        self.energy = rechit.energy
        self.rho = 0.
        self.dmin = 0.
        self.clustId = -1
        self.centerId = -1

    # func that is called with print
    def __repr__(self):
        return "(rho: %f, dmin: %f)" %(self.rho,self.dmin)

class hitcluster:
    ## simple cluster class

    def __init__(self, hit):
        self.x = hit.x
        self.y = hit.y
        self.z = hit.z
        self.layer = hit.layer
        self.energy = hit.energy
        if hasattr(hit,"nhitAll"): # for cmssw clusters
            self.nhits = hit.nhitAll
        else:
            self.nhits = 1
        #self.hits = [hit] # store all corresp hits

    # func that is called with print
    def __repr__(self):
        #return "(rho: %f, dmin: %f)" %(self.rho,self.dmin)
        return "(x: %f, y: %f, E: %f, nhits: %i )" %(self.x, self.y, self.energy, self.nhits)

    # add single hit
    def addhit(self,hit):
        self.energy += hit.energy
        self.nhits += 1
        #self.hits.append(hit)

def plotRhoDmin(hits):

    grXY = rt.TGraph()
    grRhoDmin = rt.TGraph()

    #print("Going to fill %i points into graph" % len(hits))
    for i,hit in enumerate(hits):
        grXY.SetPoint(i,hit.x,hit.y)
        grRhoDmin.SetPoint(i,hit.rho,hit.dmin)

    #print("Filled %i points into graph" %grRhoDmin.GetN())
    cname = "RhoDmin"
    canv = rt.TCanvas(cname,cname,1000,600)
    canv.Divide(2)
    canv.cd(1)
    grXY.Draw("ap*")
    canv.cd(2)
    grRhoDmin.Draw("ap*")

    return canv

def plotRhoDmin3D(hits):

    grXY = rt.TGraph2D(); rt.SetOwnership(grXY,0)
    grRhoDmin = rt.TGraph2D(); rt.SetOwnership(grRhoDmin,0)

    #print("Going to fill %i points into graph" % len(hits))
    for i,hit in enumerate(hits):
        grXY.SetPoint(i,hit.x,hit.y,hit.energy)
        grRhoDmin.SetPoint(i,hit.rho,hit.dmin,hit.energy)

    grXY.SetMarkerStyle(33)
    grRhoDmin.SetMarkerStyle(33)
    grXY.SetMarkerSize(4.0)
    grRhoDmin.SetMarkerSize(4.0)

    #print("Filled %i points into graph" %grRhoDmin.GetN())

    rt.gStyle.SetPalette(1)
    cname = "RhoDmin3D"
    canv = rt.TCanvas(cname,cname,1000,600)
    canv.Divide(2)
    canv.cd(1)
    grXY.Draw("pcol")
    canv.cd(2)
    grRhoDmin.Draw("pcol")

    return canv

def plotColRhoDmin(hits, minrho = 0.01, mindmin = 2):

    mgrXY = rt.TMultiGraph(); rt.SetOwnership(mgrXY,0)
    mgrRhoDmin = rt.TMultiGraph(); rt.SetOwnership(mgrRhoDmin,0)
    #print("Going to fill %i points into graph" % len(hits))
    for i,hit in enumerate(hits):

        grXY = rt.TGraph(1)
        grRhoDmin = rt.TGraph(1)

        grXY.SetPoint(0,hit.x,hit.y)
        grRhoDmin.SetPoint(0,hit.rho,hit.dmin)

        col = rt.gStyle.GetColorPalette(i * 255/len(hits))
        grXY.SetMarkerColor(col)
        grRhoDmin.SetMarkerColor(col)

        grXY.SetMarkerStyle(20)
        grRhoDmin.SetMarkerStyle(20)

        grXY.SetMarkerSize(1.0)
        grRhoDmin.SetMarkerSize(1.0)

        mgrXY.Add(grXY)
        mgrRhoDmin.Add(grRhoDmin)

        # dublicate possible cluster center hits
        if hit.rho > minrho and hit.dmin > mindmin:
            grXYfilt = rt.TGraph(1)
            grXYfilt.SetPoint(0,hit.x,hit.y)
            grXYfilt.SetMarkerColor(col)
            grXYfilt.SetMarkerStyle(24)
            grXYfilt.SetMarkerSize(2.5)
            mgrXY.Add(grXYfilt)

            # modify marker of rho/dmin graph
            grRhoDmin.SetMarkerStyle(24)
            grRhoDmin.SetMarkerSize(2.5)


    #print("Filled %i graphs into multigr" %mgrRhoDmin.GetListOfGraphs().GetSize())
    cname = "Hits"
    canv = rt.TCanvas(cname,cname,1000,600)
    canv.Divide(2)
    canv.cd(1)
    mgrXY.Draw("apg off")
    canv.cd(2)
    mgrRhoDmin.Draw("ap goff")

    mgrXY.GetXaxis().SetTitle("hit x")
    mgrXY.GetYaxis().SetTitle("hit y")

    mgrRhoDmin.GetXaxis().SetTitle("#rho")
    mgrRhoDmin.GetYaxis().SetTitle("#delta")

    return canv

def plotColRhoDminEneXY(hits, minrho = 0.01, mindmin = 2):

    mgrXY = rt.TMultiGraph(); rt.SetOwnership(mgrXY,0)
    mgrRhoDmin = rt.TMultiGraph(); rt.SetOwnership(mgrRhoDmin,0)

    mgrEneX = rt.TMultiGraph(); rt.SetOwnership(mgrEneX,0)
    mgrEneY = rt.TMultiGraph(); rt.SetOwnership(mgrEneY,0)

    # circles around cluster centers
    clustrs = []

    #print("Going to fill %i points into graph" % len(hits))
    for i,hit in enumerate(hits):

        grXY = rt.TGraph(1)
        grRhoDmin = rt.TGraph(1)
        grEneX = rt.TGraph(1)
        grEneY = rt.TGraph(1)
        grs = [grXY,grRhoDmin,grEneX,grEneY]

        grXY.SetPoint(0,hit.x,hit.y)
        grEneX.SetPoint(0,hit.x,hit.energy)
        grEneY.SetPoint(0,hit.energy,hit.y)
        grRhoDmin.SetPoint(0,hit.rho,hit.dmin)

        col = rt.gStyle.GetColorPalette(i * 255/len(hits))
        for gr in grs:
            gr.SetMarkerColor(col)
            gr.SetMarkerStyle(20)
            gr.SetMarkerSize(1.0)

        # highlight possible cluster center hits with big round marker
        if hit.rho > minrho and hit.dmin > mindmin:
            for gr in grs:
                gr.SetMarkerStyle(24)
                gr.SetMarkerSize(2.5)

            # add circle around cluster centre
            clustr = rt.TEllipse(hit.x,hit.y,mindmin,mindmin)
            clustr.SetLineColor(col)
            clustr.SetFillStyle(0)
            rt.SetOwnership(clustr,0)
            clustrs.append(clustr)

        mgrXY.Add(grXY)
        mgrEneX.Add(grEneX)
        mgrEneY.Add(grEneY)
        mgrRhoDmin.Add(grRhoDmin)

    #print("Filled %i graphs into multigr" %mgrRhoDmin.GetListOfGraphs().GetSize())
    cname = "Hits"
    canv = rt.TCanvas(cname,cname,800,800)
    canv.DivideSquare(4,0.001,0.001)

    canv.cd(3); mgrEneX.Draw("apg off")
    canv.cd(1); mgrXY.Draw("apg off")
    # draw also cluster centers
    #for clustr in clustrs: clustr.Draw()
    canv.cd(4); mgrRhoDmin.Draw("ap goff")
    canv.cd(2); mgrEneY.Draw("apg off")

    mgrXY.GetXaxis().SetTitle("hit x")
    mgrXY.GetYaxis().SetTitle("hit y")
    mgrEneX.GetXaxis().SetTitle("hit x")
    mgrEneX.GetYaxis().SetTitle("hit energy")
    mgrEneY.GetXaxis().SetTitle("hit energy")
    mgrEneY.GetYaxis().SetTitle("hit y")

    mgrRhoDmin.GetXaxis().SetTitle("#rho")
    mgrRhoDmin.GetYaxis().SetTitle("#delta")

    return canv

def calcDensity(hits, dcut = 2.0, minE = 0.01, layer = 10):
    "Calculate local density based on cutoff distance dcut"

    # resave hits in list
    newhits = hits#[hitpoint(hit) for hit in hits]

    # filter by energy
    newhits = [hit for hit in newhits if hit.energy > minE]

    # filter by layer
    newhits = [hit for hit in newhits if hit.layer == layer]

    # filter by +/- z
    newhits = [hit for hit in newhits if hit.z < 0]

    # sort by energy
    #newhits = sorted(newhits, key = lambda h: h.energy)
    # sort by dist from center
    newhits = sorted(newhits, key = lambda h: math.hypot(hit.x,hit.y))

    # switch back to numpy list
    newhits = np.array(newhits)

    # calc density for each point
    for i,hit1 in enumerate(newhits):
        rho = 0

        for j,hit2 in enumerate(newhits):
            #if i == j: continue
            dist = math.hypot(hit1.x-hit2.x,hit1.y-hit2.y)
            if dist < dcut:
                chi = hit2.energy
                rho += chi

        hit1.rho = rho

    # calc distance parameter (min dist to hit with higher density)
    dmax = -1

    for i,hit1 in enumerate(newhits):
        dmin = 9999
        centerId = -1

        for j,hit2 in enumerate(newhits):
            dist = math.hypot(hit1.x-hit2.x,hit1.y-hit2.y)
            if dist > dmax: dmax = dist
            if hit2.rho > hit1.rho:
                if dist < dmin:
                    dmin = dist
                    centerId = j

        if dmin == 9999: dmin = dmax
        hit1.dmin = dmin
        # store index to closest hit with higher density
        hit1.centerId = centerId

    #print(len(hits),len(newhits))
    #print(newhits[:50])
    #plotRhoDmin(newhits)

    return newhits

def makeClustPlots(hits, dcut):

    maxrho = max([hit.rho for hit in hits])
    minrho = maxrho/12
    mindmin = dcut #*1.5

    #canv = plotColRhoDmin(newhits,minrho,mindmin)
    canv = plotColRhoDminEneXY(hits,minrho,mindmin)
    #canv = plotRhoDmin3D(newhits)

    return canv

def makeClusters(hits, minrho , mindmin):

    # make initial cluster seed from cluster centers
    clusters = [hitcluster(hit) for hit in hits if (hit.rho > minrho and hit.dmin > mindmin)]
    #clusterIds = set([i for i,hit in enumerate(hits) if (hit.rho > minrho and hit.dmin > mindmin)])
    #print clusterIds

    for i,hit in enumerate(hits):

        # ignore centers
        if (hit.rho > minrho and hit.dmin > mindmin): continue
        # ignore outliers
        if (hit.dmin > mindmin): continue

        # loop over cluster centers to find closest
        min_dist = 9999
        closest_clustId = -1
        for j,clust in enumerate(clusters):
            dist = math.hypot(clust.x-hit.x,clust.y-hit.y)
            if dist < min_dist:
                min_dist = dist
                closest_clustId = j

        # add hit to closest cluster
        clusters[closest_clustId].addhit(hit)

    return clusters

def compareClusters(clusters_old, clusters_new, layer = 10):

    minE = 0.1
    minNhit = 3
    print("Old clusters:")
    for clust in sorted(clusters_old, key = lambda c: math.hypot(c.x,c.y)):
        if clust.layer != layer: continue
        if clust.energy < minE: continue
        if clust.nhits < minNhit : continue
        print clust,
    print

    print("New clusters:")
    for clust in sorted(clusters_new, key = lambda c: math.hypot(c.x,c.y)):
        if clust.layer != layer: continue
        if clust.energy < minE: continue
        if clust.nhits < minNhit : continue
        print clust,
    print
