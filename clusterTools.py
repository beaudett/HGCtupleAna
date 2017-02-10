#!/usr/bin/env python2
import math
import ROOT as rt
rt.gStyle.SetPalette(rt.kRainBow)

class hitpoint:
    ## Simple class for yield,error storing (instead of tuple)

    def __init__(self, rechit):
        self.x = rechit.x
        self.y = rechit.y
        self.z = rechit.z
        self.layer = rechit.layer
        self.energy = rechit.energy
        self.rho = 0.
        self.dmin = 0.

    # func that is called with print
    def __repr__(self):
        return "(rho: %i, dmin: %f)" %(self.rho,self.dmin)

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

        mgrXY.Add(grXY)
        mgrEneX.Add(grEneX)
        mgrEneY.Add(grEneY)
        mgrRhoDmin.Add(grRhoDmin)

    #print("Filled %i graphs into multigr" %mgrRhoDmin.GetListOfGraphs().GetSize())
    cname = "Hits"
    canv = rt.TCanvas(cname,cname,800,800)
    canv.DivideSquare(4,0.01,0.01)

    canv.cd(4); mgrEneX.Draw("apg off")
    canv.cd(2); mgrXY.Draw("apg off")
    canv.cd(3); mgrRhoDmin.Draw("ap goff")
    canv.cd(1); mgrEneY.Draw("apg off")

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
    newhits = [hitpoint(hit) for hit in hits]

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

    # calc density for each point
    for i,hit1 in enumerate(newhits):
        rho = 0

        for j,hit2 in enumerate(newhits):
            #if i == j: continue
            #if hit2.energy < minE: continue
            dist = math.hypot(hit1.x-hit2.x,hit1.y-hit2.y)
            if dist < dcut:
                chi = hit2.energy
                rho += chi

        hit1.rho = rho

    # calc distance parameter (min dist to hit with higher density)
    for i,hit1 in enumerate(newhits):
        #if hit1.energy < minE: continue

        dmin = 9999
        dmax = 0
        for j,hit2 in enumerate(newhits):
            #if i == j: continue
            #if hit2.energy < minE: continue
            if hit2.rho > hit1.rho:
                dist = math.hypot(hit1.x-hit2.x,hit1.y-hit2.y)
                if dist < dmin: dmin = dist
                if dist > dmax: dmax = dist

        if dmin == 9999: dmin = dmax
        hit1.dmin = dmin

    #print(len(hits),len(newhits))
    #print(newhits[:50])
    #plotRhoDmin(newhits)

    maxrho = max([hit.rho for hit in newhits])
    minrho = maxrho/12
    mindmin = dcut*1.5

    #canv = plotColRhoDmin(newhits,minrho,mindmin)
    canv = plotColRhoDminEneXY(newhits,minrho,mindmin)
    #canv = plotRhoDmin3D(newhits)

    return canv
