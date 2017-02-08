#!/usr/bin/env python2

import math

def calcDensity(hits, dcut = 0.1):
    "Calculate local density based on cutoff distance dcut"

    newhits = []

    for hit1 in hits:
        rho = 0

        for hit2 in hits:
            dist = math.hypot(hit1.x-hit2.x,hit1.x-hit2.x)
            if dist < dcut: rho+=1

        hist.rho = rho
