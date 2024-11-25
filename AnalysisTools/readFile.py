#!/usr/bin/env python3



import sys, os, pwd
import optparse, shlex, re
import math
from ROOT import *
import ROOT
from array import array
from numpy import sqrt


MiniTreeFile = ROOT.TFile.Open('../UFCSCRootMaker/CSC_UF_Ntuple_SegmentAlgoUF.root')
MiniTreeFile.cd()
tree = MiniTreeFile.Get("cscRootMaker/CSCTree")


for i in range( tree.GetEntries() ):
    tree.GetEntry(i)

    print(" Event :  ", tree.Event, "   nSegments  ", tree.cscSegments_nSegments)
    for s in range(tree.cscSegments_nSegments):

        segmentRing          = tree.cscSegments_ID_ring[s]
        segmentStation       = tree.cscSegments_ID_station[s]
        segmentEndcap        = tree.cscSegments_ID_endcap[s]
        segmentChamber       = tree.cscSegments_ID_chamber[s]
        print('  The chamber   ', '        E:',segmentEndcap,'S:',segmentStation,'R:',segmentRing,'C:',segmentChamber)
        print('------- seg # :', s, ' X / Y: ', tree.cscSegments_localX[s], '  /  ', tree.cscSegments_localY[s] )
