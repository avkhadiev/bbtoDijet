#! /usr/bin/env python
import os
import glob
import math
import array
import numpy as np
import sys
import time
import string
import ROOT


# import tdrstyle
# tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadTopMargin(0.09);
ROOT.gStyle.SetPadLeftMargin(0.12);
ROOT.gStyle.SetPadRightMargin(0.25);
ROOT.gStyle.SetPaintTextFormat("1.1f");
ROOT.gROOT.SetBatch()
# see https://root.cern.ch/phpBB3/viewtopic.php?t=15076	ROOT.gEnv.GetValue("Canvas.SavePrecision", -1)
ROOT.gEnv.GetValue("Canvas.SavePrecision", -1)
ROOT.gEnv.SetValue('Canvas.SavePrecision', "15")
# print sys.float_info

probes = ["HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_v2",
     "HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160_v2"
    ]
tags = ["HLT_DoubleJetsC100_p014_DoublePFJetsC100MaxDeta1p6_v2",
   # "HLT_Mu50_v3"#, # "HLT_TkMu50_v1",
     "HLT_DoubleJetsC100_p026_DoublePFJetsC160_v2"#,
   # "HLT_AK8PFJet140_v1",
   # "HLT_AK8PFJet200_v1",
   # "HLT_DiPFJetAve80_v3",
   # "HLT_PFJet80_v5"
   ]

# use leading and sub-leading jetss
jets = [0, 1]

        # PARAMETER NAME |     AXIS TITLE    | BINS,   MIN,   MAX
parameters = {'bTagCSVOnline':      ('Online CSV',            100,  0, 1.0     ),
        'bTagCSVOffline':          ('Offline CSV',              100,  0, 1.0    )#,
      #    'pfjetpt':         ('PFJet p_{T} (GeV)', 30,  0, 300    )#,
      #    'pfjetmass':        ('PFJet mass (GeV)',  24    0, 1200
      #   'pfjeteta':         ('PFJet eta',         10,   0, 5      )
         }


# open file
input_file = ROOT.TFile('/eos/uscms/store/user/aavkhadi/JetHT/eff_bTagDijetV11/hlt_bTagDijetV11.root')
tree = input_file.Get('hltana/HLTAnalysis')

# speed things up by disabling branches
tree.SetBranchStatus("*", 0)
for parameter in parameters:
    tree.SetBranchStatus(parameter, 1)
for tag in tags:
     tree.SetBranchStatus(tag, 1)
for probe in probes:
    tree.SetBranchStatus(probe, 1)

# offline cuts
# for plotting correlations with kinematics cuts:
tree.SetBranchStatus("pfjeteta", 1)
tree.SetBranchStatus("pfjetpt", 1)

# WANT OPTION PARSER
# WANT Loop over parameters, then over probes, then over tags, then jets

# declareHistograms
# input:  string probe, string tag, string parameter 
#        
# output: dictionary  of TH1F histograms with binning correposponding to the parameters dictonary
def declareHistograms(probe, tag, parameter):
    print 'Declaring histograms in %s for a tag-probe pair %s-%s' % (parameter, tag, probe)
    histograms = {}
    axisTitle, nbins, minbin, maxbin = parameters[parameter]
    for jet in jets:
                               #  HISTOGRAM NAME:  TAG  + PROBE + PARAMETER + JET | AXIS TITLE + JET | BINS,  MIN,   MAX
        numerator   = ROOT.TH1F("%s_%s_%s_pass_%d" % (tag, probe, parameter, jet),"%s -- tag + probe" % (parameter) + ";" + axisTitle + str(jet) + ";" + "tag + probe" + ";", nbins, minbin, maxbin)
        denominator = ROOT.TH1F("%s_%s_%s_all_%d"  % (tag, probe, parameter, jet),"%s -- probe"       % (parameter) + ";" + axisTitle + str(jet) + ";" + "tag"         + ";", nbins, minbin, maxbin) 
	histograms[jet] = (numerator, denominator)
	# debugging -------------------------
	print 'numerator and bins: %s, %d' % ( numerator.GetName(), numerator.GetNbinsX() )
	print 'numerator as written to the dictionary, bins: %s, %d' % (histograms[jet][0].GetName(), histograms[jet][0].GetNbinsX())
        # -----------------------------------
    return histograms

# fillHistograms
# input:  dictionary of histograms, string probe, string tag, string parameter
#         TTree tree -- contains branches with parameters
# output: the dictionary histograms with filled histograms
def fillHistograms(probe, tag,  parameter, tree):
    print 'Filling histograms for a tag-probe pair %s-%s, parameter %s' % (tag, probe, parameter)
    global histograms
    for jet in jets:
        numerator   = histograms[jet][0]
        denominator = histograms[jet][1]
	# cheap trick to get the second max if it exists, max ow
        tree.Draw("%s[%d]>>%s" % (parameter, jet, numerator.GetName()), "(%s == 1) && (%s == 1)" % (tag, probe))
        # debugging ------------------------
        # print ("Entries in %s --- %d" % (numerator.GetName(), numerator.GetEntries()))
        # print ("Entries in %s --- %d" % (histoPairs[jet][0].GetName(), histoPairs[jet][0].GetEntries()))
        # ----------------------------------
        tree.Draw("%s[%d]>>%s" % (parameter, jet, denominator.GetName()), "(%s == 1)"            % (tag       ))
        # debugging ------------------------
        # print ("Entries in %s --- %d" % (denominator.GetName(), denominator.GetEntries()))
        # print ("Entries in %s --- %d" % (histoPairs[jet][1].GetName(), histoPairs[jet][1].GetEntries()))
        # ---------------------------------
        # debugging --------------------------------
        print 'in fillHistograms:'
        print 'numerator is %s'    % (numerator.GetName())
        print 'denominator is %s' % (denominator.GetName())
        #-------------------------------------------
        # Print out numerator
	histocanvas_num   = ROOT.TCanvas(numerator.GetName(), numerator.GetName(), 800, 600) 
        histocanvas_num.SetFillColor( 19 ) 
        xAxis_num = numerator.GetXaxis()
        yAxis_num = numerator.GetYaxis()
        xAxis_num.SetTitle(parameters[parameter][0])
        yAxis_num.SetTitle("Accepted")
        numerator.Draw("LP")
        # CMS text
        CMSLine="CMS"
        CP=ROOT.TLatex(0.12,0.92, CMSLine)
        CP.SetNDC(ROOT.kTRUE)
        CP.SetTextSize(0.05)
        CP.Draw()
        # Lumi
        CMSLineLumi="#sqrt{s}=13 TeV"
        CP1=ROOT.TLatex(0.67,0.92, CMSLineLumi)
        CP1.SetNDC(ROOT.kTRUE)
        CP1.SetTextSize(0.04)
        CP1.Draw()
        # ExtraText
        CMSLineExtra="#bf{#it{Preliminary}}"
        CP2=ROOT.TLatex(0.195,0.92, CMSLineExtra)
        CP2.SetNDC(ROOT.kTRUE)
        CP2.SetTextSize(0.04)
        CP2.Draw()
        histocanvas_num.SaveAs("output/histos/" + "histo_" + numerator.GetName() + ".pdf")
        # Print out denominator
	histocanvas_denom = ROOT.TCanvas(denominator.GetName(), denominator.GetName(), 800, 600)
        histocanvas_denom.SetFillColor( 19 ) 
        xAxis_denom = denominator.GetXaxis()
        yAxis_denom = denominator.GetYaxis()
        xAxis_denom.SetTitle(parameters[parameter][0])
        yAxis_denom.SetTitle("All")
        denominator.Draw("LP")
        # CMS text
        CMSLine="CMS"
        CP=ROOT.TLatex(0.12,0.92, CMSLine)
        CP.SetNDC(ROOT.kTRUE)
        CP.SetTextSize(0.05)
        CP.Draw()
        # Lumi
        CMSLineLumi="#sqrt{s}=13 TeV"
        CP1=ROOT.TLatex(0.67,0.92, CMSLineLumi)
        CP1.SetNDC(ROOT.kTRUE)
        CP1.SetTextSize(0.04)
        CP1.Draw()
        # ExtraText
        CMSLineExtra="#bf{#it{Preliminary}}"
        CP2=ROOT.TLatex(0.195,0.92, CMSLineExtra)
        CP2.SetNDC(ROOT.kTRUE)
        CP2.SetTextSize(0.04)
        CP2.Draw()
        histocanvas_denom.SaveAs("output/histos/" + "histo_" +  denominator.GetName() + ".pdf")         
    return histograms

# TGraphAsymmErrors makeEffGraph
# input:  TH1F numerator, denominator histograms
# output: TGraphAsymmErrors corresponding turn-on curve
def makeEffGraph(numerator, denominator):
    # binning of numerator and denominator has to match!
    # x[ibin] of numerator and denominator is equal for each ibin!
    nbins = denominator.GetNbinsX()
    print nbins
    npass =  ROOT.TVector(nbins)
    ntotal =  ROOT.TVector(nbins)
    x =  ROOT.TVector(nbins)
    y =  ROOT.TVector(nbins)
    errxLow =  ROOT.TVector(nbins)
    errxHigh =  ROOT.TVector(nbins)
    erryLow =  ROOT.TVector(nbins)
    erryHigh = ROOT.TVector(nbins)
    for ibin in range(0, nbins):
        # https://root.cern.ch/doc/master/classTH1.html#a3e2be0555e806ae3276f9fcec91865c6
        # bin = 0 is the underflow bin, bin = nbins+1 is the overflow bin
        # set x & y coordinate
        x[ibin]      = float( denominator.GetXaxis().GetBinCenter( ibin + 1)  )
        # print "xbin is %f" % (x[ibin])	
        ntotal[ibin] = float( denominator.GetBinContent(           ibin + 1) )
	npass[ibin]  = float(  numerator.GetBinContent(            ibin + 1) )
	# debugging ----------------------------------------------------
	#if ibin % 10 == 0:
	#	print '------------------------------------'
	#	print 'ibin is %d'      % (ibin)
	#	print 'numerator bin content is %d (%s)'   % (numerator.GetBinContent(ibin+1),   numerator.GetName())
	#	print 'denominator bin content is %d (%s)' % (denominator.GetBinContent(ibin+1), denominator.GetName())
	#	print 'Using ntotal %d' % (ntotal[ibin]) 
	#	print 'Using npass  %d' % (npass[ibin])
	## --------------------------------------------------------------	
	if ntotal[ibin] == 0.0:
             # print 'ntotal is 0'
	     y[ibin] = 0.0
        else: 
             # print 'y[ibin] is %.2f'%(npass[ibin]/ntotal[ibin])
             y[ibin] = npass[ibin]/ntotal[ibin]
             # print '  --> y[ibin] is %.2f'%(y[ibin])
        errxLow[ibin] = 0.0
        # set x & y coordinate errors
        errxLow[ibin]  =  0.0
        errxHigh[ibin] = 0.0
        if y[ibin] == 0.0:
            erryLow[ibin]  = 0.0
            erryHigh[ibin] = 0.0
        else:
            erryLow[ibin]  = y[ibin] - ROOT.TEfficiency.ClopperPearson( ntotal[ibin], npass[ibin], 0.683, False)
            erryHigh[ibin] =           ROOT.TEfficiency.ClopperPearson( ntotal[ibin], npass[ibin], 0.683, True) - y[ibin]
        # debugging ----------------------------------------------------
	#if ibin % 10 == 0:
	#	print 'y, yerrLow, yerrHigh are:'
	#	print y[ibin], erryLow[ibin], erryHigh[ibin]
	#	print '------------------------------------'
	# --------------------------------------------------------------	
    curve = ROOT.TGraphAsymmErrors(x,y,errxLow, errxHigh, erryLow, erryHigh) 
    return curve

# makeCurves
# input:  dictionary histograms -- contain histograms by jet for a given parameter 
# output  dictionary -- parameters as keys, list of TGraphAssymetricError turn-on curves for 4 leading jets as values
def makeCurves(histograms, parameter):
    print 'Creating turn-on curves'
    curves = {}	
    for jet in jets:
    	numerator   = histograms[jet][0] 
	denominator = histograms[jet][1]
	# debugging ------------------------ 
        print ("Numerator for curve: %s, %d entries"   % (numerator.GetName(), numerator.GetEntries()))
        print ("Denominator for curve: %s, %d entries" % (denominator.GetName(), denominator.GetEntries()))
        # ----------------------------------
        curve = makeEffGraph(numerator, denominator)        
        canvas = ROOT.TCanvas("eff_%s_%s_vs_%s_jet_%d" % (tag, probe, parameter, jet), "%s / %s efficiency vs. %s, %d jet" % (probe, tag, parameter, jet), 800, 600)
        canvas.SetFillColor( 19 )
        canvas.SetGrid()
        curve.SetTitle( "turn-on curve: %s / %s" % ( tag, probe )  )
        curve.SetMarkerColor( 4  )
        curve.SetMarkerStyle( 21 )
        xAxis = curve.GetXaxis()
        yAxis = curve.GetYaxis()
        xAxis.SetTitle(parameters[parameter][0])
        yAxis.SetTitle("Efficiency")
        curve.Draw("ALP")
	# CMS text
	CMSLine="CMS"
	CP=ROOT.TLatex(0.12,0.92, CMSLine)
	CP.SetNDC(ROOT.kTRUE)
	CP.SetTextSize(0.05)
	CP.Draw()
	# Lumi
	CMSLineLumi="#sqrt{s}=13 TeV"
	CP1=ROOT.TLatex(0.67,0.92, CMSLineLumi)
	CP1.SetNDC(ROOT.kTRUE)
	CP1.SetTextSize(0.04)
	CP1.Draw()
	# ExtraText
	CMSLineExtra="#bf{#it{Preliminary}}"
	CP2=ROOT.TLatex(0.195,0.92, CMSLineExtra)
	CP2.SetNDC(ROOT.kTRUE)
	CP2.SetTextSize(0.04)
	CP2.Draw()
        canvas.SaveAs("output/plots/" + canvas.GetName() + ".pdf")
    return
    
# makeCorrelationPlots: tag, onParam, offParam, offlineCut -> TH2D correlation plot onParam vs. offParam for events accepted by tag (+ offlineCut on offParam)
def makeCorrelationPlots(tag, onParam, offParam, offlineCut):
    global tree
    print 'making correlation plots for %s vs. %s' % (onParam, offParam)
    axisTitleX, nbinsX, minbinX, maxbinX = parameters[onParam]
    axisTitleY, nbinsY, minbinY, maxbinY = parameters[offParam]
    for jet in jets:
        print 'working with jet %d' % (jet)
        print 'declaring a 2D histogram...'
        correlation = ROOT.TH2D(("%s_vs_%s_%s_%d") % (onParam, offParam, tag, jet), "%s vs. %s, offline cut %s, [%d]" % (onParam, offParam, offlineCut, jet) + ";" + axisTitleX + ";" + axisTitleY + ";" + "Events", nbinsX, minbinX, maxbinX, nbinsY, minbinY, maxbinY)
        correlation
        print 'filling the histogram %s...' % (correlation.GetName())
        tree.Draw("MaxIf$(%s, %s):Max$(%s)>>%s" % (offParam, offlineCut, onParam, correlation.GetName()), "(%s == 1)" % (tag))
        print 'Printing histograms...'
        canvas = ROOT.TCanvas(correlation.GetName(), correlation.GetTitle(), 800, 600)
        canvas.SetFillColor( 19 )
        # xAxis = correlation.GetXaxis()
        # yAxis = correlation.GetYaxis()
        zAxis = correlation.GetZaxis()
        zAxis.SetRangeUser(10, 20000)
        correlation.Draw("COLZ")
        canvas.SaveAs("output/correlations/" + "corr_" + correlation.GetName() + ".pdf")
    return

# input: string yaxis, xaxis, trigger cut
#        float ymin, ymax, xmin, xmax
#        name of the tree (from where to draw)
# output: 2D histogram for the output 
def plotCorrelation(paramY, paramX, trigger_cut):
    global tree 
    global parameters
    ROOT.gStyle.SetOptStat(0)
    for jet in jets:
	axisTitleY, nbinsY, minbinY, maxbinY = parameters[paramY]
	axisTitleX, nbinsX, minbinX, maxbinX = parameters[paramX]
	histo = ROOT.TH2D("%s_vs_%s_%d_%s" % (paramY, paramX, jet, trigger_cut),"%s vs. %s (%d)" % (paramY, paramX, jet) + ";" + axisTitleX + str(jet) + ";" + axisTitleY + str(jet) + ";" + "Events" + ";", nbinsX, minbinX, maxbinX, nbinsY, minbinY, maxbinY)	
	tree.Draw("%s[%d]:%s[%d]>>%s" % (paramY, jet, paramX, jet, histo.GetName()),"%s == 1" % (trigger_cut))
	# set up canvas
	histocanvas = ROOT.TCanvas(histo.GetName(),"%s vs %s (%d)" % (paramY, paramX, jet), 800, 600) 
        histocanvas.SetFillColor( 19 ) 
        xAxis = histo.GetXaxis()
        yAxis = histo.GetYaxis()
	zAxis = histo.GetZaxis()
        xAxis.SetTitle(parameters[paramX][0])
        yAxis.SetTitle(parameters[paramY][0])
        zAxis.SetTitle("Events")
	zAxis.SetRangeUser(50, 10000)
	histo.Draw("COLZ")
        histocanvas.SaveAs(histocanvas.GetName() + ".pdf")
    return

# loop over T&P pairs
for probe in probes:
	for parameter in parameters:
		for tag in tags:
        		histograms = declareHistograms(probe, tag, parameter)
        		fillHistograms(probe, tag, parameter, tree)
        		makeCurves(histograms, parameter)  

for tag in tags:
	plotCorrelation("bTagCSVOffline", "bTagCSVOnline", tag) 

