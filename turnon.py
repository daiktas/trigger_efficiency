import ROOT
import numpy as np
import os
import multiprocessing
from array import array

import sys
sys.path.insert(0, '/home/hep/vc1117/tdr-style')

import tdrstyle

def draw_cumul(name, x, x_bins, ctau, cut, color, path , mLLP=0, mLSP=0):
    t_MC = ROOT.TChain("Friends")
    t_MC.Add("*ctau-"+ctau+"*")
    MC_after = ROOT.TH1F("MC_after", "MC_after", len(x_bins)-1, x_bins)
    MC_before = ROOT.TH1F("MC_before", "MC_before", len(x_bins)-1, x_bins)
    t_MC.Project("MC_before", x, '(llpinfo_llp_mass == 2000)*(llpinfo_lsp_mass == 200)*({0})*(genWeight)'.format(cut))
    t_MC.Project("MC_after", x, '(llpinfo_llp_mass == 2000)*(llpinfo_lsp_mass == 200)*({0})*(genWeight)*({1})'.format(cut, path))
    ratio = MC_after.Clone("name")
    ratio.Divide(MC_after, MC_before, 1.0, 1.0 ,"B")
    ratio.Sumw2()
    ratio.SetStats(0)
    ratio.SetLineColor(color)
    return ratio

def hadd(files, process):
    string = "hadd -f " + process + ".root "
    for f in files[process]:
        print f
        string = string+f+" "
    print string
    os.system(string)

adding = input('Do you want to hadd (1) or no (0): ')
if adding:

    print "hadding"

    filePath = "friends/"
    files = {}

    for fName in os.listdir(filePath):
        print "reading in process ", fName,"..."
        process = fName.split(".")[0]
        files[process] = []
        print process
        for f in os.listdir(os.path.join(filePath,process)):
            files[process].append(os.path.join(filePath,process,f))
            print f

    jobs = []
    for process in files:
        thread = multiprocessing.Process(target=hadd, args=(files, process))
        jobs.append(thread)

    for j in jobs:
        j.start()

    for j in jobs:
        j.join()

t = ROOT.TChain("Friends")
t.Add("SingleMu*.root")
t_tt = ROOT.TChain("Friends")
t_tt.Add("TT*.root")

def turnon_1D(x, xname, x_bins, path, cut, cut_text, do_fit):
    x_bins = array('d', x_bins)
    c = ROOT.TCanvas(xname, xname)

    before = ROOT.TH1F("before", "before", len(x_bins)-1, x_bins)
    before_tt = ROOT.TH1F("before_tt", "before_tt", len(x_bins)-1, x_bins)

    after = ROOT.TH1F("after", "after", len(x_bins)-1, x_bins)
    after_tt = ROOT.TH1F("after_tt", "after_tt", len(x_bins)-1, x_bins)

    print '{0}*({1})'.format(x, cut)
    t.Project("before", '{0}*({1})'.format(x, cut))
    t_tt.Project("before_tt", '{0}*({1})*(genWeight)'.format(x, cut))

    before.Sumw2()
    before_tt.Sumw2()

    t.Project("after", '{0}*({1})*({2})'.format(x, cut, path))
    t_tt.Project("after_tt", '{0}*({1})*({2})*(genWeight)'.format(x, cut, path))

    after.Sumw2()
    after_tt.Sumw2()

    eff = after.Clone("ratio")
    eff.GetXaxis().SetTitle(xname)
    eff.GetYaxis().SetTitle("efficiency")
    eff.Divide(after, before, 1.0, 1.0 ,"B")
    eff.SetStats(0)
    eff.SetLineColor(ROOT.kOrange)
    eff.SetMarkerColor(ROOT.kOrange)
    eff.GetYaxis().SetRangeUser(0,1.2)
    eff.GetXaxis().SetNdivisions(6)

    eff_tt = after_tt.Clone("ratio")
    eff_tt.GetXaxis().SetTitle(xname)
    eff_tt.GetYaxis().SetTitle("efficiency")
    eff_tt.Divide(after_tt, before_tt, 1.0, 1.0 ,"B")
    eff_tt.SetStats(0)
    eff_tt.SetLineColor(ROOT.kBlue)
    eff_tt.SetMarkerColor(ROOT.kBlue)

    if do_fit:

        eff_fit = ROOT.TF1("eff_fit", "[4] + ([0] - [4])/(1+[3]*TMath::Exp(-[1]*(x-[2])))^([5])", min(x_bins), max(x_bins))
        eff_prefit = ROOT.TF1("eff_prefit", "[0]/(1 + TMath::Exp(-[1]*(x-[2])))", min(x_bins), max(x_bins))
        eff_prefit.SetParameters(1, 0.2, 200)
        eff_fit.SetLineColor(ROOT.kOrange)
        eff.Fit("eff_prefit", "")
        eff_fit.SetParameters(eff_prefit.GetParameter(0), eff_prefit.GetParameter(1), eff_prefit.GetParameter(2), 1, 0, 1)
        eff.Fit("eff_fit", "")

    u_0p001 = draw_cumul("u_0p001", x, x_bins, "0p001", cut, ROOT.kRed-7, path, mLLP = 2000., mLSP=200.)
    u_0p01 = draw_cumul("u_0p01", x, x_bins, "0p01", cut, ROOT.kRed-6, path, mLLP = 2000., mLSP=200.)
    u_0p1 = draw_cumul("u_0p1", x, x_bins, "0p1", cut, ROOT.kRed-5, path, mLLP = 2000., mLSP=200.)
    u_1 = draw_cumul("u_1", x, x_bins, "1", cut, ROOT.kRed-4, path, mLLP = 2000., mLSP=200.)
    u_10 = draw_cumul("u_10", x, x_bins, "10", cut, ROOT.kRed-3, path, mLLP = 2000., mLSP=200.)
    u_100 = draw_cumul("u_100", x, x_bins, "100", cut, ROOT.kRed-2, path, mLLP = 2000., mLSP=200.)
    u_1000 = draw_cumul("u_1000", x, x_bins, "1000", cut, ROOT.kRed-1, path, mLLP = 2000., mLSP=200.)
    u_10000 = draw_cumul("u_10000", x, x_bins, "10000", cut, ROOT.kRed, path, mLLP = 2000., mLSP=200.)

    prompt_c = draw_cumul("prompt_c", x, x_bins, "0p001", cut, ROOT.kViolet+4, path, mLLP = 2000., mLSP=1900.)
    blike_c = draw_cumul("blike_c", x, x_bins, "1", cut, ROOT.kViolet+2, path, mLLP = 2000., mLSP=1900.)
    displaced_c = draw_cumul("displaced_c", x, x_bins, "10000", cut, ROOT.kViolet,path, mLLP = 2000., mLSP=1900.)

    prompt_c.SetLineStyle(9)
    blike_c.SetLineStyle(9)
    displaced_c.SetLineStyle(9)

    c.Draw("")

    eff.Draw("SAME")
    eff_tt.Draw("SAME")
    
    p = ROOT.TLatex()
    p.SetTextSize(0.04)
    p.DrawLatexNDC(0.4, 0.85, cut_text)

    if do_fit:
        eff_fit.Draw("SAME")
        p.DrawLatexNDC(0.4, 0.4, "#varepsilon (200 GeV) = " + str(round(eff_fit.Eval(200), 3)) )
        p.DrawLatexNDC(0.4, 0.3, "#varepsilon (250 GeV) = " + str(round(eff_fit.Eval(250), 3)) )
        p.DrawLatexNDC(0.4, 0.2, "#varepsilon (300 GeV) = " + str(round(eff_fit.Eval(300), 3)) )

    u_0p001.Draw("HISTSAME")
    u_0p01.Draw("HISTSAME")
    u_0p1.Draw("HISTSAME")
    u_1.Draw("HISTSAME")
    u_10.Draw("HISTSAME")
    u_100.Draw("HISTSAME")
    u_1000.Draw("HISTSAME")
    u_10000.Draw("HISTSAME")

    #prompt_c.Draw("HISTSAME")
    #blike_c.Draw("HISTSAME")
    #displaced_c.Draw("HISTSAME")
    
    line = ROOT.TGraph(2)
    line.SetPoint(0, 0, 1)
    line.SetPoint(1, x_bins[-1], 1) 
    line.SetLineColor(ROOT.kRed)
    line.SetLineStyle(9)
    line.Draw("L")
    line2 = ROOT.TLine(300, 0, 300, 1)
    line2.SetLineStyle(2)
    line2.Draw("SAME")

    if not do_fit:
        l = ROOT.TLegend(0.3, 0.30, 0.5, 0.55)
    else:
        l = ROOT.TLegend(0.75, 0.45, 0.95, 0.75)

    l.AddEntry(eff, "SingleMu", "p")
    #l.AddEntry(eff_tt, "t #bar{t}", "p")
    l.AddEntry(u_0p001, "1 #mum, un.", "l")
    l.AddEntry(u_0p01, "10 #mum, un.", "l")
    l.AddEntry(u_0p1, "0.1 mm, un.", "l")
    l.AddEntry(u_1, "1 mm, un.", "l")
    l.AddEntry(u_10, "1 cm, un.", "l")
    l.AddEntry(u_100, "10 cm, un.", "l")
    l.AddEntry(u_1000, "1 m, un.", "l")
    l.AddEntry(u_10000, "10 m, un.", "l")
    #l.AddEntry(prompt_c, "1 #mum, compressed", "l")
    #l.AddEntry(blike_c, "1 mm, compresed", "l")
    #l.AddEntry(displaced_c, "10 m, compressed", "l")
    
    l.SetTextSize(0.05)
    l.Draw()

    tdrstyle.cmsPrel(35900, energy=13, simOnly=False, onLeft=True, sp=0, textScale=1, xoffset=0, thisIsPrelim=True)
    c.SaveAs(x.replace('/', '')+"_"+cut.replace('>','_').replace(' ','')+".pdf")
    c.SaveAs(x.replace('/', '')+"_"+cut.replace('>','_').replace(' ','')+".png")
    c.SaveAs(x.replace('/', '')+"_"+cut.replace('>','_').replace(' ','')+".C")

mht_bins = np.concatenate(([20], np.arange(100, 400, 20), [400, 800, 1600, 3000]))
met_bins = np.concatenate(([20], np.arange(100, 300, 10), np.arange(300, 800, 100), [800, 1600]))
R_bins = np.arange(0.5, 2.5, step=0.1)

#urnon_1D("R_NoMu_central", "H_{T}^{miss}/p_{T}^{miss} (no mu)", R_bins, "HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight || HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight || HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight ||HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", "mht_NoMucentral > 300 ", "H_{T}^{miss} > 300 GeV", False)

#turnon_1D("mht_NoMu_central", "H_{T}^{miss} (no mu) / GeV", mht_bins, "HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight || HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight || HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight ||HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", "R_NoMu_central < 1.25", "H_{T}^{miss}/p_{T}^{miss} < 1.25", True)

turnon_1D("mht_NoMu_central", "H_{T}^{miss} (no mu) / GeV", mht_bins, "HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight || HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight || HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight ||HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight || HLT_PFHT900", "R_NoMu_central < 1.25", "H_{T}^{miss}/p_{T}^{miss} < 1.25", True)
