import ROOT
import numpy as np
import os
import multiprocessing

import sys
import tdrstyle

from array import array

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

    filePath = "/vols/cms/vc1117/AN-18-199/nanoAOD_friends/trigger"
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

path_signal = "(HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight || HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight || HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight ||HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight || HLT_PFHT900)"

x_bins = array('d', [200, 210, 220, 230, 240, 250, 260, 270, 280, 300, 350, 400, 500, 600])
y_bins = array('d', [0, 1.25, 1.5, 1.7,  2., 2.5])

procs = [
        ["SingleElectron", "SingleEle*.root", "e", "(1)"],
        ["SingleMuon", "SingleMu*.root", "mu", "(1)"],
        ["t#bar{t}: electron", "TT*.root", "e_mc", "(electron_trigger && n_veto_muons == 0 && n_tight_electrons > 0)"],
        ["t#bar{t}: muon", "TT*.root", "mu_mc", "(muon_trigger && n_veto_electrons == 0 && n_loose_muons > 0)"],
        ]


f = ROOT.TFile.Open("trigger_SF.root", "RECREATE")
for proc in procs:

    t = ROOT.TChain("Friends")
    title = proc[0]
    t.Add(proc[1])
    label = proc[2]
    preselection = proc[3]

    eff_nominal = ROOT.TH2D("eff_nominal", "eff_nominal", len(x_bins)-1, x_bins, len(y_bins)-1, y_bins)
    eff_up = ROOT.TH2D("eff_up", "eff_up", len(x_bins)-1, x_bins, len(y_bins)-1, y_bins)
    eff_down = ROOT.TH2D("eff_down", "eff_down", len(x_bins)-1, x_bins, len(y_bins)-1,y_bins)

    before = ROOT.TH2D("before_"+label, "before_"+label, len(x_bins)-1, x_bins, len(y_bins)-1,y_bins)
    after = ROOT.TH2D("after_"+label, "after_"+label, len(x_bins)-1, x_bins, len(y_bins)-1,y_bins)

    t.Project("before_"+label, "R_NoMu_central:mht_NoMu_central", preselection)
    t.Project("after_"+label, "R_NoMu_central:mht_NoMu_central", preselection+"*"+path_signal)

    before.Sumw2()
    after.Sumw2()

    eff = after.Clone("ratio")
    eff = ROOT.TEfficiency(after, before)

    eff_nominal.GetXaxis().SetTitle("H_{T}^{miss} (no mu) / GeV")
    eff_nominal.GetYaxis().SetTitle("H_{T}^{miss}/ p_{T}^{miss} (no mu)")
    eff_nominal.SetTitle("eff_nominal")

    eff_up.GetXaxis().SetTitle("H_{T}^{miss} (no mu) / GeV")
    eff_up.GetYaxis().SetTitle("H_{T}^{miss}/ p_{T}^{miss} (no mu)")
    eff_up.SetTitle("eff_up")

    eff_down.GetXaxis().SetTitle("H_{T}^{miss} (no mu) / GeV")
    eff_down.GetYaxis().SetTitle("H_{T}^{miss}/ p_{T}^{miss} (no mu)")
    eff_down.SetTitle("eff_down")

    for j in range(len(y_bins)+1):
        for i in range(len(x_bins)+1):
            counter = j*(len(x_bins)+1) + i
            eff_nominal.SetBinContent(i, j, eff.GetEfficiency(counter))
            eff_up.SetBinContent(i, j, eff.GetEfficiency(counter)+eff.GetEfficiencyErrorUp(counter))
            eff_down.SetBinContent(i,j , eff.GetEfficiency(counter)-eff.GetEfficiencyErrorLow(counter))


    f.mkdir(label)
    f.cd(label)

    eff_nominal.Write()
    eff_up.Write()
    eff_down.Write()

    pTitle=ROOT.TPaveText(0.1,0.96,0.6,0.96,"NDC")
    pTitle.SetFillColor(ROOT.kWhite)
    pTitle.SetBorderSize(0)
    pTitle.SetTextFont(43)
    pTitle.SetTextSize(23)
    pTitle.SetTextAlign(31)
    pTitle.AddText(title)

    c = ROOT.TCanvas(label, label)
    c.cd()
    ROOT.gPad.SetRightMargin(0.2)
    ROOT.gPad.SetTopMargin(0.06)

    c.Draw()
    eff_nominal.Draw("COLZtext45")
    tdrstyle.cmsPrel(35900,  energy=13,  simOnly=False,  onLeft=True,  sp=0, textScale=1.3, xoffset=0., thisIsPrelim=True)
    pTitle.Draw("Same")
    c.SaveAs("eff_nominal_"+label+".pdf")

    eff_up.Draw("COLZtext45")
    tdrstyle.cmsPrel(35900,  energy=13,  simOnly=False,  onLeft=True,  sp=0, textScale=1.3, xoffset=0., thisIsPrelim=True)
    pTitle.Draw("Same")
    c.SaveAs("eff_up_"+label+".pdf")

    eff_down.Draw("COLZtext45")
    tdrstyle.cmsPrel(35900,  energy=13,  simOnly=False,  onLeft=True,  sp=0, textScale=1.3, xoffset=0., thisIsPrelim=True)
    pTitle.Draw("Same")
    c.SaveAs("eff_down_"+label+".pdf")

    before_proj = ROOT.TH1F("before_proj_"+label, "before_proj_"+label, len(x_bins)-1, x_bins)
    after_proj = ROOT.TH1F("after_proj_"+label, "after_proj_"+label, len(x_bins)-1, x_bins)

    t.Project("before_proj_"+label, "mht_NoMu_central", preselection+"*(R_NoMu_central < 1.25)")
    t.Project("after_proj_"+label, "mht_NoMu_central", preselection+"*(R_NoMu_central < 1.25)"+"*"+path_signal)

    before_proj.Sumw2()
    after_proj.Sumw2()

    after_proj.GetXaxis().SetTitle("H_{T}^{miss} (no mu) / GeV")
    after_proj.GetYaxis().SetTitle("efficiency")

    eff_proj = after_proj.Clone("eff_proj")
    eff_proj.Divide(after_proj, before_proj, 1, 1, "b")
    eff_proj.SetDirectory(0)
    proc.append(eff_proj)

f.Close()

l = ROOT.TLegend(0.60, 0.45, 0.9, 0.75)
l.SetTextSize(0.05)
canvas_proj = ROOT.TCanvas("canvas_proj", "canvas_proj")

for proc in procs:
    hist = proc[4]
    title = proc[0]
    hist.GetYaxis().SetRangeUser(0.85, 1.03)
    if title == "SingleElectron":
        hist.Draw("PLC PMC")
    else:
        hist.Draw("SAME PLC PMC")
    l.AddEntry(hist, title, "p")

p = ROOT.TLatex()
p.SetTextSize(0.04)
p.DrawLatexNDC(0.4, 0.85, "H_{T}^{miss}/p_{T}^{miss} < 1.25")

line = ROOT.TGraph(2)
line.SetPoint(0, 0, 1)
line.SetPoint(1, x_bins[-1], 1) 
line.SetLineColor(ROOT.kRed)
line.SetLineStyle(9)
line.Draw("L")
line2 = ROOT.TLine(300, 0.85, 300, 1.03)
line2.SetLineStyle(2)
line2.Draw("SAME")

l.Draw("SAME")
tdrstyle.cmsPrel(35900, energy=13, simOnly=False, onLeft=True, sp=0, textScale=1, xoffset=0, thisIsPrelim=True)
canvas_proj.SaveAs("proj.pdf")
canvas_proj.Print("proj_printed.pdf")

'''

if do_fit:

    eff_fit = ROOT.TF1("eff_fit", "[4] + ([0] - [4])/(1+[3]*TMath::Exp(-[1]*(x-[2])))^([5])", min(x_bins), max(x_bins))
    eff_prefit = ROOT.TF1("eff_prefit", "[0]/(1 + TMath::Exp(-[1]*(x-[2])))", min(x_bins), max(x_bins))
    eff_prefit.SetParameters(1, 0.2, 200)
    eff_fit.SetLineColor(ROOT.kOrange)
    eff.Fit("eff_prefit", "")
    eff_fit.SetParameters(eff_prefit.GetParameter(0), eff_prefit.GetParameter(1), eff_prefit.GetParameter(2), 1, 0, 1)
    eff.Fit("eff_fit", "")
'''


