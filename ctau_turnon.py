import ROOT
import os

import sys
sys.path.insert(0, '/home/hep/vc1117/tdr-style')
import tdrstyle

path = "friends/"

samples = os.listdir(path)

chain = ROOT.TChain("Friends")

for sample in samples:
    if "SMS" in sample and "1000" in sample:
        print sample
        friends = os.listdir(os.path.join(path, sample))
        for friend in friends:
            chain.AddFile(os.path.join(path, sample, friend))

class Variable():
    def __init__(self, name, varexp, title, xmin, xmax):
        self.name = name
        self.varexp = varexp
        self.title = title
        self.xmin = xmin
        self.xmax = xmax

variables = [
        Variable("jet1_pt_central", "jet1_pt_central", "leading jet p_{T} / GeV", 0, 1000),
        Variable("jet2_pt_central", "jet2_pt_central", "sub-leading jet p_{T} / GeV", 0, 1000),
        Variable("jet1_jetId_central", "jet1_jetId_central", "leading jet ID", 0, 4),
        Variable("jet2_jetId_central", "jet2_jetId_central", "sub-leading jet ID", 0, 4),
        Variable("jet1_CHM_central", "jet1_CHM_central", "leading jet CHM", 0, 80),
        Variable("jet1_nConstituents_central", "jet1_nConstituents_central", "leading jet number of constituents", 0, 80),
        Variable("jet1_chEmEF_central", "jet1_chEmEF_central", "leading jet chEmEF", 0, 0.05),
        Variable("jet1_chHEF_central", "jet1_chHEF_central", "leading jet chHEF", 0, 1),
        Variable("MET_NoMu", "MET_NoMu", "E_{}^{miss} / GeV", 200, 1000),
        Variable("ht_central", "ht_central", "H_{T} / GeV", 100, 1000),
        ]

for variable in variables:
    c = ROOT.TCanvas()
    baseline = "(mht_NoMu_central > 300. && R_NoMu_central < 1.25 && jet1_fromLLP_central == 0)"
    uncompressed = "1"
    uncompressed = "llpinfo_lsp_mass == 200.0 && llpinfo_llp_mass == 2400."
    compressed = "llpinfo_llp_mass == 2000. && llpinfo_lsp_mass == 1800."
    #passed = "(HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight || HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight || HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight || HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight)"
    passed = "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"
    #passed = "1"

    hist_passed = ROOT.TH1F("hist_passed"+variable.name, "hist_passed", 25, variable.xmin, variable.xmax)
    hist_notpassed = ROOT.TH1F("hist_notpassed"+variable.name, "hist_notpassed", 25, variable.xmin, variable.xmax)
    hist_compressed = ROOT.TH1F("hist_compressed"+variable.name, "hist_compressed", 25, variable.xmin, variable.xmax)
    print "projecting", variable.name, "with limits", variable.xmin, variable.xmax

    chain.Project("hist_passed"+variable.name, variable.varexp,'{0} * {1} * {2}'.format(baseline, uncompressed, passed))
    chain.Project("hist_notpassed"+variable.name, variable.varexp,'{0} * {1} * !{2}'.format(baseline, uncompressed, passed))
    chain.Project("hist_compressed"+variable.name, variable.varexp,'{0} * {1} * {2}'.format(baseline, compressed, passed))

    c.Draw("")
    
    print hist_passed.Integral(), hist_notpassed.Integral()

    hist_passed.Scale(1./hist_passed.Integral())
    hist_notpassed.Scale(1./hist_notpassed.Integral())
    hist_compressed.Scale(1./hist_compressed.Integral())

    hist_passed.SetLineColor(ROOT.kSpring)
    hist_notpassed.SetLineColor(ROOT.kRed)

    hist_passed.GetXaxis().SetTitle(variable.title)
    hist_passed.GetYaxis().SetTitle("number of entries, normalised")
    hist_passed.GetYaxis().SetRangeUser(0, 1.1*max(hist_passed.GetMaximum(), hist_notpassed.GetMaximum(), hist_compressed.GetMaximum()))

    hist_passed.Draw("HIST")
    hist_compressed.Draw("HISTSAME")
    hist_notpassed.Draw("HISTSAME")

    leg = ROOT.TLegend(0.75, 0.6, 0.95, 0.8)
    leg.AddEntry(hist_passed, "passed", "l")
    leg.AddEntry(hist_notpassed, "not passed", "l")
    leg.AddEntry(compressed, "compressed", "l")
    leg.Draw("SAME")
 
    p = ROOT.TLatex()
    p.SetTextSize(0.03)
    p.DrawLatexNDC(0.2, 0.95, "preselection applied")
    p.DrawLatexNDC(0.8, 0.4, "H_{T}^{miss} > 300 GeV")
    p.DrawLatexNDC(0.8, 0.5, "H_{T}^{miss}/p_{T}^{miss} < 1.25")
    p.DrawLatexNDC(0.5, 0.95, "uncompr., c#tau = 1m, 10m and 100m")

    c.SaveAs(variable.name+".pdf")

