import ROOT
import argparse


def main(args):
    # base plotting code take from: https://root.cern.ch/root/html/tutorials/hist/ratioplot.C.html

    osc_file = ROOT.TFile(args.osc_fname)
    osc_histo = osc_file.Get(args.histo_name)

    baseline_file = ROOT.TFile(args.baseline_fname)
    baseline_histo = baseline_file.Get(args.histo_name)

    canvas = ROOT.TCanvas("c", "canvas", 800, 800, 800, 800)
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0) # Upper and lower plot are joined
    pad1.SetGridx()
    pad1.Draw()
    pad1.cd()
    
    baseline_histo.Draw("")
    osc_histo.Draw("HIST SAME")

    # Do not draw the Y axis label on the upper plot and redraw a small
    # axis instead, in order to avoid the first label (0) to be clipped.
    osc_histo.GetYaxis().SetLabelSize(0)
    axis = ROOT.TGaxis( -5, 20, -5, 220, 20,220,510,"")
    axis.SetLabelFont(43)
    axis.SetLabelSize(15)
    axis.Draw()

    # make legend
    legend = ROOT.TLegend(0.3,0.7,0.7,0.9)
    legend.AddEntry(osc_histo, "Oscillation \\nu_{\mu}")
    legend.AddEntry(baseline_histo, "Baseline \\nu_{\mu}")
    legend.Draw()

    canvas.cd()

    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.2)
    pad2.SetGridx()
    pad2.Draw()
    pad2.cd()

    ratio_histo = osc_histo.Clone("ratio")
    ratio_histo.SetLineColor(ROOT.kBlack)
    ratio_histo.SetMinimum(0.8)
    ratio_histo.SetMaximum(1.35)
    ratio_histo.Sumw2()
    ratio_histo.SetStats(0)
    ratio_histo.Divide(baseline_histo)
    ratio_histo.SetMarkerStyle(21)
    ratio_histo.Draw("HIST")

    # osc_histo settings
    osc_histo.SetLineColor(ROOT.kBlue+1)
    osc_histo.SetLineWidth(2)
    
    # Y axis baseline_histo plot settings
    baseline_histo.GetYaxis().SetTitle("# Events")
    baseline_histo.GetYaxis().SetTitleSize(20)
    baseline_histo.GetYaxis().SetTitleFont(43)
    baseline_histo.GetYaxis().SetTitleOffset(1.55)

    baseline_histo.SetTitleSize(60)
    baseline_histo.SetTitleFont(43)
    baseline_histo.SetTitle("Ratio of Oscillation/Baseline \\nu_{\mu} Events")
    
    # baseline_histo settings
    baseline_histo.SetLineColor(ROOT.kRed)
    baseline_histo.SetLineWidth(2)
    
    # Ratio plot (ratio_histo) settings
    ratio_histo.SetTitle("") # Remove the ratio title
    
    # Y axis ratio plot settings
    ratio_histo.GetYaxis().SetTitle("ratio osc/baseline ")
    ratio_histo.GetYaxis().SetNdivisions(505)
    ratio_histo.GetYaxis().SetTitleSize(20)
    ratio_histo.GetYaxis().SetTitleFont(43)
    ratio_histo.GetYaxis().SetTitleOffset(1.55)
    ratio_histo.GetYaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
    ratio_histo.GetYaxis().SetLabelSize(15)
    
    # X axis ratio plot settings
    ratio_histo.GetXaxis().SetTitle("Truth Energy [GeV]")
    ratio_histo.GetXaxis().SetTitleSize(20)
    ratio_histo.GetXaxis().SetTitleFont(43)
    ratio_histo.GetXaxis().SetTitleOffset(4.)
    ratio_histo.GetXaxis().SetLabelFont(43) # Absolute font size in pixel (precision 3)
    ratio_histo.GetXaxis().SetLabelSize(15)

    canvas.Update()
    canvas.Show()
    raw_input("Press enter to continue...")
    canvas.SaveAs("ratio.pdf")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--baseline_fname", required=True)
    parser.add_argument("-o", "--osc_fname", required=True)
    parser.add_argument("-n", "--histo_name", default="numu_trueE")
    args = parser.parse_args()
    main(args)
