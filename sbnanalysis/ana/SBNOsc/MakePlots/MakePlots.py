## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Plotting Outputs
#
#  (from covariance and sensitivity
#  calculations)
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import ROOT
from ROOT import TFile, TCanvas, TH2D, TH1D, TGraph, TGraph2D, TStyle, TLegend, THStack, TPad
import argparse
import os


## Main function
## ~~~~~~~~~~~~~

def main(args):
    
    
    # Count plots
    
#    countfile = TFile(args.cntfile)
#    
#    numu = countfile.Get('numu_counts')
#    numu_nkg = countfile.Get('numu_bkgs')
#    
#    numu_canvas = TCanvas("numu_canvas", "", 1600, 550)
#    numu_canvas.Divide(, 1);

    ### Later might come in handy: to loop through objects in TFile: 
    # [key.GetName() for key in file.GetListOfKeys()]
    
    
    # Covariance, fractional covariance and correlation

    covfile = TFile(args.covfile)

    covcanvas = TCanvas()

    gStyle = TStyle()
    gStyle.SetPadLeftMargin(0.25); gStyle.SetPadRightMargin(0.15)
    gStyle.SetPalette(56)

    for matname in ('cov', 'fcov', 'corr'):

        mat = covfile.Get(matname)

        mat.GetYaxis().LabelsOption('v')    # doesn't work...
        mat.GetYaxis().SetLabelSize(0.07)
        mat.GetXaxis().SetLabelSize(0.07)

        mat.SetTitleSize(0.3, 't')  # doesn't work...

        if matname == 'corr': 
            mat.GetZaxis().SetRangeUser(-0.4, 1)
            mat.SetTitle("Flux Correlation Matrix")

        mat.Draw("colz")
        mat.SetStats(False)
        covcanvas.SetLeftMargin(0.25)
        covcanvas.Update()
        covcanvas.SaveAs(args.outdir + matname + "_plot.pdf")


    # Chi squareds
    
    chi2file = TFile(args.chifile)
    
    chi2 = chi2file.Get('chisq')
    
    chi2canvas = TCanvas()
    
    chi2.SetTitle('#chi^{2}; log_{10}(sin^{2}(2#theta)); log_{10}(#Delta m^{2}); #chi^{2}');
    gStyle.SetPalette(1)
    chi2.Draw('surf1')
    chi2canvas.SaveAs(args.outdir + "chisq.pdf")
    
    
    # Contours
    contcanvas = TCanvas('cont_canvas', '', 1020, 990)
    gStyle.SetPadLeftMargin(0.15); gStyle.SetPadRightMargin(0.15)
    
    colours = [30, 38, 46]
    contours = [chi2file.Get('90pct'), 
                chi2file.Get('3sigma'), 
                chi2file.Get('5sigma')]
    
    for g in range(len(contours)):
        
        contours[g].SetMarkerStyle(20)
        contours[g].SetMarkerSize(0.25)
        contours[g].SetMarkerColor(colours[g])
        contours[g].SetLineColor(colours[g])
    
    gr_range = TGraph()
    gr_range.SetPoint(0, 0.001, 0.01)
    gr_range.SetPoint(1, 1, 100)
    gr_range.SetMarkerColor(0)
    
    bestfit = TGraph()
    bestfit.SetPoint(0, 0.062, 1.7)
    bestfit.SetMarkerStyle(29)
    bestfit.SetMarkerSize(1.6)
    bestfit.SetMarkerColor(40)
    
    gr_range.SetTitle('SBN Sensitivity; sin^{2}(2#theta); #Delta m^{2} (eV^{2})')
    
    legend = TLegend()
    legend.AddEntry(contours[0], '90% CL', 'l')
    legend.AddEntry(contours[1], '3#sigma CL', 'l')
    legend.AddEntry(contours[2], '5#sigma CL', 'l')
    legend.AddEntry(bestfit, 'Best Fit Point', 'p')
    
    contcanvas.SetLogy()
    contcanvas.SetLogx()
    
    gr_range.Draw('AP')
    gr_range.GetXaxis().SetRangeUser(0.001, 1)
    gr_range.GetYaxis().SetRangeUser(0.01, 100)
    
    contours[0].Draw('P same')
    contours[1].Draw('P same')
    contours[2].Draw('P same')
    
    legend.Draw()
    bestfit.Draw('P same')
    
    contcanvas.SetLeftMargin(0.15)
    contcanvas.Update()
    contcanvas.SaveAs(args.outdir+'Sensitivity.pdf')
    
    

def compare_w_proposal(args):
    
    chi2file = TFile(args.chifile)
    
    gStyle = TStyle()
    gStyle.SetPadLeftMargin(0.15); gStyle.SetPadRightMargin(0.15)
    
    colours = [30, 38, 46]
    contours = [chi2file.Get('90pct'), 
                chi2file.Get('3sigma'), 
                chi2file.Get('5sigma')]
    
    propcontours = []
    contournames = ['90pct', '3s', '5s']
    contourtitles = ['90% Confidence Level', '3$\\sigma$ Confidence Level', '5$\\sigma$ Confidence Level']
    
    for i in range(len(contours)):

        with open('numu_'+contournames[i]+'.txt') as f:
            for line in f:
                x.append(line.split(', ')[0])
                y.append(line.split(', ')[1])
        propcontours.append(TGraph2D())
        for j in range(len(x)):
            propcontours[i].SetPoint(j, x[j], y[j])

        tempcanvas = TCanvas('temp_canvas', '', 1020, 990)

        templegend = TLegend()
        legend.AddEntry(contours[i], 'Our contour', 'l')
        legend.AddEntry(propcontours[1], 'From proposal', 'l')
        legend.AddEntry(bestfit, 'Best Fit Point', 'p')

        tempcanvas.SetLogy()
        tempcanvas.SetLogx()

        gr_range.SetTitle(contourtitles[i]+' Comparison; sin^{2}(2#theta); #Delta m^{2} (eV^{2})')

        gr_range.Draw('AP')
        gr_range.GetXaxis().SetRangeUser(0.001, 1)
        gr_range.GetYaxis().SetRangeUser(0.01, 100)
        
        for lst in (contours, propcontours):
            lst[i].SetMarkerStyle(20)
            lst[i].SetMarkerSize(0.25)
            lst[i].SetMarkerColor(colours[i])
            lst[i].SetLineColor(colours[i])
        
        contours[i].Draw('P same')
        propcontours[i].Draw('P same')

        templegend.Draw()
        bestfit.Draw('P same')

        tempcanvas.SaveAs(args.outdir+contournames+'_comparison.pdf')
    
    

if __name__ == "__main__":
    
    buildpath = os.environ['SBN_LIB_DIR']
    if not buildpath:
        print('ERROR: SBNDDAQ_ANALYSIS_BUILD_PATH not set')
        sys.exit()
    
    ROOT.gROOT.ProcessLine('.L ' + buildpath + '/libsbnanalysis_Event.so')
    ROOT.gROOT.ProcessLine('.L ' + buildpath + '/libsbnanalysis_SBNOsc_classes.so')
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-chi", "--chifile", required = True)
    parser.add_argument("-cov", "--covfile", required = True)
#    parser.add_argument("-cts", "--cntfile", required = True)
    parser.add_argument("-o", "--outdir", required = True)
    parser.add_argument("-comp", "--compare", default = False)

    main(parser.parse_args())
    
    if parser.parse_args().compare: 
        compare_w_proposal(parser_parse_args())
        with open('filename') as f:
            for line in f:
                data = [float(x) for x in line.split(",")]





