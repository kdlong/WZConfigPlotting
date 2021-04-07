#!/usr/bin/env python
import Utilities.plot_functions as plotter
import Utilities.helper_functions as helper
import argparse
import ROOT
import Utilities.config_object as config_object
import Utilities.UserInput as UserInput
import os
from Utilities.ConfigHistFactory import ConfigHistFactory 
from Utilities.prettytable import PrettyTable
import math
import sys
import array
import datetime
from Utilities.scripts import makeSimpleHtml
from IPython import embed
import logging

def getComLineArgs():
    parser = UserInput.getDefaultParser()
    parser.add_argument("-s", "--selection", type=str, required=True,
                        help="Specificy selection level to run over")
    parser.add_argument("-b", "--branches", type=str, default="all",
                        help="List (separate by commas) of names of branches "
                        "in root and config file to plot") 
    parser.add_argument("--systematics", type=str, required=True,
                        help="Specificy systematics to plot")
    return parser.parse_args()

def rebinHist(hist, rebin):
    if rebin:
        hist_name = hist.GetName()
        if len(rebin) == 1:
            hist.Rebin(int(rebin[0]))
        else:
            bins = array.array('d', rebin)
            rebin_hist = hist.Rebin(len(bins)-1, "", bins)
            #rebin_hist = hist.Rebin(len(bins)-1, "temp", bins)

        hist.SetName(hist_name)
    return hist
def main():
    args = getComLineArgs()
    ROOT.gROOT.SetBatch(True)

    branches = [x.strip() for x in args.branches.split(",")]
    
    (plot_path, html_path) = helper.getPlotPaths(args.selection, args.folder_name, True)

    meta_info = '-'*80 + '\n' + \
        'Script called at %s\n' % datetime.datetime.now() + \
        'The command was: %s\n' % ' '.join(sys.argv) + \
        '-'*80 + '\n'

    rtfile = ROOT.TFile(args.hist_file)
    #colors = [ROOT.TColor.GetColor(x) for x in ["#a6cee3","#1f78b4","#b2df8a", "#33a02c", "#fb9a99", "#fdbf6f",
    colors = [ROOT.TColor.GetColor(x) for x in ["#b2df8a", "#33a02c", "#fb9a99", "#fdbf6f",
                    "#e31a1c","#ff7f00","#cab2d6","#6a3d9a", "#ffff99", "#b15928",]
    ]

    for branch in args.branches.split(","):
        with open("temp.txt", "w") as mc_file:
            mc_file.write(meta_info)
            mc_file.write("Selection: %s" % args.selection)
            mc_file.write("\nPlotting branch: %s\n" % branch)

        systematics = args.systematics.split(",")
        systs = "_".join(systematics)
        plot_name = "_".join([branch, systs]) + args.append_to_name 
        hist_stack = ROOT.THStack("stack_"+branch, "stack_"+branch)
        file_names = args.files_to_plot.split(",")
        if plot_name.count(file_names[0]) > 1:
            names = plot_name.split("_"+file_names[0])
            plot_name = "".join([names[0], "_"+file_names[0]]+names[1:])

        # To avoid having to paste the file name again and again if you want
        # to plot multiple variations
        while len(file_names) < len(systematics):
            file_names.append(file_names[0]) 
            
        for i, file_name in enumerate(file_names):
            systematic = systematics[i]
            central_hist = 0
            for chan in args.channels.split(","):
                hist_name = file_name.replace("_standalone", "") + "/"+ "_".join([branch, chan])
                print hist_name
                names = hist_name.rsplit("_"+chan, 1)
                uphist_name = "_".join([names[0], systematic+"Up", chan]) if systematic != "None" else hist_name
                downhist_name = "_".join([names[0],systematic+"Down", chan]) if systematic != "None" else hist_name

                if not central_hist:
                    central_hist = rtfile.Get(hist_name)
                    central_hist.SetName(hist_name+"_%i" % i)
                    if not central_hist:
                        raise ValueError("Failed to find hist %s in file %s" % (hist_name, file_name))
                    up_hist = rtfile.Get(uphist_name)
                    print uphist_name
                    up_hist.SetName(uphist_name+"_%i" % i)
                    down_hist = rtfile.Get(downhist_name)
                    down_hist.SetLineColor(ROOT.kYellow)
                    down_hist.SetName(downhist_name+"_%i" % i)
                    if not up_hist or not down_hist:
                        raise ValueError("Failed to find hist %s for variation %s in file %s" % (uphist_name, systematic, file_name))
                    central_hist = rebinHist(central_hist, args.rebin)
                    up_hist = rebinHist(up_hist, args.rebin)
                    down_hist = rebinHist(down_hist, args.rebin)
                    central_hist.Draw()
                else:
                    central_hist_chan = rtfile.Get(hist_name)
                    central_hist_chan.SetName(hist_name+"_%i" % i)
                    central_hist_chan = rebinHist(central_hist_chan, args.rebin)
                    central_hist.Add(central_hist_chan)
                    up_hist_chan = rtfile.Get(uphist_name)
                    down_hist_chan = rtfile.Get(downhist_name)
                    up_hist_chan = rebinHist(up_hist_chan, args.rebin)
                    down_hist_chan = rebinHist(down_hist_chan, args.rebin)
                    up_hist.Add(up_hist_chan)
                    down_hist.Add(down_hist_chan)
            path = "/cms/kdlong" if "hep.wisc.edu" in os.environ['HOSTNAME'] else \
                "/afs/cern.ch/user/k/kelong/work"
            config_factory = ConfigHistFactory(
                "%s/AnalysisDatasetManager" % path,
                args.selection.split("_")[0],
            )
            scale_fac = 1
            scale = True
            # Normalize to data in CR
            if file_name in ["wz", "wz_standalone"]:
                scale_fac = 0.8224 
            elif file_name in ["wz-powheg", "wz-powheg_standalone"]:
                scale_fac =(191-(28.78+3.84+18.35+22.27))/139.43
            if scale:
                central_hist.Scale(scale_fac)
                up_hist.Scale(scale_fac)
                down_hist.Scale(scale_fac)
            with open("temp.txt", "a") as mc_file:
                mc_file.write("\nYield for %s is %0.2f (%i entries)" % (file_name, central_hist.Integral(), central_hist.GetEntries()))
                mc_file.write("\nYield up for var %s for %s is %0.2f (%i entries)" % (file_name, systematic, up_hist.Integral(), up_hist.GetEntries()))
                mc_file.write("\nYield down for var %s for %s is %0.2f (%i entries)" % (file_name, systematic, down_hist.Integral(), down_hist.GetEntries()))
            
            branch_name = branch.replace("_Fakes", "")
            config_factory.setHistAttributes(central_hist, branch_name, file_name)
            config_factory.setHistAttributes(up_hist, branch_name, file_name)
            config_factory.setHistAttributes(down_hist, branch_name, file_name)
            
            central_hist.SetMinimum(0.001)
            central_hist.SetFillColor(0)
            up_hist.SetFillColor(0)
            down_hist.SetFillColor(0)
            
            central_hist.SetLineStyle(1)
            central_hist.SetLineWidth(1)
            up_hist.SetLineWidth(1)
            down_hist.SetLineWidth(1)
            up_hist.SetLineStyle(5)
            down_hist.SetLineStyle(3)
            if i != 0 and file_names[i] == file_names[i-1]:
                up_hist.SetLineColor(colors[i])
                down_hist.SetLineColor(colors[i])

            hist_stack.Add(central_hist)
            hist_stack.Add(up_hist)
            hist_stack.Add(down_hist)

        canvas_dimensions = [800, 800] if "unrolled" not in branch else [1200, 800]
        canvas = ROOT.TCanvas("canvas", "canvas", *canvas_dimensions)
        hist_stack.Draw("nostack hist")
        hist_stack.SetMinimum(central_hist.GetMinimum()*args.scaleymin)
        hist_stack.SetMaximum(central_hist.GetMaximum()*args.scaleymax)
        hist_stack.GetHistogram().GetYaxis().SetTitle("Events / bin")
        hist_stack.GetHistogram().GetYaxis().SetTitleOffset(1.05)

        text_box = ROOT.TPaveText(0.65, 0.75-(0.1*len(systematics))*args.scalelegy, 0.9, 0.85, "NDCnb")
        text_box.SetFillColor(0)
        text_box.SetTextFont(42)
        text_box.AddText("Systematic variation")
        for s in systematics:
            text_box.AddText(s)
        text_box.AddText("; ".join(set(file_names)))
        for line, h in zip(text_box.GetListOfLines()[1:], hist_stack.GetHists()[1::3]):
            line.SetTextColor(h.GetLineColor())
        text_box.Draw()
        if args.logy:
            canvas.SetLogy()

        if not args.no_ratio:
            canvas = plotter.splitCanvas(canvas, canvas_dimensions,
                    "syst./cent.",
                    [float(i) for i in args.ratio_range],
                    args.uncertainties == "none"
            )
            ratioPad = canvas.FindObject("ratioPad")
            hist = ratioPad.GetListOfPrimitives().FindObject("canvas_central_ratioHist")
            hist.GetXaxis().SetLabelSize(0.15)
            hist.GetXaxis().SetLabelOffset(0.03)
            hist.GetXaxis().SetTitleOffset(1.15)
            canvas.Update()
        helper.savePlot(canvas, plot_path, html_path, plot_name, True, args)
        makeSimpleHtml.writeHTML(html_path.replace("/plots",""), args.selection)

if __name__ == "__main__":
    main()

