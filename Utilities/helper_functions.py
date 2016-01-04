import ROOT
import plot_functions as plotter
import Utilities.WeightInfo as WeightInfo
import Utilities.WeightedHistProducer as WeightedHistProducer
from Utilities.ConfigHistFactory import ConfigHistFactory 
from collections import OrderedDict
import os
import glob
import logging
import datetime
import shutil

def makePlot(hist_stack, data_hist, branch_name, args):
    canvas = ROOT.TCanvas("%s_canvas" % branch_name, branch_name, 800, 600) 
    hists = hist_stack.GetHists()
    hist_stack.Draw("nostack hist" if args.nostack else "hist")
    if data_hist:
        data_hist.Draw("e1 same")
    if not args.no_decorations:
        ROOT.dotrootImport('nsmith-/CMSPlotDecorations')
        ROOT.CMSlumi(canvas, 0, 11, "%0.2f fb^{-1} (13 TeV)" % (args.luminosity/1000.))
    hist_stack.GetYaxis().SetTitleSize(hists[0].GetYaxis().GetTitleSize())    
    hist_stack.GetYaxis().SetTitleOffset(hists[0].GetYaxis().GetTitleOffset())    
    hist_stack.GetYaxis().SetTitle(
        hists[0].GetYaxis().GetTitle())
    hist_stack.GetHistogram().GetXaxis().SetTitle(
        hists[0].GetXaxis().GetTitle())
    hist_stack.GetHistogram().SetLabelSize(0.04)
    if hist_stack.GetMaximum() < hists[0].GetMaximum():
        hist_stack.SetMaximum(hists[0].GetMaximum())
        hist_stack.SetMinimum(hists[0].GetMinimum())
    else:
        new_max = 1.1*hist_stack.GetMaximum() if not data_hist else \
                1.1*max(data_hist.GetMaximum(), hist_stack.GetMaximum()) 
        hist_stack.SetMaximum(new_max)
        hist_stack.SetMinimum(0.001)
    if not args.no_errors:
        histErrors = getHistErrors(hist_stack, args.nostack)
        for error_hist in histErrors:
            error_hist.Draw("E2 same")
    legend = getPrettyLegend(hist_stack, data_hist, args.legend_left)
    legend.Draw()
    if not args.no_ratio:
        canvas = plotter.splitCanvas(canvas, "stack", ("#Sigma MC", "Data"))
    return canvas
def getHistErrors(hist_stack, separate):
    histErrors = []
    for hist in hist_stack.GetHists():
        error_hist = plotter.getHistErrors(hist)
        if separate:
            error_hist.SetFillColor(hist.GetFillColor())
            histErrors.append(error_hist)
        else:
            error_hist.SetFillColor(ROOT.kBlack)
            if len(histErrors) == 0:
                histErrors.append(error_hist)
            else:
                histErrors[0].Add(error_hist)
    return histErrors
def getPrettyLegend(hist_stack, data_hist, left):
    hists = hist_stack.GetHists()
    xcoords = [.15, .35] if left else [.78, .93]
    unique_entries = len(set([x.GetFillColor() for x in hists]))
    ycoords = [.9, .9 - 0.06*unique_entries]
    legend = ROOT.TLegend(xcoords[0], ycoords[0], xcoords[1], ycoords[1])
    legend.SetFillColor(0)
    if data_hist:
        legend.AddEntry(data_hist, data_hist.GetTitle(), "lp")
    hist_names = []
    for hist in reversed(hists):
        if hist.GetTitle() not in hist_names:
            legend.AddEntry(hist, hist.GetTitle(), "f")
        hist_names.append(hist.GetTitle())
    return legend
def buildChain(filelist, treename):
    chain = ROOT.TChain(treename)
    for filename in glob.glob(filelist):
        filename = filename.strip()
        if ".root" not in filename or not os.path.isfile(filename):
            raise IOError("%s is not a valid root file!" % filename)
        chain.Add(filename)
    return chain
def getHistFactory(config_factory, selection, filelist, luminosity=1, cut_string=""):
    mc_info = config_factory.getMonteCarloInfo()
    all_files = config_factory.getFileInfo()
    hist_factory = OrderedDict() 
    for name in filelist:
        if name not in all_files.keys():
            logging.warning("%s is not a valid file name (must match a definition in FileInfo/%s.json)" % \
                (name, selection))
            continue
        hist_factory[name] = dict(all_files[name])
        if "data" not in name:
            metaTree = buildChain(hist_factory[name]["file_path"],
                    "eee/metaInfo")
            kfac = 1. if 'kfactor' not in mc_info[name].keys() else mc_info[name]['kfactor']
            weight_info = WeightInfo.WeightInfoProducer(metaTree, 
                    mc_info[name]['cross_section']*kfac,
                    "summedWeights").produce()
            histProducer = WeightedHistProducer.WeightedHistProducer(weight_info, "GenWeight")  
            histProducer.setLumi(luminosity)
            histProducer.setCutString(cut_string)
        else:
            histProducer = WeightedHistProducer.WeightedHistProducer(
                    WeightInfo.WeightInfo(1, 1,), "")  
        hist_factory[name].update({"histProducer" : histProducer})
    return hist_factory
def getConfigHist(config_factory, plot_group, selection, branch_name, 
    states, luminosity=1, cut_string=""):
    try:
        filelist = config_factory.getPlotGroupMembers(plot_group)
    except ValueError as e:
        logging.warning(e.message)
        logging.warning("Treating %s as file name" % plot_group)
        filelist = [plot_group]
    hist_info = getHistFactory(config_factory, selection, filelist, luminosity, cut_string)
    bin_info = config_factory.getHistBinInfo(branch_name)
    hist = ROOT.TH1F(plot_group, plot_group, bin_info['nbins'], bin_info['xmin'], bin_info['xmax'])
    log_info = ""
    for name, entry in hist_info.iteritems():
        producer = entry["histProducer"]
        log_info += "\n" + "-"*70 +"\nName is %s entry is %s" % (name, entry)
        for state in states:
            log_info += "\nFor state %s" % state
            config_factory.setProofAliases(state)
            draw_expr = config_factory.getHistDrawExpr(branch_name, name, state)
            logging.debug("Draw expression was %s" % draw_expr)
            proof_name = "-".join([name, "WZAnalysis-%s#/%s/final/Ntuple" % (selection, state)])
            try:
                state_hist = producer.produce(draw_expr, proof_name)
                log_info += "\nNumber of events: %f" % state_hist.Integral()
                hist.Add(state_hist)
            except ValueError as error:
                logging.warning(error)
                log_info += "\nNumber of events: 0.0" 
        log_info += "total number of events: %f" % hist.Integral()
        config_factory.setHistAttributes(hist, branch_name, plot_group)
    logging.debug(log_info)
    logging.debug("Hist has %i entries" % hist.GetEntries())
    return hist
def getPlotPaths(selection, write_log_file):
    if "hep.wisc.edu" in os.environ['HOSTNAME']:
        storage_area = "/nfs_scratch/kdlong"
        html_area = "/afs/hep.wisc.edu/home/kdlong/public_html"
    else:
        storage_area = "/data/kelong"
        html_area = "/afs/cern.ch/user/k/kelong/www"
    base_dir = "%s/WZAnalysisData/PlottingResults" % storage_area
    plot_path = "/".join([base_dir, selection, 
        '{:%Y-%m-%d}'.format(datetime.datetime.today()),
        '{:%Hh%M}'.format(datetime.datetime.today())])
    makeDirectory(plot_path)
    if write_log_file:
        makeDirectory("/".join([plot_path, "logs"]))
    html_path = plot_path.replace(storage_area, html_area)
    return (plot_path, html_path)
def savePlot(canvas, plot_path, html_path, branch_name, write_log_file, args):
    if args.output_file != "":
        canvas.Print(args.output_file)
        return
    if write_log_file:
        log_file = "/".join([plot_path, "logs", "%s_event_info.log" % branch_name])
        shutil.move("temp.txt", log_file) 
    canvas.Print("/".join([plot_path, branch_name + ".pdf"]))
    canvas.Print("/".join([plot_path, branch_name + ".root"]))
    if not args.no_html:
        makeDirectory(html_path)
        makeDirectory("/".join([html_path, "logs"]))
        canvas.Print("/".join([html_path, branch_name + ".pdf"]))
        if write_log_file:
            shutil.copy(log_file, log_file.replace(plot_path, html_path))
def makeDirectory(path):
    '''
    Make a directory, don't crash
    '''
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: 
            raise
