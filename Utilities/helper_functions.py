import ROOT
import plot_functions as plotter
import Utilities.WeightInfo as WeightInfo
from Utilities.WeightedHistProducer import WeightedHistProducer
from Utilities.FromFileHistProducer import FromFileHistProducer
from Utilities.ConfigHistFactory import ConfigHistFactory 
from collections import OrderedDict
import os
import subprocess
import glob
import logging
import datetime
import shutil
import errno
import math
import array
from IPython import embed

def makePlots(hist_stacks, data_hists, name, args, signal_stacks=[0], errors=[]):
    canvas_dimensions = [800, 800] if "unrolled" not in name else [1200, 800]
    canvas = ROOT.TCanvas("%s_canvas" % name, name, *canvas_dimensions) 
    first = True
    for hist_stack, data_hist, signal_stack in zip(hist_stacks, data_hists, signal_stacks):
        makePlot(hist_stack, data_hist, name, args, signal_stack, 
            same=(" same" if not first else ""))
        first = False
    offset = ROOT.gPad.GetLeftMargin() - 0.04 if args.legend_left else \
        ROOT.gPad.GetRightMargin() - 0.04 
    if hasattr(args, "selection"):
        width = .2 if "WZxsec" in args.selection else 0.33
    else: 
        width = .33
    width *= args.scalelegx
    xdist = 0.1 if args.legend_left else 0.91
    xcoords = [xdist+offset, xdist+width+offset] if args.legend_left \
        else [xdist-width-offset, xdist-offset]
    unique_entries = min(len(hist_stacks[0].GetHists()), 8)
    ymax = 0.8 if args.legend_left else 0.9
    ycoords = [ymax, ymax - 0.08*unique_entries*args.scalelegy]
    coords = [xcoords[0], ycoords[0], xcoords[1], ycoords[1]]
    
    if "none" not in args.uncertainties:
        histErrors = getHistErrors(hist_stacks[0], args.nostack) if not errors else errors
        for error_hist,signal_stack,data_hist in zip(histErrors, signal_stacks, data_hists):
            ROOT.SetOwnership(error_hist, False)
            error_hist.SetLineWidth(1)
            ROOT.gStyle.SetHatchesLineWidth(1)
            ROOT.gStyle.SetHatchesSpacing(0.75)
            error_hist.Draw("same e2")
            if signal_stack:
                signal_stack.Draw("nostack same hist")
            if not args.data == 'none':
                data_hist.Draw("e0 same")
            error_title = "Stat. unc."
            if "all" in args.uncertainties:
                error_title = "Stat.#oplusSyst."
            elif "scale" in args.uncertainties:
                error_title = "Stat.#oplusScale"
            error_hist.SetTitle(error_title)
    else:
        histErrors = []
    legend = getPrettyLegend(hist_stacks[0], data_hists[0], signal_stacks[0], histErrors, coords)
    legend.Draw()

    if not args.no_decorations and hasattr(ROOT, "CMSlumi"):
        scale_label = "Normalized to Unity" if args.luminosity < 0 else \
            "%0.1f fb^{-1}" % args.luminosity
        
        lumi_text = ""
        if args.thesis:
            lumi_text = "Thesis" 
        elif args.preliminary:
            lumi_text = "Preliminary" 
        if args.simulation:
            lumi_text = "Simulation " + lumi_text 

        ROOT.CMSlumi(canvas, 0, 11, "%s (13 TeV)" % scale_label, lumi_text)
    if args.extra_text != "":
        lines = [x.strip() for x in args.extra_text.split(";")]
        ymax = coords[3] - 0.02
        box_size = 0.05*len(lines)*args.scalelegy
        if args.extra_text_above:
            ymax = coords[1] 
            coords[1] -= box_size
            coords[3] -= box_size
        ymin = ymax - box_size
        text_box = ROOT.TPaveText(coords[0], ymin, coords[2], ymax, "NDCnb")
        text_box.SetFillColor(0)
        text_box.SetFillStyle(0)
        text_box.SetLineColor(0)
        text_box.SetTextFont(42)
        for i, line in enumerate(lines):
            text_box.AddText(line)
        text_box.Draw()
        ROOT.SetOwnership(text_box, False)
    if args.logy:
        canvas.SetLogy()
    if not args.no_ratio:
        canvas = plotter.splitCanvas(canvas, canvas_dimensions,
                "#scale[0.85]{Data / Pred.}" if data_hists[0] else args.ratio_text,
                [float(i) for i in args.ratio_range],
                args.noRatioUnc
        )
    return canvas
def makePlot(hist_stack, data_hist, name, args, signal_stack=0, same=""):
    canvas = ROOT.gROOT.FindObject("%s_canvas" % name)
    ROOT.SetOwnership(canvas, False)
    stack_signal = signal_stack != 0 and args.stack_signal
    stack_drawexpr = " ".join(["hist"] + 
        ["nostack" if args.nostack else ""]
    )
    hists = hist_stack.GetHists()
    first_stack = hist_stack

    hist_stack.Draw(stack_drawexpr + (same if "same" not in stack_drawexpr else ""))
    if signal_stack:
        if stack_signal:
            sum_stack = hists[0].Clone()
            for hist in hists[1:]:
                skip = False
                if args.exclude_from_sigstack:
                    for name in args.exclude_from_sigstack.split(","):
                        if name in hist.GetName():
                            skip = True
                if not skip:
                    sum_stack.Add(hist)
            for i,hist in enumerate(signal_stack.GetHists()):
                hist.Add(sum_stack)
        signal_stack.Draw("hist nostack same")
        signal_stack.GetHistogram().GetXaxis().SetTitle(
            hists[0].GetXaxis().GetTitle())

    if data_hist:
        if not "yield" in name.lower():
            data_hist.Sumw2(False)
            data_hist.SetBinErrorOption(ROOT.TH1.kPoisson)
        data_hist.Draw("e0 same")
    first_stack.GetYaxis().SetTitleSize(hists[0].GetYaxis().GetTitleSize())    
    first_stack.GetYaxis().SetTitleOffset(hists[0].GetYaxis().GetTitleOffset())    
    first_stack.GetYaxis().SetTitle(
        hists[0].GetYaxis().GetTitle())

    if not args.no_ratio and float(ROOT.gROOT.GetVersion().split("/")[0]) > 6.07:
        # Remove first bin label to avoid overlap of canvases
        if hists[0].GetMinimum() == 0.0:
            first_stack.GetYaxis().ChangeLabel(1, -1.0, 0)
    first_stack.GetHistogram().GetXaxis().SetTitle(
        hists[0].GetXaxis().GetTitle())
    first_stack.GetHistogram().SetLabelSize(0.04)
    first_stack.SetMinimum(hists[0].GetMinimum()*args.scaleymin)
    ## Adding an arbirary factor of 100 here so the scaling doesn't cut off info from
    ## Hists. Not applied when lumi < 1. This should be fixed
    scale = args.luminosity/100 if args.luminosity > 0 else 1
    first_stack.SetMaximum(hists[0].GetMaximum()*args.scaleymax*scale)
def getHistErrors(hist_stack, separate):
    histErrors = []
    for hist in hist_stack.GetHists():
        error_hist = plotter.getHistErrors(hist)
        if separate:
            error_hist.SetFillColor(hist.GetLineColor())
            error_hist.SetLineColor(hist.GetLineColor())
            histErrors.append(error_hist)
        else:
            if len(histErrors) == 0:
                histErrors.append(error_hist)
            else:
                histErrors[0].Add(error_hist)
    return histErrors
def getPrettyLegend(hist_stack, data_hist, signal_stack, error_hists, coords):
    hists = hist_stack.GetHists()
    if signal_stack != 0:
        #hists += signal_stack.GetHists()
        hists = signal_stack.GetHists() + hists
    legend = ROOT.TLegend(coords[0], coords[1], coords[2], coords[3])
    ROOT.SetOwnership(legend, False)
    legend.SetName("legend")
    legend.SetFillStyle(0)
    if data_hist:
        legend.AddEntry(data_hist, data_hist.GetTitle(), "lp")
    hist_names = []
    for hist in reversed(hists):
        if hist.GetTitle() not in hist_names:
            legend.AddEntry(hist, hist.GetTitle(), "f")
        hist_names.append(hist.GetTitle())
    for error_hist in error_hists:
        legend.AddEntry(error_hist, error_hist.GetTitle(), "f")
    return legend
def getHistFactory(config_factory, selection, filelist, luminosity=1, hist_file=None):
    if "Gen" not in selection:
        metaTree_name = "metaInfo/metaInfo"
        sum_weights_branch = "summedWeights"
        weight_branch = "genWeight"
    else:
        metaTree_name = "analyze%s/MetaData" % ("ZZ" if "ZZ" in selection else "WZ")
        sum_weights_branch = "initSumWeights"
        weight_branch = "weight"
    mc_info = config_factory.getMonteCarloInfo()
    all_files = config_factory.getFileInfo()
    hist_factory = OrderedDict() 
    for name in filelist:
        subtract = False
        if "-" in name[0]:
            subtract = True
            name = name[1:]
        base_name = name.split("__")[0]
        if name not in all_files.keys():
            if not hist_file is None and not hist_file.Get(name):
                logging.warning("%s is not a valid file name (must match a definition in FileInfo/%s.json)" % \
                    (name, selection))
                continue
            hist_factory[name] = dict(all_files[base_name])
        else:
            hist_factory[name] = dict(all_files[name])
        hist_factory[name]["subtract"] = subtract
        if "data" not in name.lower() and name != "nonprompt":
            ref_info = mc_info[name] if name in mc_info.keys() else mc_info[base_name]
            kfac = 1. if 'kfactor' not in ref_info.keys() else ref_info['kfactor']
            if not hist_file:
                metaTree = ROOT.TChain(metaTree_name)
                metaTree.Add(hist_factory[name]["file_path"])
                weight_info = WeightInfo.WeightInfoProducer(metaTree, 
                        ref_info['cross_section']*kfac,
                        sum_weights_branch).produce()
            else:
                # Not the most elegant way to go about it, but better to read
                # sum_of_weights from the histogram created by the anlysis code than to recalculate
                # it from the metadata just in case the file paths have changed
                sumweights_hist = hist_file.Get(str("/".join([name, "sumweights"])))
                if not sumweights_hist:
                    logging.warning("Using sumWeights from %s for file %s" % (base_name, name))
                    sumweights_hist = hist_file.Get(str("/".join([base_name, "sumweights"])))
                ROOT.SetOwnership(sumweights_hist, False)
                weight_info = WeightInfo.WeightInfo(
                        mc_info[base_name]['cross_section']*kfac,
                        sumweights_hist.Integral(0, sumweights_hist.GetNbinsX()+1) if sumweights_hist else 0
                )
        else:
            weight_info = WeightInfo.WeightInfo(1, 1)
            weight_branch = ""
        if not hist_file:
            histProducer = WeightedHistProducer(weight_info, weight_branch)  
        else:
            histProducer = FromFileHistProducer(weight_info, hist_file)  
        if "data" not in name.lower() and name != "nonprompt":
            histProducer.setLumi(luminosity)
        hist_factory[name].update({"histProducer" : histProducer})
        hist_factory[name].update({"configFactory" : config_factory})
        hist_factory[name].update({"fromFile" : hist_file is not None})
    return hist_factory
def getConfigHist(hist_factory, branch_name, bin_info, plot_group, selection, states, 
        uncertainties="none", addOverflow=False, rebin=0, cut_string="", removeNegatives=True, scaleType=False):
    hist_name = "_".join([plot_group, selection.replace("/", "_"), branch_name])
    # TODO: Understand why this is broken in newer ROOT versions
    #rootdir = "gProof" if hasattr(ROOT, "gProof") else "gROOT"
    #hist = getattr(ROOT, rootdir).FindObject(str(hist_name))
    hist = ROOT.gROOT.FindObject(hist_name)
    if hist:
        hist.Delete()
    if not hist_factory.itervalues().next()["fromFile"]:
        hist = ROOT.TH1D(hist_name, hist_name, bin_info['nbins'], bin_info['xmin'], bin_info['xmax'])
    log_info = "" 
    final_counts = {c : 0 for c in states}
    for name, entry in hist_factory.iteritems():
        log_info += "_"*80 + "\n"
        log_info += "Results for file %s in plot group %s\n" % (name, plot_group)
        producer = entry["histProducer"]
        config_factory = entry["configFactory"]
        weight = config_factory.getPlotGroupWeight(plot_group)
        scale = config_factory.getPlotGroupScale(plot_group)
        if weight != 1 and not entry["fromFile"]:
            producer.addWeight(weight)
        for state in states:
            if entry["fromFile"]:
                chan = state
                state_hist_name = str("%s/%s_%s" % (name, branch_name, state))
                if "dress" in plot_group: 
                    state_hist_name = state_hist_name.replace("_lhe", "")
                    state_hist_name = state_hist_name.replace("_barelep", "")
                    state_hist_name = state_hist_name.replace("_born", "")
                elif "bare" in plot_group: 
                    state_hist_name = state_hist_name.replace("_lhe", "_barelep")
                    state_hist_name = state_hist_name.replace("_born", "_barelep")
                elif "prefsr" in plot_group: 
                    state_hist_name = state_hist_name.replace("_lhe", "_prefsr")
                    state_hist_name = state_hist_name.replace("_born", "_prefsr")
                elif "lhe" in plot_group: 
                    state_hist_name = state_hist_name.replace("_prefsr", "_lhe")
                    state_hist_name = state_hist_name.replace("_born", "_prefsr")
                if "nonprompt" in str(plot_group).lower():
                    splitName = state_hist_name.rsplit(branch_name, 1)
                    state_hist_name = ''.join(splitName[:-1] +[branch_name, "_Fakes", splitName[-1]])
                scaleFac = config_factory.getPlotGroupScale(plot_group)
                args = [state_hist_name, addOverflow, rebin, scaleType, scale]
            else:
                chan = state.split("/")[0] if "ntuple" in state else ""
                config_factory.setProofAliases(chan)
                draw_expr = config_factory.getHistDrawExpr(branch_name, name, chan)
                proof_name = "_".join([name, "%s#/%s" % (selection.replace("/", "_"), state)])
                producer.setCutString(cut_string)
                args = [draw_expr, proof_name, addOverflow, rebin]
            
            try:
                state_hist = producer.produce(*args)
                if entry["subtract"]:
                    state_hist.Scale(-1.)
            except ValueError as error:
                logging.warning(error)
                logging.warning("Error for file %s: %s" % (name, error))
                log_info += "Number of events in %s channel: 0.0\n"  % state
                continue
            if not hist:
                hist = state_hist
                hist.SetName(hist_name)
                if state_hist.InheritsFrom("TH1"):
                    hist.SetTitle(hist_name)
            else:
                hist.Add(state_hist)
            final_counts[state] += state_hist.Integral() 
            log_info += "Number of events in %s channel: %0.2f\n" % (state, state_hist.Integral())
            log_info += "Number of entries is %i\n" % (state_hist.GetEntries() - addOverflow)
        log_info += "Total number of events: %0.2f\n" % (hist.Integral() if hist and hist.InheritsFrom("TH1") else 0)
        log_info += "Cross section is %0.4f\n" % producer.getCrossSection()
        log_info += "Sum of weights is %0.2f\n" % producer.getSumOfWeights() 
        log_info += "--> Scaling raw events by %0.4E*(scaleFac=%0.2f)\n" % (producer.getHistScaleFactor(), config_factory.getPlotGroupScale(plot_group))
    logging.debug(log_info)
    logging.debug("Hist has %i entries" % (hist.GetEntries() if hist and hist.InheritsFrom("TH1") else 0) )
    log_info += "*"*80 + "\n"
    log_info += "    Summary for plot group %s\n" % plot_group
    log_info += "    Total entries in all states: %0.2f\n" % (hist.Integral() if hist and hist.InheritsFrom("TH1") else 0)
    for state in states:
        log_info += "    %0.2f events in state %s\n" % (final_counts[state], state)
    log_info += "*"*80 + "\n"
    with open("temp-verbose.txt", "a") as log_file:
        log_file.write(log_info)
    if not hist or not hist.InheritsFrom("TH1"):
        raise ValueError("Invalid histogram %s for selection %s" % (branch_name, selection))
    if removeNegatives and "data" not in name:
        removeZeros(hist)
    return hist

def getConfigHistFromFile(filename, config_factory, plot_group, selection, branch_name, channels,
        luminosity=1, addOverflow=False, rebin=0, uncertainties="none", removeNegatives=True, scaleType=False):
    try:
        filelist = config_factory.getPlotGroupMembers(plot_group)
    except ValueError as e:
        logging.warning(e.message)
        logging.warning("Treating %s as file name" % plot_group)
        filelist = [plot_group]

    if branch_name not in config_factory.getListOfPlotObjects():
        raise ValueError("Invalid histogram %s for selection %s" % (branch_name, selection))
    
    # If reading from file, the weighted 1D hists need to be precomputed
    # Look for them stored in the file with the plot_group name
    weight = config_factory.getPlotGroupWeight(plot_group)
    if weight != 1:
        filelist = [str(plot_group)]
        logging.warning("Treating %s as file name" % plot_group)

    hist_file = ROOT.TFile(filename)
    ROOT.SetOwnership(hist_file, False)
    hist_factory = getHistFactory(config_factory, selection, filelist, luminosity, hist_file)

    bin_info = config_factory.getHistBinInfo(branch_name)
    states = channels.split(",")
    hist = getConfigHist(hist_factory, branch_name, bin_info, plot_group, selection, states, uncertainties, addOverflow, rebin, scaleType=scaleType)
    config_factory.setHistAttributes(hist, branch_name, plot_group)

    return hist

def getConfigHistFromTree(config_factory, plot_group, selection, branch_name, channels, blinding=[],
    luminosity=1, addOverflow=False, rebin=0, cut_string="", uncertainties="none"):
    if "Gen" not in selection:
        states = [x.strip() for x in channels.split(",")]
        scale_weight_expr = "scaleWeights/scaleWeights[0]"
        trees = ["%s/ntuple" % state for state in states]
    else:
        trees = ["analyze%s/Ntuple" % ("ZZ" if "ZZ" in selection else "WZ")]
        scale_weight_expr = "LHEweights/LHEweights[0]"
        if channels != "eee,mmm,eem,emm":
            chan_cuts = []
            for chan in channels.split(","): 
                chan_cuts.append(getGenChannelCut(chan))
            cut_string = appendCut(cut_string, " || ".join(chan_cuts))
    try:
        filelist = config_factory.getPlotGroupMembers(plot_group)
    except ValueError as e:
        logging.warning(e.message)
        logging.warning("Treating %s as file name" % plot_group)
        filelist = [plot_group]
    hist_factory = getHistFactory(config_factory, selection, filelist, luminosity)
    bin_info = config_factory.getHistBinInfo(branch_name)
    
def histWithScaleUnc(scale_hist2D, entries, name):
    if not isinstance(scale_hist2D, ROOT.TH2):
        raise ValueError("Scale uncertainties require 2D histogram")
    scale_hist = scale_hist2D.ProjectionY("temp", 1, 1, "e")
    scale_hist.SetName(name)
    ROOT.SetOwnership(scale_hist, False)
    hists = []
    for i in range(2, entries+1):
        hist = scale_hist2D.ProjectionY(name+"_weight%i" % i, i, i, "e")
        hist.Sumw2()
        hists.append(hist)
    # Choose max variation between scale choices by bin
    for i in range(1, scale_hist.GetNbinsX()+1):
        try:
            maxScale = max([h.GetBinContent(i) for h in hists \
                    if not h.GetBinContent(i) == 0])
            minScale = min([h.GetBinContent(i) for h in hists \
                    if not h.GetBinContent(i) == 0])
            scaleUp_diff = maxScale - scale_hist.GetBinContent(i)
            scaleDown_diff = scale_hist.GetBinContent(i) - minScale
            maxScaleErr = max(scaleUp_diff, scaleDown_diff)
        except:
            maxScaleErr = 0
        # Just symmetric errors for now
        err = math.sqrt(scale_hist.GetBinError(i)**2 + maxScaleErr**2)
        scale_hist.SetBinError(i, err)
    return scale_hist 

def appendCut(cut_string, add_cut):
    if cut_string != "" and add_cut not in cut_string:
        append_cut = lambda x: "*(%s)" % x if x not in ["", None] else x
        return "(" + cut_string + ")" + append_cut(add_cut)
    else:
        return add_cut

def getPlotPaths(selection, folder_name, write_log_file=False):
    if "hep.wisc.edu" in os.environ['HOSTNAME']:
        storage_area = "/nfs_scratch/kdlong"
        html_area = "/afs/hep.wisc.edu/home/kdlong/public_html"
    else:
        storage_area = "/eos/user/k/%s" % os.environ["USER"]
        html_area = "/eos/user/k/%s/www" % os.environ["USER"]
    base_dir = "%s/PlottingResults" % storage_area
    plot_path = "/".join([base_dir, selection] +
       (['{:%Y-%m-%d}'.format(datetime.datetime.today()),
        '{:%Hh%M}'.format(datetime.datetime.today())] if folder_name == "" \
            else [folder_name])
    )
    makeDirectory(plot_path + "/plots")
    if write_log_file:
        makeDirectory(plot_path + "/logs")
    html_path = plot_path.replace(storage_area, html_area)
    return (plot_path, html_path)

def getGenChannelCut(channel):
    cut_string = ""
    if channel == "eem":
        cut_string = "((abs(l1pdgId) == 11 && abs(l2pdgId) == 11 && abs(l3pdgId) == 13)" \
            " || (abs(l1pdgId) == 11 && abs(l2pdgId) == 13 && abs(l3pdgId) == 11)" \
            " || (abs(l1pdgId) == 13 && abs(l2pdgId) == 11 && abs(l3pdgId) == 11))"
    elif channel == "emm":
        cut_string = "((abs(l1pdgId) == 13 && abs(l2pdgId) == 13 && abs(l3pdgId) == 11)" \
            " || (abs(l1pdgId) == 13 && abs(l2pdgId) == 11 && abs(l3pdgId) == 13)" \
            " || (abs(l1pdgId) == 11 && abs(l2pdgId) == 13 && abs(l3pdgId) == 13))"
    elif channel == "eee":
        cut_string = "(abs(l1pdgId) == 11 && abs(l2pdgId) == 11 && abs(l3pdgId) == 11)"
    elif channel == "mmm":
        cut_string = "(abs(l1pdgId) == 13 && abs(l2pdgId) == 13 && abs(l3pdgId) == 13)"
    return cut_string

def savePlot(canvas, plot_path, html_path, branch_name, write_log_file, args):
    canvas.Update()
    if args.output_file != "":
        canvas.Print(args.output_file)
        return
    if write_log_file:
        log_file = "/".join([plot_path, "logs", "%s_event_info.log" % branch_name])
        verbose_log = log_file.replace("event_info", "event_info-verbose")
        shutil.move("temp.txt", log_file) 
        if os.path.isfile("temp-verbose.txt"):
            shutil.move("temp-verbose.txt", verbose_log) 
    output_name ="/".join([plot_path, "plots", branch_name]) 
    canvas.Print(output_name + ".root")
    canvas.Print(output_name + ".C")
    if not args.no_html:
        coverteps = subprocess.call(["command", "-v", "epstopdf"])
        makeDirectory(html_path + "/plots")
        output_name ="/".join([html_path, "plots", branch_name])
        canvas.Print(output_name + ".png")
        # Because ROOT pdf output doesn't display latex properly.
        # Unfortunately you can't take epstopdf for granted
        if coverteps == 0:
            canvas.Print(output_name + ".eps")
            subprocess.call(["epstopdf", "--outfile=%s" % output_name+".pdf", output_name+".eps"])
            os.remove(output_name+".eps")
        else:
            canvas.Print(output_name + ".pdf")
        if write_log_file:
            makeDirectory(html_path + "/logs")
            shutil.copy(log_file, log_file.replace(plot_path, html_path))
            if os.path.isfile(verbose_log):
                shutil.copy(verbose_log, verbose_log.replace(plot_path, html_path))
    del canvas

# Leave zero events untouched, fix negative yields (for nonprompt)
def removeZeros(hist):
    for i in range(hist.GetNbinsX()+2):
        err = hist.GetBinError(i)
        if hist.GetBinContent(i) < 0:
            if "Up" in hist.GetName():
                hist.SetBinContent(i, 0.0001)
            elif "Down" in hist.GetName():
                hist.SetBinContent(i, 0.00001)
            else: 
                hist.SetBinContent(i, 0.001)
            hist.SetBinError(i, err)

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
