import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptFit(111)

#Pads
ROOT.gStyle.SetPadBorderMode(0)
ROOT.gStyle.SetCanvasBorderMode(0)

#Margins
ROOT.gStyle.SetPadTopMargin(0.09)
ROOT.gStyle.SetPadBottomMargin(0.13)
ROOT.gStyle.SetPadLeftMargin(0.13)
ROOT.gStyle.SetPadRightMargin(0.05)

#Text size and font
#ROOT.gStyle.SetTitleSize(0.06)
#ROOT.gStyle.SetTitleFont(42)
#ROOT.gStyle.SetLabelSize(0.06)
#ROOT.gStyle.SetOptTitle(0)
#ROOT.gStyle.SetTitleFontSize(0.05)

#legends
topRight_legend_relEresol = ROOT.TLegend(0.45, 0.75, 0.9, 0.9)
topRight_legend_relEresol.SetBorderSize(0)
topRight_legend_relEresol.SetFillStyle(0)

topRight_legend = ROOT.TLegend(0.65, 0.75, 1.15, 0.9)
topRight_legend.SetBorderSize(0)
topRight_legend.SetFillStyle(0)

#Colors
ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)
colors = [40, 41, 46, 30, 38, 43, 37, 29] #range(40, 50)
