from ROOT import *
from utils import *
import os,sys
#gROOT.SetBatch(1)
gStyle.SetOptStat(0)
gROOT.ForceStyle()

lumi = 1.0#/fb
signal_mag = 10000
try: infile = sys.argv[1]
except: infile = 'hists-Summer16.TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1AOD_103*_RA2AnalysisTree.root'

binningAnalysis = binning
redoBinning = binningAnalysis


fileA = TFile(infile)
keys = fileA.GetListOfKeys()


newfile = TFile('plots.root','recreate')


bkghistlists = {}
sighistlists = {}
canv = mkcanvas('canv')
#canv.SetLogy()

for key in keys:
	name = key.GetName()
	if not 'G_' in name: continue
	kinvar = name.split('_')[1].replace('Observed','')
	if '_' in kinvar: continue
	histg = fileA.Get(name)
	histoStyler(histg, kCyan-8)
	histg.SetFillColor(histg.GetLineColor())
	histg.SetFillStyle(1001)

	histq = fileA.Get(name.replace('G_','Q_'))
	histoStyler(histq, kGreen-5)
	histq.SetFillColor(histq.GetLineColor())
	histq.SetFillStyle(1001)
	histq.SetLineColor(kGray)
			
	histq.Add(histg)
	

	leg = mklegend(x1=.53, y1=.76, x2=.88, y2=.88, color=kWhite)
	
	histq.GetYaxis().SetTitle('# of jets')
	
	histq.GetYaxis().SetRangeUser(0.00001,1.5*histq.GetMaximum())
	histq.SetTitle('')
	histq.GetXaxis().SetTitle(kinvar)
	histq.Draw('hist')	
	histg.Draw('hist same')	
	histq.Draw('axis same')		
	leg.AddEntry(histg, 'gluon jets (t#bar{t}+jets)', 'lpf')
	leg.AddEntry(histq, 'quark jets (t#bar{t}+jets)', 'lpf')	
	leg.Draw()
	tl.SetTextSize(0.7*tl.GetTextSize())
	tl.DrawLatex(0.55,0.7,'passes W#rightarrow jj selection')
	tl.SetTextSize(1.0/0.7*tl.GetTextSize())
	stamp('n/a')
	newfile.cd()
	canv.Update()
	pause()
	#canv.Write('c_'+name)	
	#canv.Print('pdfs/'+name[1:]+'.pdf')
	

print 'just created', newfile.GetName()
newfile.Close()





exit(0)




















