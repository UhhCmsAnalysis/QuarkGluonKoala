#Welcome to the industrial age of Sam's rebalance and smear code. You're going to have a lot of fun!
import os,sys
from ROOT import *
from array import array
from glob import glob
from utils import *
#from ra2blibs import *
import time

##read in command line arguments
defaultInfile_ = "/pnfs/desy.de/cms/tier2/store/user/vormwald/NtupleHub/ProductionRun2v3/Summer16.TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8AOD_0_RA2AnalysisTree.root"
defaultInfile_ = "/pnfs/desy.de/cms/tier2/store/user/vormwald/NtupleHub/ProductionRun2v3/Summer16.TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1AOD_103*_RA2AnalysisTree.root"
#T2qqGG.root
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbosity", type=int, default=1,help="analyzer script to batch")
parser.add_argument("-printevery", "--printevery", type=int, default=1000,help="short run")
parser.add_argument("-fin", "--fnamekeyword", type=str,default=defaultInfile_,help="file")
parser.add_argument("-quickrun", "--quickrun", type=bool, default=False,help="short run")
args = parser.parse_args()
fnamekeyword = args.fnamekeyword
inputFiles = glob(fnamekeyword)
verbosity = args.verbosity
printevery = args.printevery
quickrun = args.quickrun
if quickrun: n2process = 10000
else: n2process = 9999999999999

if 'Run2016' in fnamekeyword or 'Summer16' in fnamekeyword: 
    BTAG_CSVv2 = 0.8484
    BTAG_deepCSV = 0.6324

#Dictionary list of region selection sets
regionCuts = {}
pi = 3.14159
Inf = 9999
#varlist                    = ['Ht',      'Met',    'NLeptons'   'NJets', 'BTags'   'Mjj'        'PartonFlav']
regionCuts['NoCuts']        = [[0,Inf],   [0,Inf],  [1,Inf],     [0,Inf], [0,Inf],  [-Inf,Inf], [-Inf,Inf]]
regionCuts['qqbarControl']  = [[100,Inf], [50,Inf], [1,1],       [4,4],   [2,2],    [60,100],   [-Inf,Inf]]
regionCuts['qqbarControlQ'] = [[100,Inf], [50,Inf], [1,1],       [4,4],   [2,2],    [60,100],   [1,5]]
regionCuts['qqbarControlG'] = [[100,Inf], [50,Inf], [1,1],       [4,4],   [2,2],    [60,100],   [21,21]]
    
    
##declare and load a tree
c = TChain('TreeMaker2/PreSelection')
for fname in inputFiles: 
	print 'adding', fname, 'to tree'
	c.Add(fname)
nentries = c.GetEntries()
c.Show(0)
n2process = min(n2process, nentries)
print 'n(entries) =', n2process

varlist = ['Ht', 'Met','NLeptons', 'NJets', 'BTags', 'Mjj', 'PartonFlav']
indexVar = {}
for ivar, var in enumerate(varlist): indexVar[var] = ivar
indexVar[''] = -1
nmain = len(varlist)

histoStructDict = {}
for region in regionCuts:
	for var in varlist:
		histname = region+'_'+var
		histoStructDict[histname] = mkHistoStruct(histname, binning)

def selectionFeatureVector(fvector, regionkey='', omitcuts='', omitcuts_dphi=''):
	iomits, iomits_dphi = [], []  
	for cut in omitcuts.split('Vs'): iomits.append(indexVar[cut])
	for i, feature in enumerate(fvector):
		if i==nmain: break
		if i in iomits: continue
		if not (feature>=regionCuts[regionkey][i][0] and feature<=regionCuts[regionkey][i][1]): 
			return False
	return True

##Create output file
infileID = fnamekeyword.split('/')[-1].replace('.root','')
newname = 'hists-'+infileID+'.root'
print 'creating file', newname
fnew = TFile(newname, 'recreate')

hHt = TH1F('hHt','hHt',120,0,2500)
hHt.Sumw2()
hHtWeighted = TH1F('hHtWeighted','hHtWeighted',120,0,2500)
hHtWeighted.Sumw2()

t0 = time.time()
for ientry in range(n2process):

	if ientry%printevery==0:
		print "processing event", ientry, '/', n2process, 'time', time.time()-t0
	c.GetEntry(ientry)

	#br = 0.33
	weight = c.CrossSection
	
	fillth1(hHt, c.madHT, 1)
	fillth1(hHtWeighted, c.madHT, weight)	
	
	recomuons = []
	#build up the vector of jets using TLorentzVectors; 
	for imu, mu in enumerate(c.Muons):
		if not mu.Pt()>30: continue
		if not abs(mu.Eta())<2.4: continue
		tlvmu = TLorentzVector()
		tlvmu.SetPtEtaPhiE(mu.Pt(), mu.Eta(), mu.Phi(), mu.E())
		recomuons.append(tlvmu)	
		
	recoelectrons = []
	#build up the vector of jets using TLorentzVectors; 
	for imu, mu in enumerate(c.Electrons):
		if not mu.Pt()>30: continue
		if not abs(mu.Eta())<2.4: continue
		tlvmu = TLorentzVector()
		tlvmu.SetPtEtaPhiE(mu.Pt(), mu.Eta(), mu.Phi(), mu.E())
		recoelectrons.append(tlvmu)		
		
	#build up the vector of jets using TLorentzVectors; 
	#this is where you have to interface with the input format you're using
	recojets = []
	lightjets = []
	nb = 0
	for ijet, jet in enumerate(c.Jets):
		if not jet.Pt()>20: continue
		if not abs(jet.Eta())<2.4: continue
		if c.Jets_bJetTagDeepCSVBvsAll[ijet]>BTAG_deepCSV: nb+=1
		elif c.Jets_bJetTagDeepCSVBvsAll[ijet]<0.5: lightjets.append([jet, ijet])
		tlvjet = TLorentzVector()
		tlvjet.SetPtEtaPhiE(jet.Pt(), jet.Eta(), jet.Phi(), jet.E())
		recojets.append(tlvjet)
	
	if len(lightjets)>1: 
		mll = (lightjets[0][0]+lightjets[1][0]).M()
		partonflav = abs(c.Jets_partonFlavor[lightjets[0][1]])
	elif len(lightjets)>0:
		mll = -101.0
		partonflav = abs(c.Jets_partonFlavor[lightjets[0][1]])
	else: 
		mll = -101.0
		partonflav = -101.0
	
	fv = [c.HT, c.MET, len(recoelectrons)+len(recomuons), len(recojets), nb, mll, partonflav]
	for regionkey in regionCuts:
		for ivar, varname in enumerate(varlist):
			hname = regionkey+'_'+varname
			if selectionFeatureVector(fv,regionkey,varname,''): 
				if regionkey!='NoCuts' and varname=='PartonFlav': print 'filling ', regionkey, 'with', varname, fv[ivar]
				fillth1(histoStructDict[hname].Observed, fv[ivar], weight)


fnew.cd()
hHt.Write()
hHtWeighted.Write()
writeHistoStruct(histoStructDict)

print 'just created', fnew.GetName()
fnew.Close()








