import sys
import time
from ROOT import *
from utils import *
from glob import glob
from random import shuffle
import numpy as np
datamc = 'mc'

try: fnames = sys.argv[1]
except: 
	fnames = '/pnfs/desy.de/cms/tier2/store/user/*/NtupleHub/ProductionRun2v3/Summer16.TTJets_TuneCUETP8M1*70*'
try: pickprob_ = bool(sys.argv[2])
except: pickprob_ = False

doGluonTagging = True      

if pickprob_:
	fPickProbability = TFile('ScaleFactors/PickProbability.root')
	hPickProbability = fPickProbability.Get('hqog')
	import random    

files = glob(fnames)
c = TChain('TreeMaker2/PreSelection')
if datamc == 'mc': 
	for fname in files[:30]:
		accessname = fname.replace('/eos/uscms','root://cmsxrootd.fnal.gov/')
		print 'adding', accessname
		c.Add(accessname)

shortinputname = (fnames.split('/')[-1]+'.root').replace('.root.root','.root')
fnew_ = TFile('treeQG_'+shortinputname,'recreate')


hHt = TH1F('hHt','hHt',100,0,500)
hHtWeighted = TH1F('hHtWeighted','hHtWeighted',100,0,500)

var_pt = np.zeros(1,dtype=float)
var_eta = np.zeros(1,dtype=float)
var_qgLikelihood = np.zeros(1,dtype=float)
var_qgPtD = np.zeros(1,dtype=float)
var_qgMult = np.zeros(1,dtype=int)
var_qgAxisMaj = np.zeros(1,dtype=float)
var_qgAxisMin = np.zeros(1,dtype=float)
var_partonFlavor___ = np.zeros(1,dtype=int)
var_partonFlavor___[0] = -99

tEvent = TTree('tEvent','tEvent')
tEvent.Branch('pt', var_pt,'pt/D')
tEvent.Branch('eta', var_eta,'eta/D')
tEvent.Branch('qgLikelihood', var_qgLikelihood,'qgLikelihood/D')
tEvent.Branch('qgMult', var_qgMult,'qgMult/I')
tEvent.Branch('qgPtD', var_qgPtD,'qgPtD/D')
tEvent.Branch('qgAxisMaj', var_qgAxisMin,'qgAxisMaj/D')
tEvent.Branch('qgAxisMin', var_qgAxisMin,'qgAxisMin/D')
tEvent.Branch('partonFlavor', var_partonFlavor___,'partonFlavor/I')

c.Show(0)
nEvents = c.GetEntries()
verbosity = round(10000)

import time
t1 = time.time()
for ientry_ in range(c.GetEntries()):

	if ientry_<2: 
		import time
		t0 = time.time()    
	if ientry_%verbosity==0:
		import time
		t1 = time.time()
		rtime1 = (time.time()-t0)*nEvents/(ientry_+1)*(1-1.0*ientry_/nEvents), 'sec'
		print 'processing', ientry_,'/',nEvents, ' = ', 1.0*ientry_/nEvents,'remaining time estimate: ', rtime1 		

	c.GetEntry(ientry_) 
	hHt.Fill(c.HT)

	weight = c.CrossSection
	hHtWeighted.Fill(c.HT,weight)

	if not c.JetID: continue
	if not len(c.Electrons)==0: continue
	if not len(c.Muons)==0: continue     
	if not c.BTags==2: continue 
	for ijet, jet in enumerate(c.Jets):
		if not jet.Pt()>30: continue
		if not abs(jet.Eta())<2.4: continue
		#if c.Jets_bDiscriminatorCSV[ijet]>0.8484: continue
		if not c.Jets_bDiscriminatorCSV[ijet]<0.5: continue  
		var_qgLikelihood[0] = c.Jets_qgLikelihood[ijet]
		var_qgPtD[0] =  c.Jets_ptD[ijet]
		var_qgMult[0] =  c.Jets_multiplicity[ijet]
		var_qgAxisMaj[0] =  c.Jets_axismajor[ijet]
		var_qgAxisMin[0] =  c.Jets_axisminor[ijet]		
		var_partonFlavor___[0] = c.Jets_partonFlavor[ijet]
		if doGluonTagging:
			if not c.Jets_partonFlavor[ijet]==21: continue
		else:
			if not abs(c.Jets_partonFlavor[ijet])<5: continue		
		tEvent.Fill()



fnew_.cd()
hHt.Write()
hHtWeighted.Write()
tEvent.Write()
print 'just created', fnew_.GetName()
fnew_.Close()
