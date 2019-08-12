import sys
import time
from ROOT import *
from utils import *
from glob import glob
from random import shuffle
import numpy as np
datamc = 'mc'
mW = 80.4
DM = 1000

#arg1: filename #arg2: True if doing weights

try: fnames = sys.argv[1]
except: 
    fnames = '/pnfs/desy.de/cms/tier2/store/user/sbein/CommonNtuples/Summer16.QCD*'
    fnames = '/pnfs/desy.de/cms/tier2/store/user/sbein/CommonNtuples/Summer16.TTJets_TuneCUETP8M1*'
try: pickprob_ = bool(sys.argv[2])
except: pickprob_ = False

needsMetTrigger = False
if  'TTJets' in fnames:
    doQuarkTagging = True
    doGluonTagging = False
elif 'MET' in fnames: 
    doQuarkTagging = True
    doGluonTagging = False
    needsMetTrigger = True
elif 'QCD' in fnames:
    doQuarkTagging = False        
    doGluonTagging = True   
elif 'JetHT' in fnames:
    doQuarkTagging = False        
    doGluonTagging = True      
      

if pickprob_:
    fPickProbability = TFile('ScaleFactors/PickProbability.root')
    hPickProbability = fPickProbability.Get('hqog')
    import random    

files = glob(fnames)
c = TChain('TreeMaker2/PreSelection')
if datamc == 'mc': 
    for fname in files[:1]:
        accessname = fname.replace('/eos/uscms','root://cmsxrootd.fnal.gov/')
        print 'adding', accessname
        c.Add(accessname)

shortinputname = (fnames.split('/')[-1]+'.root').replace('.root.root','.root')
fnew = TFile('treeTagProbe_'+shortinputname,'recreate')


hHt = TH1F('hHt','hHt',100,0,500)
hHtWeighted = TH1F('hHtWeighted','hHtWeighted',100,0,500)

var_ht_ = np.zeros(1,dtype=float)
var_pt = np.zeros(1,dtype=float)
var_eta = np.zeros(1,dtype=float)
var_mjj = np.zeros(1,dtype=float)
var_dphi = np.zeros(1,dtype=float)
var_deta = np.zeros(1,dtype=float)
var_qgLikelihood = np.zeros(1,dtype=float)
var_qgPtD = np.zeros(1,dtype=float)
var_qgMult = np.zeros(1,dtype=int)
var_qgAxis2 = np.zeros(1,dtype=float)
var_partonFlavor___ = np.zeros(1,dtype=int)
var_partonFlavor___[0] = -99

tEvent = TTree('tEvent','tEvent')
tEvent.Branch('ht', var_ht_,'ht/D')
tEvent.Branch('pt', var_pt,'pt/D')
tEvent.Branch('eta', var_eta,'eta/D')
tEvent.Branch('mjj', var_mjj,'mjj/D')
tEvent.Branch('dphi', var_dphi,'dphi/D')
tEvent.Branch('deta', var_deta,'deta/D')
tEvent.Branch('qgLikelihood', var_qgLikelihood,'qgLikelihood/D')
tEvent.Branch('qgMult', var_qgMult,'qgMult/I')
tEvent.Branch('qgPtD', var_qgPtD,'qgPtD/D')
tEvent.Branch('qgAxis2', var_qgAxis2,'qgAxis2/D')      
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

    #hHtWeighted.Fill(c.HT,c.CrossSection)
    #weight = c.CrossSection

    probes = []
    probes_indices = []

    if doQuarkTagging:
        if needsMetTrigger:
            if not c.MET>200: continue
            for itrig, trigname in enumerate(c.TriggerNames):
                if not (bool(c.TriggerPass[46]) or bool(c.TriggerPass[47]) or bool(c.TriggerPass[48])):
                    continue
        #if not c.BTags==2: continue
        electrons_ = []
        for electron in c.Electrons: 
            if not electron.Pt()>30: continue
            if not abs(electron.Eta())<2.4: continue
            electrons_.append(electron)
        muons = []
        for muon in c.Muons: 
            if not muon.Pt()>30: continue
            if not abs(muon.Eta())<2.4: continue	
            muons.append(muon)
        if not (len(electrons_)+len(muons)==1 and len(electrons_)*len(muons)==0): continue
        lepton = (electrons_+muons)[0]
        lightjets = []
        lightjets_indices = []       
        bjets = [] 
        for ijet, jet in enumerate(c.Jets):
            if not jet.Pt()>30: continue
            if not abs(jet.Eta())<2.4: continue
            if not jet.DeltaR(lepton)>0.4: continue
            if c.Jets_bDiscriminatorCSV[ijet]>0.5:
                bjets.append(jet)
            else: continue
            lightjets.append(jet)
            lightjets_indices.append(ijet)
        if not len(bjets)==2: continue    
        if not len(lightjets)==2: continue            
        for ij in range(len(lightjets)):
            for jj in range(ij):
                p = (lightjets[ij]+lightjets[jj]).P()
                if not lightjets[ij].Angle(lightjets[jj].Vect())<1.3*TMath.ACos(1 - 2/pow(TMath.Sqrt(pow(p,2)+pow(mW,2))/mW,2)): continue
                if not abs((lightjets[ij]+lightjets[jj]).M()-mW)<DM: continue                
                if not lightjets_indices[ij] in probes_indices:
                    probes.append(lightjets[ij])
                    probes_indices.append(lightjets_indices[ij])
                if not lightjets_indices[jj] in probes_indices:
                    probes.append(lightjets[jj])
                    probes_indices.append(lightjets_indices[jj])
        if not len(probes)==2: continue                                        


    if doGluonTagging:
        if not c.JetID: continue
        if not c.MET<50: continue
        if not len(c.Electrons)==0: continue
        if not len(c.Muons)==0: continue       
        if not c.BTags==2: continue 
        bjets = []
        for ijet, jet in enumerate(c.Jets):
            if not jet.Pt()>30: continue
            if not abs(jet.Eta())<2.4: continue
            if c.Jets_bDiscriminatorCSV[ijet]>0.8484: bjets.append(jet)
            if not c.Jets_bDiscriminatorCSV[ijet]<0.5: continue  
            probes.append(jet.Clone())
            probes_indices.append(ijet)  
        if not len(probes)==1: continue    
        lightJetIsCloseToB = False
        for bjet in bjets:
            if bjet.DeltaR(probes[0])<0.9:
                lightJetIsCloseToB = True
                break
        if not lightJetIsCloseToB: continue



    for ijet, jet in enumerate(probes):
        if not jet.Pt()>150: continue
        if not jet.Pt()<500: continue
        if pickprob_:
            interp = hPickProbability.Interpolate(jet.Pt(),abs(jet.Eta()))
            #print 'interpolated value', interp
            if not random.random()<interp: continue
        var_pt[0] = jet.Pt()
        var_eta[0] = jet.Eta()
        try: 
            var_mjj[0] = (probes[0]+probes[1]).M()
            var_dphi[0] = probes[0].DeltaPhi(probes[1]) 
            var_deta[0] = abs(probes[0].Eta()-probes[1].Eta())                               
        except: 
            var_mjj[0] = -1.0
            var_dphi[0] = -1.0
            var_deta[0] =  -1.0
        var_qgLikelihood[0] = c.Jets_qgLikelihood[probes_indices[ijet]]
        var_qgPtD[0] =  c.Jets_qgPtD[probes_indices[ijet]]
        var_qgMult[0] =  c.Jets_qgMult[probes_indices[ijet]]
        var_qgAxis2[0] =  c.Jets_qgAxis2[probes_indices[ijet]]
        var_ht_[0] = c.HT
        try:var_partonFlavor___[0] = c.Jets_partonFlavor[probes_indices[ijet]]
        except: pass
        tEvent.Fill()



fnew.cd()
hHt.Write()
hHtWeighted.Write()
tEvent.Write()
print 'just created', fnew.GetName()
fnew.Close()
