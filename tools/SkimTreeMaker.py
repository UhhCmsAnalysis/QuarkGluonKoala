#! /usr/bin/env python
# script to create trees with track variables
# created May 3, 2017 -Sam Bein 

from ROOT import *
from utils import *
import os, sys
import numpy as np
from glob import glob

smdir = '/nfs/dust/cms/user/beinsam/CommonNtuples/MC_SM/'
smdir = '/pnfs/desy.de/cms/tier2/store/user/sbein/CommonNtuples/'
try: infilenames = sys.argv[1]
except: infilenames = smdir+'Summer16.WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1_97_RA2AnalysisTree.root'

newfilename = 'tree'+(infilenames.split('/')[-1]).replace('*','')
fnew = TFile(newfilename,'recreate')

hHt = TH1F('hHt','hHt',100,0,3000)
hHtWeighted = TH1F('hHtWeighted','hHtWeighted',100,0,3000)
histoStyler(hHt,kBlack)

var_Met = np.zeros(1,dtype=float)
var_Mht = np.zeros(1,dtype=float)
var_Mt2 = np.zeros(1,dtype=float)
var_St = np.zeros(1,dtype=float)
var_Ht = np.zeros(1,dtype=float)
var_MinDeltaPhiMetJets = np.zeros(1,dtype=float)
var_NJets = np.zeros(1,dtype=int)
var_BTags = np.zeros(1,dtype=int)
var_NLeptons = np.zeros(1,dtype=int)
var_NPhotons = np.zeros(1,dtype=int)
var_NShortTags = np.zeros(1,dtype=int)
var_NMediumTags = np.zeros(1,dtype=int)
var_NLongTags = np.zeros(1,dtype=int)
var_NTags = np.zeros(1,dtype=int)
var_DPhiMetSumTags = np.zeros(1,dtype=float)
var_SumTagPtOverMet = np.zeros(1,dtype=float)
var_CrossSection = np.zeros(1,dtype=float)

tEvent = TTree('tEvent','tEvent')
tEvent.Branch('Met', var_Met,'Met/D')
tEvent.Branch('Mht', var_Mht,'Mht/D')
tEvent.Branch('Mt2', var_Mt2,'Mt2/D')
tEvent.Branch('St', var_St,'St/D')
tEvent.Branch('Ht', var_Ht,'Ht/D')
tEvent.Branch('MinDeltaPhiMetJets', var_MinDeltaPhiMetJets,'MinDeltaPhiMetJets/D')
tEvent.Branch('NJets', var_NJets,'NJets/I')
tEvent.Branch('BTags', var_BTags,'BTags/I')
tEvent.Branch('NLeptons', var_NLeptons,'NLeptons/I')
tEvent.Branch('NPhotons', var_NPhotons,'NPhotons/I')
tEvent.Branch('NShortTags', var_NShortTags,'NShortTags/I')
tEvent.Branch('NMediumTags', var_NMediumTags,'NMediumTags/I')
tEvent.Branch('NLongTags', var_NLongTags,'NLongTags/I')
tEvent.Branch('NTags', var_NTags,'NTags/I')
tEvent.Branch('DPhiMetSumTags', var_DPhiMetSumTags,'DPhiMetSumTags/D')
tEvent.Branch('SumTagPtOverMet', var_SumTagPtOverMet,'SumTagPtOverMet/D')
tEvent.Branch('CrossSection', var_CrossSection,'CrossSection/D')


c = TChain('TreeMaker2/PreSelection')
filenamelist = glob(infilenames)
print 'adding', filenamelist
filelist = []
for filename in filenamelist:
    fname = filename.strip()
    c.Add(fname)

c.Show(0)
nentries = c.GetEntries()
print 'will analyze', nentries
verbosity = 1000

for ientry in range(nentries):

    if ientry%verbosity==0: 
        print 'analyzing event %d of %d' % (ientry, nentries)+ '....%f'%(100.*ientry/nentries)+'%'

    c.GetEntry(ientry)
    weight = c.CrossSection
    hHt.Fill(c.madHT)
    hHtWeighted.Fill(c.madHT, weight)


    if not (c.MET>150): continue
    if not (c.NJets>0): continue

    if 'TTJets_TuneCUET' in infilenames:
        if not c.madHT<600: continue
    elif 'TTJets_HT' in infilenames:
        if not c.madHT>600: continue        


    var_Met[0] = c.MET
    var_Mht[0] = c.MHT    
    var_Mt2[0] = c.MT2	
    var_Ht[0] = c.HT

    metvec = TLorentzVector()
    metvec.SetPtEtaPhiE(c.MET, 0, c.METPhi, c.MET)
    mindphi = 9999
    for jet in c.Jets[:2]:
        if not (abs(jet.Eta())<2.4 and jet.Pt()>30): continue
        if abs(jet.DeltaPhi(metvec))<mindphi:
            mindphi = abs(jet.DeltaPhi(metvec))
    var_MinDeltaPhiMetJets[0] = mindphi

    var_NJets[0] = c.NJets
    var_BTags[0] = c.BTags
    var_NLeptons[0] = len(c.Electrons)+len(c.Muons)

    sumtagvec = TLorentzVector()
    taggedtracks = []
    nshort = 0
    nmedium = 0
    nlong = 0

    for ichi, chiCand in enumerate(c.chiCands):

        if not (chiCand.Pt()>15 and abs(chiCand.Eta())<2.4): continue            
        phits = c.chiCands_nValidPixelHits[ichi]
        thits = c.chiCands_nValidTrackerHits[ichi]
        tlayers = c.chiCands_trackerLayersWithMeasurement[ichi]			
        short = phits>0 and thits==phits
        medium = tlayers< 7 and (thits-phits)>0
        long   = tlayers>=7 and (thits-phits)>0        
        if (medium or long):
            if not (c.chiCands_nMissingOuterHits[ichi]>=2): continue
        if short: passesDXY = abs(c.chiCands_dxyVtx[ichi])<0.02
        else: passesDXY = abs(c.chiCands_dxyVtx[ichi])<0.01
        if not passesDXY: continue		
        if not abs(c.chiCands_dzVtx[ichi])<0.05: continue
        neutralIso = c.chiCands_neutralPtSum[ichi]/c.chiCands[ichi].Pt()
        if not (c.chiCands_neutralPtSum[ichi]<=10 and neutralIso<=0.1): continue
        chargedIso = c.chiCands_chargedPtSum[ichi]/c.chiCands[ichi].Pt()
        if not (c.chiCands_chargedPtSum[ichi]<=10 and chargedIso<=0.1): continue
        if not c.chiCands_passPFCandVeto[ichi]: continue												
        if not c.chiCands_trkRelIso[ichi]<0.2: continue
        if not c.chiCands[ichi].Pt()*c.chiCands_trkRelIso[ichi]<10: continue				
        nhits = c.chiCands_nValidTrackerHits[ichi]
        nlayers = c.chiCands_trackerLayersWithMeasurement[ichi]
        if not (nlayers>=2 and nhits>=2): continue
        if not (c.chiCands_nMissingInnerHits[ichi]==0): continue
        pterr = c.chiCands_ptError[ichi]/(chiCand.Pt()*chiCand.Pt())
        if short: 
            if not (pterr<0.2): continue
        if medium: 
            if not (pterr<0.05): continue
        if long: 
            if not (pterr<0.005): continue								
        if not (c.chiCands_trackQualityHighPurity[ichi]): continue

        sumtagvec+=chiCand
        taggedtracks.append(chiCand)

        if short: nshort+=1
        if medium: nmedium+=1
        if long: nlong+=1 		


    var_St[0] = c.HT
    for trk in taggedtracks: var_St[0]+=trk.Pt()

    var_NShortTags[0] = nshort
    var_NMediumTags[0] = nmedium
    var_NLongTags[0] = nlong
    var_NTags[0] = nshort+nmedium+nlong
    if var_NTags[0]>0: 
        var_DPhiMetSumTags[0] = abs(metvec.DeltaPhi(sumtagvec))
        var_SumTagPtOverMet[0] = sumtagvec.Pt()/metvec.Pt()
    else:
        var_DPhiMetSumTags[0] = -10.0
        var_SumTagPtOverMet[0] = -10.0
    var_CrossSection[0] = c.CrossSection
    if var_NTags[0]>0:
    	tEvent.Fill()


fnew.cd()
tEvent.Write()
print 'just created', fnew.GetName()
hHt.Write()
hHtWeighted.Write()
fnew.Close()
