from ROOT import *
from utils import *
gStyle.SetOptStat(0)
import os, sys
datamc = 'MC'

try:
    fnameQ = sys.argv[1]
    fnameG = sys.argv[2]    
except:
    print 'doing this shit'
    fnameQ = 'TagProbeTreeQ.root'
    fnameG = 'TagProbeTreeGclever.root'
    #fg = TFile('TagProbeTreeG.root')

fq = TFile(fnameQ)
fg = TFile(fnameG)
tq = fq.Get('tEvent')
tg = fg.Get('tEvent')

tq.Show(0)
    
print 'n(q)=', tq.GetEntries(), 'n(g)=', tg.GetEntries()
branches = tq.GetListOfBranches()

c1 = mkcanvas_wide('c1')
pad1 = c1.GetPad(1)
#pad1.SetLogy()
pad2 = c1.GetPad(2)

hrocbackground = TH2F('','',20,0,1,20,0,1)
hrocbackground.GetXaxis().SetTitle("#epsilon_{s}")
hrocbackground.GetYaxis().SetTitle("#epsilon_{b}")

for branch in branches:
    
    name = branch.GetName()
    if name in ['ht']: continue
    #if not name=='qgLikelihood': continue
    print 'processing', name
    pad1.cd()
    constraint = 'pt>150 && mjj<1000'
    constraint = ''
    tg.Draw(name)
    hg = tg.GetHistogram().Clone('hg_'+name)    
    histoStyler(hg,kRed+1)            
    hg.SetFillColor(hg.GetLineColor())    
    hg.SetFillStyle(3004)            
    tg.Draw(name, constraint+' && partonFlavor>6')
    hgT = tg.GetHistogram().Clone('hgT_'+name)    
    histoStyler(hgT,kRed-4)            
    hgT.SetFillColor(hgT.GetLineColor())    
    hgT.SetFillStyle(3244)            
        
    tq.Draw(name,'','same')
    hq = tq.GetHistogram().Clone('hq_'+name)    
    scale = 1.0#*tg.GetEntries()/tq.GetEntries()
    hq.Scale(scale)
    histoStyler(hq,kAzure-1)    
    hq.SetFillColor(hq.GetLineColor())    
    hq.SetFillStyle(1001)    
    tq.Draw(name, constraint+' && abs(partonFlavor)<6','same')
    hqT = tq.GetHistogram().Clone('hqT_'+name)    
    hqT.Scale(scale)
    histoStyler(hqT,kAzure-2)    
    hqT.SetFillColor(hqT.GetLineColor())    
    hqT.SetFillStyle(1001)    

    hq.SetTitle('')
    hq.GetYaxis().SetRangeUser(0.1,1.7*hq.GetMaximum())    
    hq.Draw('hist')
    hqT.Draw('hist same')    
    hg.Draw('hist same')        
    hgT.Draw('hist same')    
    leg = mklegend(x1=.21, y1=.68, x2=.95, y2=.9)
    leg.AddEntry(hq,'quark enriched (pure=%.2f)'%(1.0*hqT.GetEntries()/hq.GetEntries()),'lf')
    leg.AddEntry(hqT,'quark enriched (matched q)','lf')    
    leg.AddEntry(hg,'gluon enriched (pure=%.2f)'%(1.0*hgT.GetEntries()/hg.GetEntries()),'lf')    
    leg.AddEntry(hgT,'gluon enriched (matched g)','lf')        
    leg.Draw()
    pad1.Update()

    pad2.cd()
    roc = mkroc(name, hg, hq, kGreen+1)
    hrocbackground.Draw()
    roc.Draw('same')
    pad2.Update()
    pause()
    
    
