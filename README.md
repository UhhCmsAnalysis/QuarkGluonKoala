# Uhhhh
This is a package for analyzing events with HH-like final states on the NAF.
## Set up code

```
cmsrel CMSSW_10_1_0
cd CMSSW_10_1_0/src
cmsenv
git clone https://github.com/UhhCmsAnalysis/QuarkGluonKoala.git/
cd QuarkGluonKoala.git/
mkdir output output/smallchunks output/mediumchunks output/bigchunks
mkdir jobs
mkdir pdfs
```

## run event analyzer script to generate histograms
### examples:

run over a semi-leptonic tt-jets file:
```
python tools/EventAnalyzer.py --fnamekeyword "/pnfs/desy.de/cms/tier2/store/user/vormwald/NtupleHub/ProductionRun2v3/Summer16.TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8AOD_40*_RA2AnalysisTree.root"
```

run over inclusive tt-jets file:
```
python tools/EventAnalyzer.py --fnamekeyword  "/pnfs/desy.de/cms/tier2/store/user/vormwald/NtupleHub/ProductionRun2v3/Summer16.TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8AOD_10*_RA2AnalysisTree.root"
```

run over w+jets:
```
python tools/EventAnalyzer.py --fnamekeyword  "/pnfs/desy.de/cms/tier2/store/user/vormwald/NtupleHub/ProductionRun2v3/Summer16.WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8AOD_0_RA2AnalysisTree.root"
```

run over signal:

```
python tools/EventAnalyzer.py --fnamekeyword  "/nfs/dust/cms/user/nissanuv/x1x2x1/signal/ntuples_sum/higgsino_mu100_dm3p28Chi20Chipm_2.root"
```

run over z+jets:
```
python tools/EventAnalyzer.py --fnamekeyword  "/pnfs/desy.de/cms/tier2/store/user/vormwald/NtupleHub/ProductionRun2v3/Summer16.ZJetsToNuNu_HT-100To200_13TeV-madgraphAOD_21*_RA2AnalysisTree.root"
```

## parallelize running of scripts with condor jobs

before submitting jobs, it's good to clean out the jobs/ directory where the output/error files are cached:

```
bash tools/CleanBird.sh
```

### examples:

submit semi-leptonic tt+jets (one job per input file)

```
python tools/SubmitJobs_condor.py --analyzer tools/EventAnalyzer.py --fnamekeyword  "/pnfs/desy.de/cms/tier2/store/user/vormwald/NtupleHub/ProductionRun2v3/Summer16.TTJets_SingleLept*.root"
```

submit QCD (one job per input file)

```
python tools/SubmitJobs_condor.py --analyzer tools/EventAnalyzer.py --fnamekeyword  "/pnfs/desy.de/cms/tier2/store/user/vormwald/NtupleHub/ProductionRun2v3/Summer16.QCD_HT*.root"
```
 
files will show up in the folder output/smallchunks. You can check the status of jobs and interact with jobs using
```
condor_q
condor_rm <job number>
condor_release <username>
```

## combine the files 
when the jobs have finished, you can combine the many files into just a few, while applying the denominator of the event weights (nevents):
```
python tools/MergeHistsPropWeights.py output/smallchunks/
```

This produces a set of files with histograms weighted to lumi=1/pb in output/bigchunks

## draw stacked histograms of background overlay signal:

```
python tools/EventPlotter.py
```

This produces a file, probably plots.root, which you can open and examine the canvases within. 


