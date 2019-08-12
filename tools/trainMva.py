#!/usr/bin/env python

#----------------------------------------------------------------------
# File: train.py
# Description: regression with TMVA
# Created: 04-August-2014, Sandwich Convention 2014, Chicago, IL
#----------------------------------------------------------------------
import os, sys, re
from ROOT import *
from math import *
from string import *
from time import sleep
from array import array

gROOT.SetBatch(1)

def formatOptions(options):
    from string import joinfields, strip, split
    options = joinfields(map(strip, split(strip(options), '\n')), ':')
    return options
#----------------------------------------------------------------------
#----------------------------------------------------------------------
def main():
    print "\n", "="*80
    print "\tclassification with TMVA"
    print "="*80

    
    try:
        infileS, infileB = sys.argv[1], sys.argv[2]
    except:
        print "please specify signal and bkg input files as first and second arguments"
        exit(0)

    FileS = TFile(infileS)
    FileB = TFile(infileB)    
    SigTree = FileS.Get('tEvent')
    BkgTree = FileB.Get('tEvent')

    if 'TagProbe' in infileS:
        cuts = ''
        SigTree.Draw('mjj',cuts)
        hcountS = SigTree.GetHistogram().Clone('hcountS')
        BkgTree.Draw('mjj',cuts)
        hcountB = BkgTree.GetHistogram().Clone('hcountB')
    if 'IdealSample' in infileS:
        cuts = 'trueResponse>0.6 && trueResponse<1.4 '
        SigTree.Draw('mjj',cuts)
        hcountS = SigTree.GetHistogram().Clone('hcountS')
        BkgTree.Draw('mjj',cuts)
        hcountB = BkgTree.GetHistogram().Clone('hcountB')
            
    nevents =min(min(hcountS.GetEntries(),hcountB.GetEntries()), 3000000)-4#/2

    methods = 'Mlp'
    methods = 'Bdt'
    #methods = 'Dnn'
    
    outFileName = 'Basic'+methods+'_'+FileS.GetName()[FileS.GetName().rfind('/')+1:]
    print 'creating file', outFileName
    outputFile = TFile(outFileName, "RECREATE")
    
    factory = TMVA.Factory("TMVAClassification", outputFile,
                            ":".join([
                                "!V",
                                "!Silent",
                                "Color",
                                "DrawProgressBar",
                                "!V:Transformations=I;D;P;G,D",
                                "AnalysisType=Classification"]
                                     ))
    
    #dataloader=TMVA.DataLoader ( "dataset" )
    factory.AddVariable("pt", 'F')
    factory.AddVariable("eta", 'F')
    factory.AddVariable("qgMult", 'F')
    factory.AddVariable("qgPtD", 'F')
    factory.AddVariable("qgAxis2", 'F')    

    factory.AddSignalTree(SigTree)
    factory.AddBackgroundTree(BkgTree)

    factory.PrepareTrainingAndTestTree(TCut(cuts),   # signal events
                                   TCut(cuts),    # background events
                                   ":".join([
                                        "nTrain_Signal="+str(nevents/2),
                                        "nTrain_Background="+str(nevents/2),
                                        "nTest_Signal="+str(nevents/2),
                                        "nTest_Background="+str(nevents/2),
                                        "SplitMode=Random",
                                        "NormMode=NumEvents",
                                        "!V"
                                       ]))
    
    #    factory.BookMethod( TMVA.Types.kMLP, "MLP_orig", 
    #    "!H:!V:NeuronType=tanh:NCycles=100:HiddenLayers=N:TestRate=5:BPMode=batch:BatchSize=-1" 
    # );

    if 'Bnn' in methods:
        factory.BookMethod( TMVA.Types.kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" )


    if 'Mlp' in methods:#xxx
        factory.BookMethod( TMVA.Types.kMLP,
                        "MLP",
                        "!H:!V:"\
                        "VarTransform=N:"\
                        "LearningRate=0.8:"\
                        "HiddenLayers=19:"\
                        "TrainingMethod=BFGS:"\
                        "UseRegulator=false")
        #factory.BookMethod( TMVA.Types.kMLP,"MLP","H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:!UseRegulator")
        

    if 'Bdt' in methods:
        factory.BookMethod(  TMVA.Types.kBDT,
                        "BDT",
                        "!V:"\
                        "BoostType=AdaBoost:"\
                        "NTrees=200:"\
                        "nEventsMin=100:"\
                        "nCuts=50")    

    if 'Dnn' in methods:
        layoutString = "Layout=TANH|128,TANH|128,TANH|128,LINEAR"
        training0 = "LearningRate=1e-1,Momentum=0.9,Repetitions=1, ConvergenceSteps=20,BatchSize=256,TestRepetitions=10, WeightDecay=1e-4,Regularization=L2, DropConfig=0.0+0.5+0.5+0.5, Multithreading=True"
        training1 = "LearningRate=1e-2,Momentum=0.9,Repetitions=1, ConvergenceSteps=20,BatchSize=256,TestRepetitions=10, WeightDecay=1e-4,Regularization=L2, DropConfig=0.0+0.0+0.0+0.0, Multithreading=True"
        training2 = "LearningRate=1e-3,Momentum=0.0,Repetitions=1, ConvergenceSteps=20,BatchSize=256,TestRepetitions=10, WeightDecay=1e-4,Regularization=L2, DropConfig=0.0+0.0+0.0+0.0, Multithreading=True"
        trainingStrategyString = "TrainingStrategy=" + training0 + "|" + training1 + "|" + training2
        stdOptions = "!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:WeightInitialization=XAVIERUNIFORM:Architecture=STANDARD"+':'+trainingStrategyString+':'+layoutString
        factory.BookMethod( TMVA.Types.kDNN, "DNN", stdOptions);

    #harry                        
  
    factory.TrainAllMethods()  
    factory.TestAllMethods()
    factory.EvaluateAllMethods()
    
    outputFile.Close()

    del factory
    #os.system('mv weights weightfiles/'+outFileName.replace('.root','_weights').replace('../treesoflife/',''))
#----------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "\nciao"
