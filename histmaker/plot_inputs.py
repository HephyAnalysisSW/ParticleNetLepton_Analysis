#!/usr/bin/env python

import os

import ROOT, os, copy, pickle
ROOT.gROOT.SetBatch(True)
from math import sqrt, cos, sin, pi, cosh, isnan
from RootTools.core.standard import *
import array as arr

from ROOT import TLorentzVector, TVector3, TRandom3, TH1D, TH2D, TFile

#arguments 

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--small',       action='store_true',                                            help="Run the file on a small sample (for test purpose), bool flag set to True if used" )
argParser.add_argument('--sample',             action='store',      default='ZH_fast',   type=str,      help="Which sample?")
argParser.add_argument('--nJobs',              action='store',      nargs='?', type=int, default=1,     help="Maximum number of simultaneous jobs.")
argParser.add_argument('--job',                action='store',      nargs='?', type=int, default=0,     help="Run only job i")
argParser.add_argument('--maxN',               action='store',      default=10000,   type=int,          help="Maximum number of events (for small run)?")
argParser.add_argument('--logLevel',       action='store',      default='INFO',      nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")

args = argParser.parse_args()

def getbinid(val,array):
	if (val<array[0]): 
		return -2;
	for ix in range(len(array)):
		if (val < array[ix]):
			return ix-1;                                                                                                                                                                        
	return -3;


#logging

import logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

#load sample

exec("from ParticleNetLepton_Analysis.samples.ParticleNet_trees_UL2018 import  %s"%(args.sample))
data_sample = eval(args.sample)
org_sample = eval(args.sample)

#split into jobs
if args.nJobs>1:
    n_files_before = len(data_sample.files)
    data_sample = data_sample.split(args.nJobs)[args.job]
    n_files_after  = len(data_sample.files)
    logger.info( "Running job %i/%i over %i files from a total of %i.", args.job, args.nJobs, n_files_after, n_files_before)

# output directory
from ParticleNetLepton_Analysis.Tools.user import histogram_output_directory
skim_output_directory               =  histogram_output_directory 
output_directory = os.path.join(skim_output_directory, org_sample.name)

if args.small:
    data_sample.reduceFiles( to = 1 )

#declaring variables for histograms

PFCand_list=[
['PFCand_eta_rel', 100, -0.5, 0.5],
['PFCand_phi_rel',100, -0.5, 0.5],
['PFCand_phiAtVtx_rel', 100, -0.5, 0.5],
['PFCand_deltaR', 200, 0, 0.5],
['PFCand_caloFraction', 100, 0, 2],
['PFCand_hcalFraction', 90, 0, 1], 
['PFCand_hcalFractionCalib', 50, 0, 1],
['PFCand_puppiWeight', 100, 0, 1],
['PFCand_puppiWeightNoLep', 100, 0, 1],
['PFCand_dz', 100, -10, 10], 
['PFCand_dzSig', 100, -500, 500],
['PFCand_dxy', 100, -2, 2],
['PFCand_dxySig', 100, -50, 50],
['PFCand_charge', 3, -1, 2],
['PFCand_pvAssocQuality', 8, 0, 8], 
['PFCand_status', 60, 0, 2200],
['PFCand_pixelhits', 12, 0, 12], 
['PFCand_nTrackerLayers', 22, 0, 22]
]

SV_list=[
['SV_eta_rel', 100, -0.5, 0.5],
['SV_phi_rel', 100, -0.5, 0.5],
['SV_deltaR', 100, 0, 0.5],
['SV_chi2', 100, 0, 20],
['SV_ndof', 17, 0, 17],
['SV_pt_rel', 100, 0, 10],
['SV_mass', 100, 0, 25],
['SV_dxy', 100, 0, 5],
['SV_d3d', 100, 0, 20],
['SV_ntracks', 10, 0, 10],
['SV_dxySig', 90, 0, 3000],
['SV_d3dSig', 100, 0, 250],
['SV_cospAngle', 200, -1, 1]
]


lepton_list=[
['lepton_pt', 50, 0, 500],
['lepton_eta', 90, -2.5, 2.5],
['lepton_phi', 90, -3.5, 3.5],
['lepton_mass', 50, 0, 0.5],
['lepton_charge', 3, -1.5 , 1.5],
['lepton_tightcharge', 3, 0, 3],
['lepton_pdgId', 30, -15, 15],
['lepton_dxy', 60, -4, 4],
['lepton_dz', 60, -50, 50],
['lepton_dxyError', 50, 0, 0.1],
['lepton_dzError', 30, 0, 5],
['lepton_dxySig', 75, -500, 500],
['lepton_dzSig', 50, -2000, 2000],
['lepton_ip3d', 60, -25, 25],
['lepton_sip3d', 80, -1000, 1000],
['lepton_dxy_sv', 40, -10, 10],
['lepton_genPartFlav', 60, 0, 25],
['lepton_chi2', 100, 0, 2000],
['lepton_ndof', 100, 0, 100],
['lepton_hit', 60, 0, 60],
['lepton_pixhit', 10, 0, 10],
['lepton_nTrackerLayers', 20, 0, 20],
['lepton_lostHits', 5, 0, 5],
['lepton_e_ECAL', 100, 0, 1],
['lepton_e_HCAL', 100, 0, 1],
['lepton_hoe', 90, 0, 150],
['lepton_minisoch', 100, 0, 5],
['lepton_minisonh', 100, 0, 2],
['lepton_minisoph', 100, 0, 5],
['lepton_minisoall', 100, 0, 10],
['lepton_pfRelIso03_drcor', 100, 0, 10],
['lepton_pfRelIso03_ChargedHadron', 50, 0, 5],
['lepton_pfRelIso03_NeutralHadron', 50, 0, 2],
['lepton_pfRelIso03_Photon', 50, 0, 5],
['lepton_pfRelIso03_PileUp', 50, 0, 2],
['lepton_tkRelIso', 50, 0, 5],
['lepton_trkKink', 50, 0, 50],
['lepton_nStations', 6, 0, 6],
['lepton_segmentComp', 60, 0, 1],
['lepton_posmatch', 100, 0, 5000],
['lepton_isPFCand', 2, -0.5, 1.5],
['lepton_isGlobal', 2, -0.5, 1.5],
['lepton_isTracker', 2, -0.5, 1.5],
['lepton_isLoose', 2, -0.5, 1.5],
['lepton_isGoodGlobal', 2, -0.5, 1.5],
['lepton_isMedium', 2, -0.5, 1.5],
['lepton_isMedPr', 2, -0.5, 1.5],
['lepton_isTight', 2, -0.5, 1.5],
['lepton_isHighPt', 2, -0.5, 1.5],
['lepton_isHighPttrk', 2, -0.5, 1.5],
['lepton_jetPtRelv2', 50, 0, 0.5],
['lepton_jetPtRelv2_log', 50, -5, 0],
['lepton_jetRelIso', 90, 0, 10],
['lepton_jetbtag', 50, 0, 1] 
] 

#'PFCand_trackHighPurity', 'PFCand_isElectron', 'PFCand_isMuon', 'PFCand_isChargedHadron', 'PFCand_fromPV']

label_list=['label_Muon_Prompt', 'label_Muon_fromTau', 'label_Muon_fromHFHadron', 'label_Muon_fromLFHadron', 'label_Muon_unknown', 'label_Muon_fromPhoton']

#create output file

output_filename =  os.path.join(output_directory, org_sample.name + str(args.job) + '.root')
if not os.path.exists( output_directory ):
    try:
        os.makedirs( output_directory )
    except OSError:
        pass
    logger.info( "Created output directory %s", output_directory )

ref_ptbins = arr.array('d', [20,50,100,10000])
ref_etabins = arr.array('d', [0,1.4,2.4])

output_file = TFile( output_filename, 'recreate')

#define histograms

PFcand_hist=[[[[ROOT.TH1F('h_'+var[0]+label+'_ptbin'+str(ipt+1)+'_etabin'+str(ieta+1), var[0]+'_'+label+'_ptbin'+str(ipt+1)+'_etabin'+str(ieta+1), var[1], var[2], var[3]) for var in PFCand_list] for ieta in range(len(ref_etabins)-1)] for ipt in range(len(ref_ptbins)-1)] for label in label_list]
lepton_hist=[[[[ROOT.TH1F('h_'+var[0]+label+'_ptbin'+str(ipt+1)+'_etabin'+str(ieta+1), var[0]+'_'+label+'_ptbin'+str(ipt+1)+'_etabin'+str(ieta+1), var[1], var[2], var[3]) for var in lepton_list] for ieta in range(len(ref_etabins)-1)] for ipt in range(len(ref_ptbins)-1)] for label in label_list]
SV_hist=[[[[ROOT.TH1F('h_'+var[0]+label+'_ptbin'+str(ipt+1)+'_etabin'+str(ieta+1), var[0]+'_'+label+'_ptbin'+str(ipt+1)+'_etabin'+str(ieta+1), var[1], var[2], var[3]) for var in SV_list] for ieta in range(len(ref_etabins)-1)] for ipt in range(len(ref_ptbins)-1)] for label in label_list]

#load files
#filenames=['/eos/vbc/experiments/cms/store/user/agruber/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/ParticleNetLepton_UL2018_TT_Semilep/220527_142442/0000/tree_'+str(256)+'.root'] # for i in range (1000)] 
filenames=data_sample.files
print filenames

files = [ROOT.TFile.Open(x) for x in filenames]
fl=filter(None, files)
print "# of files: ",len(fl)

TreeName = "deepntuplizer"
pt_cut = 20
eta_cut = 2.4

for y in fl:
    dirc = y.Get(TreeName)
    tree=dirc.Get("tree")
    for ievt in range(tree.GetEntries()):
        if args.small:
            if ievt>args.maxN:
                break
        tree.GetEntry(ievt)
        if tree.lepton_pt>pt_cut and abs(tree.lepton_eta)<eta_cut:
	    ipt = getbinid(tree.lepton_pt,ref_ptbins)
	    ieta = getbinid(abs(tree.lepton_eta),ref_etabins)
            for label in range(len(PFcand_hist)):
                if getattr(tree, label_list[label])>0:
                    for ihx in range(len(PFcand_hist[label])):
                        x=getattr(tree, PFCand_list[ihx][0])
                        weight = getattr(tree, 'Generator_weight')
                        for i in range(len(x)):
                            PFcand_hist[label][ipt][ieta][ihx].Fill(x[i], weight)
            for label in range(len(lepton_hist)):
                if getattr(tree, label_list[label])>0:
                    for ihx in range(len(lepton_hist[label])):
                        x=getattr(tree, lepton_list[ihx][0])
                        weight = getattr(tree, 'Generator_weight')
                        lepton_hist[label][ipt][ieta][ihx].Fill(x, weight)
            for label in range(len(SV_hist)):
                if getattr(tree, label_list[label])>0:
                    for ihx in range(len(SV_hist[label])):
                        x=getattr(tree, SV_list[ihx][0])
                        weight = getattr(tree, 'Generator_weight')
                        for i in range(len(x)):
                            SV_hist[label][ipt][ieta][ihx].Fill(x[i], weight)
        #print ('Entry: ', ievt)
    y.Close()

ROOT.gROOT.SetBatch(True)

output_file.cd()
output_file.Write()
output_file.Close()
logger.info( "Written output file %s", output_filename )
