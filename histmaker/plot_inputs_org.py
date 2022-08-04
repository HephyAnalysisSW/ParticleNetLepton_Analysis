#!/usr/bin/env python

import os

import ROOT, os, copy, pickle
ROOT.gROOT.SetBatch(True)
from math import sqrt, cos, sin, pi, cosh, isnan
from RootTools.core.standard import *

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
['PFCand_eta_rel', -0.8, 0.8, 100],
['PFCand_phi_rel', -0.8, 0.8, 100],
['PFCand_phiAtVtx_rel', -3.2, 3.2, 60],
['PFCand_deltaR', 0, 0.6, 200],
['PFCand_caloFraction', 0, 3, 100],
['PFCand_hcalFraction', 0, 1, 90], 
['PFCand_hcalFractionCalib', -1.5, 1.5, 60],
['PFCand_puppiWeight', 0, 1, 50],
['PFCand_puppiWeightNoLep',  0, 1, 50],
['PFCand_dz', -500, 500, 90], 
['PFCand_dzSig', -40000, 40000, 60],
['PFCand_dxy', -130, 130, 90],
['PFCand_dxySig', -4000, 4000, 60],
['PFCand_charge', -1, 2, 3],
['PFCand_pvAssocQuality', 0, 8, 90], 
['PFCand_status', 0, 2200, 60],
['PFCand_pixelhits', 0, 12, 12], 
['PFCand_nTrackerLayers', 0, 22, 22]
]

SV_list=[
['SV_eta_rel', 90, -0.6, 0.6],
['SV_phi_rel', 90, -0.6, 0.6],
['SV_deltaR', 90, 0, 0.6],
['SV_chi2', 60, -2000, 2000],
['SV_ndof', 17, 0, 17],
['SV_pt_rel', 90, 0, 11],
['SV_mass', 60, 0, 35],
['SV_dxy', 90, 0, 35],
['SV_d3d', 60, 0, 140],
['SV_ntracks', 10, 1, 11],
['SV_dxySig', 90, 0, 3000],
['SV_d3dSig', 90, 0, 3500],
['SV_cospAngle', 60, -1, 1]
]


lepton_list=[
['lepton_pt', 90, 0, 800],
['lepton_eta', 90, -3, 3],
['lepton_phi', 90, -3.5, 3.5],
['lepton_mass', 60, -1, 1],
['lepton_charge', 3, -1 , 2],
['lepton_tightcharge', 3, 0, 3],
['lepton_pdgId', 30, -15, 15],
['lepton_dxy', 90, -4, 16]
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

output_file = TFile( output_filename, 'recreate')

#define histograms

PFcand_hist=[[ROOT.TH1F('h_'+var[0]+label, var[0]+'_'+label, var[3], var[1], var[2]) for var in PFCand_list] for label in label_list]
lepton_hist=[[ROOT.TH1F('h_'+var[0]+label, var[0]+'_'+label, var[1], var[2], var[3]) for var in lepton_list] for label in label_list]
SV_hist=[[ROOT.TH1F('h_'+var[0]+label, var[0]+'_'+label, var[1], var[2], var[3]) for var in SV_list] for label in label_list]

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
			for label in range(len(PFcand_hist)):
				if getattr(tree, label_list[label])>0:
					for ihx in range(len(PFcand_hist[label])):
						x=getattr(tree, PFCand_list[ihx][0])
						for i in range(len(x)):
							PFcand_hist[label][ihx].Fill(x[i])
			for label in range(len(lepton_hist)):
                                if getattr(tree, label_list[label])>0:
                                        for ihx in range(len(lepton_hist[label])):
                                                x=getattr(tree, lepton_list[ihx][0])
                                                lepton_hist[label][ihx].Fill(x)
			for label in range(len(SV_hist)):
                                if getattr(tree, label_list[label])>0:
                                        for ihx in range(len(SV_hist[label])):
                                                x=getattr(tree, SV_list[ihx][0])
                                                for i in range(len(x)):
                                                        SV_hist[label][ihx].Fill(x[i])
		#print ('Entry: ', ievt)
	y.Close()

ROOT.gROOT.SetBatch(True)

output_file.cd()
output_file.Write()
output_file.Close()
logger.info( "Written output file %s", output_filename )
