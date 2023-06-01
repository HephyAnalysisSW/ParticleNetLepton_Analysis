#!/usr/bin/env python
''' script for producing histograms and flat ntuples 
'''
#
# Standard imports and batch mode
#
import os

import ROOT, os, copy, pickle
ROOT.gROOT.SetBatch(True)
from math import sqrt, cos, sin, pi, cosh, isnan
from RootTools.core.standard import *

from ROOT import TLorentzVector, TVector3, TRandom3, TH1D, TH2D, TFile

from Analysis.Tools.HyperPoly   import HyperPoly
from Analysis.Tools.metFilters  import getFilterCut
from Analysis.Tools.BTagEfficiency  import BTagEfficiency

from Histo_EFT import *

import operator                                                                                                                                                                                             
from operator import mul
import array as arr
import importlib

import itertools
from itertools import combinations_with_replacement

from tqdm import tqdm

from Analysis.Tools.WeightInfo                   import WeightInfo
#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO',      nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',       action='store_true',                                            help="Run the file on a small sample (for test purpose), bool flag set to True if used" )
argParser.add_argument('--noData',         action='store_true', default=False, help='also plot data?')
argParser.add_argument('--sample',             action='store',      default='TT_Dilep_0j_SMEFTatNLO_LO_MG3X',   type=str,      help="Which sample?")
argParser.add_argument('--nJobs',              action='store',      nargs='?', type=int, default=1,     help="Maximum number of simultaneous jobs.")
argParser.add_argument('--job',                action='store',      nargs='?', type=int, default=0,     help="Run only job i")
argParser.add_argument('--maxN',               action='store',      default=1000,   type=int,          help="Maximum number of events (for small run)?")

args = argParser.parse_args()

#
# Logger
#
import VH.Tools.logger as logger
#import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
#logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

#Samples
exec("from ParticleNetLepton_Analysis.samples.ParticleNet_trees_UL2018 import  %s"%(args.sample))
#
#
data_sample = eval(args.sample)
org_sample = eval(args.sample)

if args.nJobs>1:
        n_files_before = len(data_sample.files)
        data_sample = data_sample.split(args.nJobs)[args.job]
        n_files_after  = len(data_sample.files)
        logger.info( "Running job %i/%i over %i files from a total of %i.", args.job, args.nJobs, n_files_after, n_files_before)

# output directory
from ParticleNetLepton_Analysis.Tools.user import histogram_output_directory
skim_output_directory               =  histogram_output_directory+'/Normalization' 
output_directory = os.path.join(skim_output_directory, org_sample.name)

# make small
if args.small:
    data_sample.reduceFiles( to = 1 )

# define binnings (for plots) here:

# preselection condition:
#
preSelection = "1>0"#Generator_weight>-1000"
#
# Read variables: 
#
#read MC variables: 
read_variables_MC = ["Generator_weight/F"]

#list of variables to be stored in flat ntuple

all_read_variables = []
all_read_variables.extend(read_variables_MC)

data_sample.read_variables = all_read_variables
data_sample.texName = data_sample.texName + " (SM)"
data_sample.setSelectionString(preSelection)
data_sample.weight = lambda event, sample: event.weight if sample.isData else 1#event.reweightPU*event.reweightL1Prefire*lumi_scale 
data_sample.style = styles.lineStyle( ROOT.kBlack, width=2)

dataMCScale = 1

tmp_dir     = ROOT.gDirectory

output_filename =  os.path.join(output_directory, org_sample.name + str(args.job) + '.root')

if not os.path.exists( output_directory ): 
    try:
        os.makedirs( output_directory )
    except OSError:
        pass
    logger.info( "Created output directory %s", output_directory )

output_file = TFile( output_filename, 'recreate')
#output_file.cd()

# define histos ##
print "Defining histos"

# histogram array (1D) for LHE-level info #

hist_norm = ROOT.TH1D("h_norm","",1,0.5,1.5)
hist_count = ROOT.TH1D("h_count","",1,0.5,1.5)

filenames=data_sample.files

files = [ROOT.TFile.Open(x) for x in filenames]
fl=filter(None, files)

TreeName = "deepntuplizer"

for y in fl:
	dirc = y.Get(TreeName)
	tree=dirc.Get("tree")
	for ievt in range(tree.GetEntries()):
		tree.GetEntry(ievt)
		hist_norm.Fill(1,getattr(tree, 'Generator_weight'))
		hist_count.Fill(1,1)

print "end of it"

output_file.cd()
output_file.Write()
output_file.Close()

#del output_file

logger.info("Written output file %s", output_filename)
