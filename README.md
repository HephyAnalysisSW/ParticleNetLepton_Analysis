# ParticleNetLepton_Analysis

## Follow the instructions below to set up the framework:

- export SCRAM_ARCH=slc7_amd64_gcc700

- cmsrel CMSSW_10_6_0

- cd CMSSW_10_6_0/src

- cmsenv

- git clone https://github.com/HephyAnalysisSW/ParticleNetLepton_Analysis.git

- ./ParticleNetLepton_Analysis/setup_106X.sh

## For a test run:

- voms-proxy-init -rfc -voms cms -valid 48:00

- cd ParticleNetLepton_Analysis/histmaker

- cmsenv

- python plot_inputs.py --sample TT_Semilep --small

## For submitting batch job:

- voms-proxy-init -rfc -voms cms -valid 48:00

- cd ParticleNetLepton_Analysis/histmaker

- cmsenv

- ./submit_jobs.sh
