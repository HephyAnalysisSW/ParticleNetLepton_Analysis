# ParticleNetLepton_Analysis

## Follow the instructions below to set up the framework:

export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_6_0

cd CMSSW_10_6_0/src

cmsenv

git cms-init

git clone https://github.com/HephyAnalysisSW/ParticleNetLepton_Analysis.git

./ParticleNetLepton_Analysis/setup_106X.sh

For a test run:

#Set up voms proxy

cd ParticleNetLepton_Analysis/histmaker

python plot_inputs.py --sample TT_Semilep --small
