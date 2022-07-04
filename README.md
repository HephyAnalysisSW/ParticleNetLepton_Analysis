# ParticleNetLepton_Analysis

cmsrel CMSSSW_10_6_0

cd CMSSW_10_6_0/src

cmsenv

git cms-init

git clone https://github.com/HephyAnalysisSW/ParticleNetLepton_Analysis.git

./setup_106X.sh

For a test run:

#Set up voms proxy

cd ParticleNetLepton_Analysis/histmaker

python plot_inputs.py --sample TT_Semilep --small
