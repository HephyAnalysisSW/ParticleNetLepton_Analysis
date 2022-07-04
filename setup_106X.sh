#
eval `scram runtime -sh`
cd $CMSSW_BASE/src

#RootTools (for plotting, sample handling, processing)
git clone https://github.com/HephyAnalysisSW/RootTools.git
cd $CMSSW_BASE/src

# Shared samples (miniAOD/nanoAOD)
git clone https://github.com/HephyAnalysisSW/Samples.git
cd $CMSSW_BASE/src

# Shared analysis tools and data
git clone https://github.com/HephyAnalysisSW/Analysis.git
cd $CMSSW_BASE/src

#compile
eval `scram runtime -sh`
cd $CMSSW_BASE/src 
scram b -j10
