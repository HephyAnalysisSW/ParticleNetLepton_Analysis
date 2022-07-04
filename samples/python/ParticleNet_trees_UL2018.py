import copy, os, sys
from RootTools.core.Sample import Sample
import ROOT

import logging
logger = logging.getLogger(__name__)

try:
	directory_ = sys.modules['__main__'].directory_
except:
	directory_ = "/eos/vbc/group/cms/suman.chatterjee/ParticleNetLepton_Files/MC/"

logger.info("Loading MC samples from directory %s", directory_)

def make_dirs( dirs ):
	return [ os.path.join( directory_, dir_ ) for dir_ in dirs ]

dirs = {}

dirs['TT_Semilep']         = ["TT_Semilep"]
TT_Semilep = Sample.fromDirectory(name="TT_Semilep", redirector = None, treeName="Events", isData=False, color=ROOT.kBlack, texName="ZH", directory=make_dirs( dirs['TT_Semilep']),xSection=1)
TT_Semilep.normalization = 1
