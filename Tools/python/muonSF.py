import ROOT
import os
from math import sqrt

from VH.Tools.correctionlib import _core

LepID_keys       = { 
                      "Medium":"NUM_MediumID_DEN_TrackerMuons",
		      "Loose":"NUM_LooseID_DEN_TrackerMuons",
		      "Tight":"NUM_TightID_DEN_TrackerMuons", 
                    }

class muonSF:

    def __init__(self, year="2016postVFP_UL", ID=None):

        if year not in [ "2016preVFP_UL","2016postVFP_UL","2017_UL","2018_UL" ]:
            raise Exception("Lepton SF for year %i not known"%year)

        self.dataDir = "$CMSSW_BASE/src/jsonpog-integration/POG/MUO/"+year+"/"
        self.year    = year

        if not ID in LepID_keys.keys():
            raise Exception("Don't know ID %s"%ID)

	#os.path.expandvars(os.path.join(self.dataDir, file)
	filename = os.path.expandvars(os.path.join(self.dataDir, "muon_Z.json.gz")) 

	if filename.endswith(".json.gz"):
		import gzip
		with gzip.open(filename,'rt') as file:
			data = file.read().strip()
			evaluator = _core.CorrectionSet.from_string(data)
	else:							
		evaluator = _core.CorrectionSet.from_file(filename)

	key = LepID_keys[ID]
	self.lepSF = evaluator[key]
    
    def getPartialSF( self, year, pt, eta, var):
	sf = self.lepSF.evaluate(year, abs(eta), pt, var)
	return sf
   
    def mult( self, list ):
        res = list[0]
        for i in list[1:]: res = res*i
        return res

    def getSF(self, pdgId, pt, eta, unc="nominal"):

	sf = float(1)
	pt, eta = float(pt), float(eta)

        if abs(pdgId)!=13:
            raise Exception("Lepton SF for PdgId %i not known"%pdgId)

        if not unc in ["nominal", "systup", "systdown"]:
            raise Exception("Don't know uncertainty %s"%unc)

	var = "sf"
        if unc == "nominal":
            var = "sf"
        elif unc == "systup":
            var = "systup"
       	elif unc == "systdown":
	    var = "systdown"

	sf = self.getPartialSF( self.year, pt, eta, var) 

        return sf


if __name__ == "__main__":

    print "2016postVFP_UL, Medium"
    LSF = muonSF(year="2016postVFP_UL", ID="Medium")
    print LSF.getSF(13, 200, 1.0)
    print LSF.getSF(13, 200, 1.0, unc="systup")
