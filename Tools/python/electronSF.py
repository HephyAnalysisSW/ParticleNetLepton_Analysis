import ROOT
import os
from math import sqrt

from VH.Tools.correctionlib import _core

LepID_keys       = { 
                      "Medium":"Medium",
		      "MVA_WP80NoIso":"wp80noiso",
		      "MVA_WP90NoIso":"wp90noiso",
		      "MVA_WP80Iso":"wp80iso",
		      "MVA_WP90Iso":"wp90iso",
                    }

class electronSF:

    def __init__(self, year="2016postVFP_UL", ID=None):

        if year not in [ "2016preVFP_UL","2016postVFP_UL","2017_UL","2018_UL" ]:
            raise Exception("Lepton SF for year %i not known"%year)

        self.dataDir = "$CMSSW_BASE/src/jsonpog-integration/POG/EGM/"+year+"/"
        self.year    = year

        if not ID in LepID_keys.keys():
            raise Exception("Don't know ID %s"%ID)

	filename = os.path.expandvars(os.path.join(self.dataDir, "electron.json.gz"))
	#self.dataDir+"electron.json.gz"
	evaluator = _core.CorrectionSet.from_file(filename)

	self.key = LepID_keys[ID]
	self.lepSF = evaluator["UL-Electron-ID-SF"]
    
    def getPartialSF( self, year, pt, eta, var):
	yearnum = year.split("_")[0]
	sf = self.lepSF.evaluate(yearnum, var, self.key, eta, pt)
	return sf
   
    def mult( self, list ):
        res = list[0]
        for i in list[1:]: res = res*i
        return res

    def getSF(self, pdgId, pt, eta, unc="nominal"):

	sf = float(1)
	pt, eta = float(pt), float(eta)

        if abs(pdgId)!=11:
            raise Exception("Lepton SF for PdgId %i not known"%pdgId)

        if not unc in ["nominal", "systup", "systdown"]:
            raise Exception("Don't know uncertainty %s"%unc)

	var = "sf"
        if unc == "nominal":
            var = "sf"
        elif unc == "systup":
            var = "sfup"
       	elif unc == "systdown":
	    var = "sfdown"

	sf = self.getPartialSF( self.year, pt, eta, var) 

        return sf


if __name__ == "__main__":

    print "2016postVFP_UL, Medium"
    LSF = electronSF(year="2016postVFP_UL", ID="MVA_WP90NoIso")
    print LSF.getSF(11, 200, 1.0)
    print LSF.getSF(11, 200, 1.0, unc="systup")
    print LSF.getSF(11, 200, 1.0, unc="systdown")
