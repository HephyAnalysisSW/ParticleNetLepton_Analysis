''' Implementation of b-tagging reweighting
'''

# Standard imports
import ROOT, pickle, itertools, os
from operator import mul

from VH.Tools.correctionlib import _core

# Logging
import logging
logger = logging.getLogger(__name__)

#binning in pt and eta
ptBorders = [20, 30, 50, 70, 100, 140, 200, 300, 600, 1000]
ptBins    = [ [ptBorders[i], ptBorders[i+1]] for i in range(len(ptBorders)-1) ]
ptBins   += [ [ptBorders[-1], -1] ]

etaBins2016 = [[0,2.4]]
etaBins2017 = [[0,2.5]]
etaBins2018 = [[0,2.5]]

def toFlavourKey(pdgId):
    if abs(pdgId)==5: return ROOT.BTagEntry.FLAV_B
    if abs(pdgId)==4: return ROOT.BTagEntry.FLAV_C
    return ROOT.BTagEntry.FLAV_UDSG

#Method 1ab

effFile2016DeepCSV = 'TTLep_pow_2016_2j_2l_DeepB_eta.pkl'
effFile2017DeepCSV = 'TTLep_pow_2017_2j_2l_DeepB_eta.pkl'
effFile2018DeepCSV = 'TTLep_pow_2018_2j_2l_DeepB_eta.pkl'

effFile2016DeepJet = 'TTLep_pow_2016_2j_2l_DeepFlavB_eta_v2.pkl'
effFile2017DeepJet = 'TTLep_pow_2017_2j_2l_DeepFlavB_eta_v2.pkl'
#effFile2018DeepJet = 'TTLep_pow_2018_2j_2l_DeepFlavB_eta_v2.pkl'

#UL files

effFile2018DeepJet = 'TTSemiLep_pow_UL2018_2j_1l_DeepFlavB_eta_WP_L_v1.pkl'

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
sfFile  = 'btagging.json.gz' #Moriond2019

class BTagEfficiency:

    @staticmethod
    def getWeightDict_1b(effs, maxMultBTagWeight):
        '''Make Weight dictionary for jets
        '''
        zeroTagWeight = 1.

        for e in effs:
            zeroTagWeight*=(1-e)

        tagWeight={}
        for i in range(min(len(effs), maxMultBTagWeight)+1):
            tagWeight[i]=zeroTagWeight
            twfSum = 0.
            for tagged in itertools.combinations(effs, i):
                twf=1.
                for fac in [x/(1-x) for x in tagged]:
                    twf*=fac
                twfSum+=twf
            tagWeight[i]*=twfSum

        for i in range(maxMultBTagWeight+1):
            if not tagWeight.has_key(i):
                tagWeight[i] = 0.

        return tagWeight

    def getBTagSF_1a(self, var, bJets, nonBJets):
        if var not in self.btagWeightNames:
            raise ValueError( "Don't know what to do with b-tag variation %s" %var )
        if var != 'MC':
            ref = reduce(mul, [j['beff']['MC'] for j in bJets] + [1-j['beff']['MC'] for j in nonBJets], 1 )
            if ref>0:
                return reduce(mul, [j['beff'][var] for j in bJets] + [1-j['beff'][var] for j in nonBJets], 1 )/ref
            else:
                logger.warning( "getBTagSF_1a: MC efficiency is zero. Return SF 1. MC efficiencies: %r "% (  [j['beff']['MC'] for j in bJets] + [1-j['beff']['MC'] for j in nonBJets] ) )
                return 1


    def __init__( self, WP='M', fastSim=False, year='2018_UL', tagger='DeepJet' ):

        if year not in [ '2016preVFP_UL', '2016postVFP_UL', '2017_UL', '2018_UL']:
            raise Exception("Lepton SF for entered year is not known")

	self.SFDir   = '$CMSSW_BASE/src/jsonpog-integration/POG/BTV/'
        #self.dataDir = "$CMSSW_BASE/src/Analysis/Tools/data/btagEfficiencyData/"
	self.dataDir = "$CMSSW_BASE/src/VH/Tools/data/btagEfficiencyData"
        self.year = year
        self.tagger = tagger
	self.WP = WP

        # Whether or not FS SF are to be used
        self.fastSim = fastSim

        # All btag weight names per jet
        self.btagWeightNames = [ 'MC', 'SF', 'SF_b_Down', 'SF_b_Up', 'SF_l_Down', 'SF_l_Up' ]
        if self.fastSim:
            self.btagWeightNames += [ 'SF_FS_Up', 'SF_FS_Down']

        # Input files

	self.scaleFactorFile   = os.path.expandvars( os.path.join( self.SFDir, year, sfFile ) )
	self.scaleFactorFileFS = os.path.expandvars( os.path.join( self.SFDir, year, sfFile ) )

        if year == '2016preVFP_UL':
            self.etaBins = etaBins2016
            if tagger == 'DeepCSV':
                self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, effFile2016DeepCSV ) )
		self.measure = 'deepCSV'
            elif tagger == 'DeepJet':
                self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, effFile2016DeepJet ) )
		self.measure = 'deepJet'

	if year == '2016postVFP_UL':
	    self.etaBins = etaBins2016
	    if tagger == 'DeepCSV':
		self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, effFile2016DeepCSV ) )
		self.measure = 'deepCSV'
	    elif tagger == 'DeepJet':
		self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, effFile2016DeepJet ) )
		self.measure = 'deepJet'

        if year == '2017_UL':
            self.etaBins = etaBins2017
            if tagger == 'DeepCSV':
                self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, effFile2017DeepCSV ) )
		self.measure = 'deepCSV'
            elif tagger == 'DeepJet':
                self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, effFile2017DeepJet ) )
		self.measure = 'deepJet'

	if year == '2018_UL': 
            self.etaBins = etaBins2018
            if tagger == 'DeepCSV':
                self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, effFile2018DeepCSV ) )
		self.measure = 'deepCSV'
            elif tagger == 'DeepJet':
                self.mcEfficiencyFile  = os.path.expandvars( os.path.join( self.dataDir, effFile2018DeepJet ) )
		self.measure = 'deepJet'

	# Load SF 
        logger.info ( "Loading scale factors from %s", self.scaleFactorFile )
	
	if self.scaleFactorFile.endswith(".gz"):
		import gzip
		with gzip.open(self.scaleFactorFile, "rt") as f:
			data = f.read().strip()
		self.correction = _core.CorrectionSet.from_string(data)
	else:
		self.correction = _core.CorrectionSet.from_file(self.scaleFactorFile)
	
	#self.correction = _core.CorrectionSet.from_file(self.scaleFactorFile)
	print "Measurement type: ", self.measure," WP:",self.WP

        # Load MC efficiency
        logger.info( "Loading MC efficiency %s", self.mcEfficiencyFile )
        self.mcEff = pickle.load( file( self.mcEfficiencyFile ) )

    def getMCEff(self, pdgId, pt, eta):
        ''' Get MC efficiency for jet
        '''
	ptin = pt
	if pt<ptBins[0][0]:
		ptin = ptBins[0][0]
	elif pt>ptBins[len(ptBins)-1][0]:
		ptin = ptBins[len(ptBins)-1][0]

        for ptBin in ptBins:
            if ptin>=ptBin[0] and (ptin<ptBin[1] or ptBin[1]<0):
                for etaBin in self.etaBins:
                    if abs(eta)>=etaBin[0] and abs(eta)<etaBin[1]:
                        if abs(pdgId)==5:      return  self.mcEff[tuple(ptBin)][tuple(etaBin)]["b"]
                        elif abs(pdgId)==4:    return  self.mcEff[tuple(ptBin)][tuple(etaBin)]["c"]
                        else:                  return  self.mcEff[tuple(ptBin)][tuple(etaBin)]["other"]

        logger.debug( "No MC efficiency for pt %f eta %f pdgId %i", pt, eta, pdgId)
        return 1

    def getSF(self, pdgId, pt, eta):
        # BTag SF Not implemented below 20 GeV
        if pt<20: 
            if self.fastSim:
                return (1,1,1,1,1,1,1)
            else:
                return (1,1,1,1,1)

        # BTag SF Not implemented above absEta 2.4 (2016) & 2.5 (2017-18)
        if abs(eta)>=self.etaBins[0][1]: 
            if self.fastSim:
                return (1,1,1,1,1,1,1)
            else:
                return (1,1,1,1,1)

        #autobounds are implemented now, no doubling of uncertainties necessary anymore
        flavKey = toFlavourKey(pdgId)
        
        #FastSim SFs
        sf_fs   = 1 #if not self.fastSim else self.readerFS.eval_auto_bounds('central', flavKey, eta, pt)
        sf_fs_d = 1 #if not self.fastSim else self.readerFS.eval_auto_bounds('down',    flavKey, eta, pt)
        sf_fs_u = 1 #if not self.fastSim else self.readerFS.eval_auto_bounds('up',      flavKey, eta, pt)
        if sf_fs == 0:  # should not happen, however, if pt=1000 (exactly) the reader will return a sf of 0.
            sf_fs = 1
            sf_fs_u = 1
            sf_fs_d = 1
        
        #FullSim SFs (times FSSF)
        if abs(pdgId)==5 or abs(pdgId)==4: #SF for b/c
	    measure = self.measure + '_comb'
	    self.evaltr = self.correction[measure]
            sf      = sf_fs*self.evaltr.evaluate('central', self.WP, abs(pdgId) , abs(eta), pt)
            sf_b_d  = sf_fs*self.evaltr.evaluate('down',    self.WP, abs(pdgId) , abs(eta), pt)
            sf_b_u  = sf_fs*self.evaltr.evaluate('up',      self.WP, abs(pdgId) , abs(eta), pt)
            sf_l_d  = 1.
            sf_l_u  = 1.
        else: #SF for light flavours
	    measure = self.measure + '_incl'
	    self.evaltr = self.correction[measure]
            sf      = sf_fs*self.evaltr.evaluate('central', self.WP, abs(pdgId) , abs(eta), pt)
            sf_b_d  = 1.
            sf_b_u  = 1.
            sf_l_d  = sf_fs*self.evaltr.evaluate('down', self.WP, abs(pdgId) , abs(eta), pt)
            sf_l_u  = sf_fs*self.evaltr.evaluate('up', self.WP, abs(pdgId) , abs(eta), pt)

        if self.fastSim:
            return (sf, sf_b_d, sf_b_u, sf_l_d, sf_l_u, sf*sf_fs_u/sf_fs, sf*sf_fs_d/sf_fs)
        else:
            return (sf, sf_b_d, sf_b_u, sf_l_d, sf_l_u)

    def addBTagEffToJet(self, j):
        mcEff = self.getMCEff(j['hadronFlavour'], j['pt'], j['eta'])
        sf =    self.getSF(j['hadronFlavour'], j['pt'], j['eta'])
        if self.fastSim:
            j['beff'] =  {'MC':mcEff, 'SF':mcEff*sf[0], 'SF_b_Down':mcEff*sf[1], 'SF_b_Up':mcEff*sf[2], 'SF_l_Down':mcEff*sf[3], 'SF_l_Up':mcEff*sf[4], 'SF_FS_Up':mcEff*sf[5], 'SF_FS_Down':mcEff*sf[6]}
        else:
            j['beff'] =  {'MC':mcEff, 'SF':mcEff*sf[0], 'SF_b_Down':mcEff*sf[1], 'SF_b_Up':mcEff*sf[2], 'SF_l_Down':mcEff*sf[3], 'SF_l_Up':mcEff*sf[4]}

#Method 1d
#https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration
#https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
sfFile_1d = '$CMSSW_BASE/src/Analysis/Tools/data/btagEfficiencyData/ttH_BTV_CSVv2_13TeV_2015D_20151120.csv'
flavourSys_1d = {
    5:{'central', 'up_jes', 'down_jes', 'up_lf', 'down_lf', 'up_hfstats1', 'down_hfstats1', 'up_hfstats2', 'down_hfstats2'},
    4:{'central', 'up_cferr1', 'down_cferr1', 'up_cferr2', 'down_cferr2'},
    0:{'central', 'up_jes', 'down_jes', 'up_hf', 'down_hf', 'up_lfstats1', 'down_lfstats1', 'up_lfstats2', 'down_lfstats2'},
}
from operator import or_
'''
class btagEfficiency_1d:

    def addBTagEffToJet(self, j):
        j['beff'] = {sys: 1. if sys not in flavourSys_1d[abs(j['hadronFlavour'])] else self.readers[sys].eval(toFlavourKey(j['hadronFlavour']), j['eta'], j['pt'], j['btagCSV']) for sys in self.btagWeightNames}

    def __init__(self,  WP = ROOT.BTagEntry.OP_MEDIUM):
        self.btagWeightNames = reduce(or_, flavourSys_1d.values())

        self.scaleFactorFile = sfFile_1d
        logger.info( "Loading scale factors from %s", self.scaleFactorFile )
        self.calib = ROOT.BTagCalibration("csvv2", self.scaleFactorFile )
        self.readers = {sys: ROOT.BTagCalibrationReader(self.calib, ROOT.BTagEntry.OP_RESHAPING, "iterativefit", sys) for sys in self.btagWeightNames}

if __name__ == "__main__":
    print "2016"
    BTagEff = BTagEfficiency( year=2016,tagger="DeepCSV" )
    print BTagEff.getSF(5, 100, 1.5)[0]
    print BTagEff.getSF(5, 100, -1.5)[0]
    print BTagEff.getSF(5, 100, 2)[0]
    print BTagEff.getSF(5, 100, -2)[0]
    print BTagEff.getSF(5, 400, 1.5)[0]
    print BTagEff.getSF(5, 400, -1.5)[0]
    print BTagEff.getSF(5, 400, 2)[0]
    print BTagEff.getSF(5, 400, -2)[0]
    del BTagEff

    print "2017"
    BTagEff = BTagEfficiency( year=2017,tagger="DeepCSV" )
    print BTagEff.getSF(5, 100, 1.5)[0]
    print BTagEff.getSF(5, 100, -1.5)[0]
    print BTagEff.getSF(5, 100, 2)[0]
    print BTagEff.getSF(5, 100, -2)[0]
    print BTagEff.getSF(5, 400, 1.5)[0]
    print BTagEff.getSF(5, 400, -1.5)[0]
    print BTagEff.getSF(5, 400, 2)[0]
    print BTagEff.getSF(5, 400, -2)[0]
    del BTagEff

    print "2018"
    BTagEff = BTagEfficiency( year=2018,tagger="DeepCSV" )
    print BTagEff.getSF(5, 100, 1.5)[0]
    print BTagEff.getSF(5, 100, -1.5)[0]
    print BTagEff.getSF(5, 100, 2)[0]
    print BTagEff.getSF(5, 100, -2)[0]
    print BTagEff.getSF(5, 400, 1.5)[0]
    print BTagEff.getSF(5, 400, -1.5)[0]
    print BTagEff.getSF(5, 400, 2)[0]
    print BTagEff.getSF(5, 400, -2)[0]
    del BTagEff
'''
