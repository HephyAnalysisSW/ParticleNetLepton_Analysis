#Standard imports
import  ROOT
from ROOT import TVector3, TLorentzVector 
from    math    import pi, sqrt, cos, sin, sinh, log, cosh, acos, acosh
from    array   import array

def neutrino_mom(vec_lep, MET_pt, MET_phi, random):

	W_mass = 80.4

	pnu = TLorentzVector()

	if vec_lep.E()<0:
		pnu.etPtEtaPhiM(0,-100,-100,0)
	else:
		mT2 = 2*vec_lep.Pt()*MET_pt*(1-cos(deltaPhi(vec_lep.Phi(),MET_phi)))
		Delta2 = (W_mass*W_mass - mT2)*1./(2*MET_pt*vec_lep.Pt())
		if (Delta2>=0):
			try:
				nueta = (vec_lep.Eta() + abs(acosh(1+(Delta2)))) if (random>=0.5) else (vec_lep.Eta() - abs(acosh(1+(Delta2))))
			except Exception:
				pass
				nueta = -100
			pnu.SetPtEtaPhiM(MET_pt,nueta,MET_phi,0)
		else:
			#pnu.SetPtEtaPhiM(0,-100,-100,0)
			nueta = vec_lep.Eta()
			pnu.SetPtEtaPhiM(MET_pt,nueta,MET_phi,0)
	
	return pnu

def neutrino_mom_both(vec_lep, MET_pt, MET_phi, W_mass=80.4):

	#W_mass = 80.4
	pnu = [TLorentzVector(), TLorentzVector()]
	nueta = [-100,-100]

	if vec_lep.E()<0:
		pnu[0].SetPtEtaPhiM(0,-100,-100,0)
		pnu[1].SetPtEtaPhiM(0,-100,-100,0)
	else:
		mT2 = 2*vec_lep.Pt()*MET_pt*(1-cos(deltaPhi(vec_lep.Phi(),MET_phi)))
		Delta2 = (W_mass*W_mass - mT2)*1./(2*MET_pt*vec_lep.Pt())
		if (Delta2>=0):
			try:
				nueta[0] = (vec_lep.Eta() + abs(acosh(1+(Delta2))))
				nueta[1] = (vec_lep.Eta() - abs(acosh(1+(Delta2))))
			except Exception:
				pass
				pass
		else:
			nueta[0] = vec_lep.Eta()
			nueta[1] = vec_lep.Eta()

	pnu[0].SetPtEtaPhiM(MET_pt,nueta[0],MET_phi,0)
	pnu[1].SetPtEtaPhiM(MET_pt,nueta[1],MET_phi,0)

	return pnu
