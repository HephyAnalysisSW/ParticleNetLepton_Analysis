#Standard imports
import  ROOT
from ROOT import TVector3, TLorentzVector 
from    math    import pi, sqrt, cos, sin, sinh, log, cosh, acos, acosh

def getTheta(lep1, lep2, H):

	beam = TLorentzVector()

	tmp_lep1, tmp_lep2, tmp_H = TLorentzVector(), TLorentzVector(), TLorentzVector()

	tmp_lep1.SetPtEtaPhiM(lep1.Pt(),lep1.Eta(),lep1.Phi(),lep1.M())
	tmp_lep2.SetPtEtaPhiM(lep2.Pt(),lep2.Eta(),lep2.Phi(),lep2.M())
	tmp_H.SetPtEtaPhiM(H.Pt(),H.Eta(),H.Phi(),H.M())

	if(lep1.Eta()<-10 or lep2.Eta()<-10 or tmp_H.Eta()<-10):
		return -100
	
	beam.SetPxPyPzE(0,0,6500,6500)
	      			
	V_mom, bVH = TLorentzVector(), TVector3()
	V_mom = tmp_lep1+tmp_lep2
	bVH = (tmp_lep1+tmp_lep2+tmp_H).BoostVector()

	V_mom.Boost(-bVH)

	Theta = float('nan')

	try:
		Theta = (V_mom.Vect().Unit()).Angle(beam.Vect().Unit())
		#Theta  = acos((V_mom.Vect().Unit()).Dot(beam.Vect().Unit()))
	except Exception:
		pass
		Theta = -100

	return Theta;

def gettheta(lep1, lep2, H):

	tmp_lep1, tmp_lep2, tmp_H = TLorentzVector(), TLorentzVector(), TLorentzVector()

	tmp_lep1.SetPtEtaPhiM(lep1.Pt(),lep1.Eta(),lep1.Phi(),lep1.M())
	tmp_lep2.SetPtEtaPhiM(lep2.Pt(),lep2.Eta(),lep2.Phi(),lep2.M())
	tmp_H.SetPtEtaPhiM(H.Pt(),H.Eta(),H.Phi(),H.M())

	if(lep1.Eta()<-10 or lep2.Eta()<-10 or tmp_H.Eta()<-10):
		return -100
	
	V_mom, bVH, bV = TLorentzVector(), TVector3(), TVector3()

	bVH = (tmp_lep1 + tmp_lep2 + tmp_H).BoostVector()
	V_mom = (tmp_lep1 + tmp_lep2)

	V_mom.Boost(-bVH)
	tmp_lep1.Boost(-bVH)

	bV = V_mom.BoostVector()
	tmp_lep1.Boost(-bV)

	theta = float('nan')
	try:
		theta = (V_mom).Angle(tmp_lep1.Vect())
	except Exception:
		pass
		theta = -100

	return theta

def getphi(lep1, lep2, H):

	beam = TLorentzVector()

	tmp_lep1, tmp_lep2, tmp_H = TLorentzVector(), TLorentzVector(), TLorentzVector()

	tmp_lep1.SetPtEtaPhiM(lep1.Pt(),lep1.Eta(),lep1.Phi(),lep1.M())
	tmp_lep2.SetPtEtaPhiM(lep2.Pt(),lep2.Eta(),lep2.Phi(),lep2.M())
	tmp_H.SetPtEtaPhiM(H.Pt(),H.Eta(),H.Phi(),H.M())

	if(lep1.Eta()<-10 or lep2.Eta()<-10 or tmp_H.Eta()<-10):
		return -100

	beam.SetPxPyPzE(0,0,6500,6500)

	V_mom, bVH, n_scatter, n_decay = TLorentzVector(), TVector3(), TVector3(), TVector3()
	bVH = (tmp_lep1+tmp_lep2+tmp_H).BoostVector()
	V_mom = tmp_lep1+tmp_lep2

	tmp_lep1.Boost(-bVH)
	tmp_lep2.Boost(-bVH)
	V_mom.Boost(-bVH)
	#beam.Boost(-bVH)

	n_scatter = ((beam.Vect().Unit()).Cross(V_mom.Vect())).Unit()
	n_decay   = (tmp_lep1.Vect().Cross(tmp_lep2.Vect())).Unit()

	sign_flip =  -1 if ( ((n_scatter.Cross(n_decay))*(V_mom.Vect())) < 0 ) else +1

	try:
		phi = sign_flip*acos(n_scatter.Dot(n_decay))	
	except Exception:
		pass
		phi = -100

	return phi


