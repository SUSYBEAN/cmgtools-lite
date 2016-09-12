import math
import ROOT
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi

class Pair(object):
    def __init__(self,leg1,leg2,pdg = 0):
        self.l1=ROOT.TLorentzVector()
        self.l1.SetPtEtaPhiM(leg1.pt(),leg1.eta(),leg1.phi(),leg1.mass())
        self.l2=ROOT.TLorentzVector()
        self.l2.SetPtEtaPhiM(leg2.p4().pt(),leg2.p4().eta(),leg2.p4().phi(),leg2.p4().mass())

	self.tunepl1=ROOT.TLorentzVector()
        self.useTuneP=3
        if hasattr(leg1,'TuneP_pt'):
            self.tunepl1.SetPtEtaPhiM(leg1.TuneP_pt(),leg1.TuneP_eta(),leg1.TuneP_phi(),leg1.TuneP_m())
        else: 
            self.tunepl1.SetPtEtaPhiM(leg1.pt(),leg1.eta(),leg1.phi(),leg1.mass())
            self.useTuneP-=1
        self.tunepl2=ROOT.TLorentzVector()
        if hasattr(leg2,'TuneP_pt'):
            self.tunepl2.SetPtEtaPhiM(leg2.TuneP_pt(),leg2.TuneP_eta(),leg2.TuneP_phi(),leg2.TuneP_m())
        else: 
            self.tunepl2.SetPtEtaPhiM(leg2.pt(),leg2.eta(),leg2.phi(),leg2.mass())
            self.useTuneP-=2

        self.leg1 = leg1
        self.leg2 = leg2
        self.pdg = pdg
        self.LV = self.l1+self.l2
        self.TuneP_LV = self.tunepl1+self.tunepl2
        et1 = math.sqrt(self.l1.M()*self.l1.M()+self.l1.Pt()*self.l1.Pt())
        et2 = math.sqrt(self.l1.M()*self.l1.M()+self.l2.Pt()*self.l2.Pt())
        TuneP_et1 = math.sqrt(self.tunepl1.M()*self.tunepl1.M()+self.tunepl1.Pt()*self.tunepl1.Pt())
        TuneP_et2 = math.sqrt(self.tunepl1.M()*self.tunepl1.M()+self.tunepl2.Pt()*self.tunepl2.Pt())

        self.MT  =math.sqrt(self.l1.M()*self.l1.M()+self.l1.M()*self.l1.M()+2*(et1*et2-self.l1.Px()*self.l2.Px()-self.l1.Py()*self.l2.Py()))
        self.TuneP_MT  =math.sqrt(2.0*self.tunepl1.M()*self.tunepl1.M()+2.0*(TuneP_et1*TuneP_et2-self.tunepl1.Px()*self.tunepl2.Px()-self.tunepl1.Py()*self.tunepl2.Py()))


    def rawP4(self):
        return self.l1+self.l2

    def p4(self):
        return self.LV

    def TuneP_p4(self):
        return self.TuneP_LV
    
    def m(self):
        return self.LV.M()

    def mass(self):
        return self.LV.M()
    
    def TuneP_m(self):
        return self.TuneP_LV.M()

    def TuneP_mass(self):
        return self.TuneP_LV.M()

    def pdgId(self):
        return self.pdg
    
    def mt2(self):
        return self.MT*self.MT

    def mt(self):
        return self.MT

    def TuneP_mt(self):
        return self.TuneP_MT

    def pt(self):
        return self.LV.Pt()

    def px(self):
        return self.LV.Px()

    def py(self):
        return self.LV.Py()

    def pz(self):
        return self.LV.Pz()

    def eta(self):
        return self.LV.Eta()

    def rapidity(self):
        return self.LV.Rapidity()

    def phi(self):
        return ROOT.TVector2.Phi_mpi_pi(self.LV.Phi())

    def TuneP_pt(self):
        return self.TuneP_LV.Pt()

    def TuneP_px(self):
        return self.TuneP_LV.Px()

    def TuneP_py(self):
        return self.TuneP_LV.Py()

    def TuneP_pz(self):
        return self.TuneP_LV.Pz()

    def TuneP_eta(self):
        return self.TuneP_LV.Eta()

    def TuneP_rapidity(self):
        return self.TuneP_LV.Rapidity()

    def TuneP_phi(self):
        return ROOT.TVector2.Phi_mpi_pi(self.TuneP_LV.Phi())

    def AbsdeltaPhi(self):
        return abs(deltaPhi(self.l1.Phi(),self.l2.Phi()))

    def AbsdeltaR(self):
        return abs(deltaR(self.l1.Eta(),self.l1.Phi(),self.l2.Eta(),self.l2.Phi()))

    def deltaPhi(self):
        return deltaPhi(self.l1.Phi(),self.l2.Phi())

    def TuneP_deltaPhi(self):
        return deltaPhi(self.l1.Phi(),self.l2.Phi())

    def deltaR(self):
        return deltaR(self.l1.Eta(),self.l1.Phi(),self.l2.Eta(),self.l2.Phi())

    def TuneP_deltaR(self):
        return deltaR(self.l1.Eta(),self.l1.Phi(),self.l2.Eta(),self.l2.Phi())

    def __getattr__(self, name):
        return getattr(self.LV,name)

