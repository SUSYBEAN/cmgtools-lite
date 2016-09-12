from CMGTools.XZZ2l2nu.analyzers.AutoFillTreeProducer  import * 
import math


LLType = NTupleObjectType("LLType", baseObjectTypes=[fourVectorType], variables = [
    #NTupleVariable("TuneP_usage",   lambda x : x.useTuneP, int), 
    #NTupleVariable("TuneP_pt",   lambda x : x.TuneP_pt(), float),               
    #NTupleVariable("TuneP_eta",   lambda x : x.TuneP_eta(), float),               
    #NTupleVariable("TuneP_phi",   lambda x : x.TuneP_phi(), float),               
    #NTupleVariable("TuneP_mass",   lambda x : x.TuneP_m(), float),               
    NTupleVariable("mt",   lambda x : x.mt(), float),       
    NTupleVariable("deltaPhi",   lambda x : x.deltaPhi(), float),       
    NTupleVariable("deltaR",   lambda x : x.deltaR(), float),       
    #NTupleVariable("TuneP_mt",   lambda x : x.TuneP_mt(), float),       
    #NTupleVariable("TuneP_deltaPhi",   lambda x : x.TuneP_deltaPhi(), float),       
    #NTupleVariable("TuneP_deltaR",   lambda x : x.TuneP_deltaR(), float),       
])



VVType = NTupleObjectType("VVType", baseObjectTypes=[], variables = [
  #NTupleSubObject("LV",  lambda x : x['pair'],fourVectorType),
  #NTupleVariable("TuneP_LV_pt",   lambda x : x['pair'].TuneP_pt(), float),       
  #NTupleVariable("TuneP_LV_eta",   lambda x : x['pair'].TuneP_eta(), float),       
  #NTupleVariable("TuneP_LV_phi",   lambda x : x['pair'].TuneP_phi(), float),       
  #NTupleVariable("TuneP_LV_mass",   lambda x : x['pair'].TuneP_m(), float),       
  NTupleVariable("deltaPhi",   lambda x : x['pair'].deltaPhi(), float),       
  #NTupleVariable("TuneP_deltaPhi",   lambda x : x['pair'].TuneP_deltaPhi(), float), 
  NTupleVariable("deltaR",   lambda x : x['pair'].deltaR(), float),       
  #NTupleVariable("TuneP_deltaR",   lambda x : x['pair'].TuneP_deltaR(), float),       
  NTupleVariable("mt",   lambda x : x['pair'].mt(), float),       
  #NTupleVariable("TuneP_mt",   lambda x : x['pair'].TuneP_mt(), float),       
])


LLNuNuType = NTupleObjectType("LLNuNuType", baseObjectTypes=[VVType], variables = [
    NTupleSubObject("l1",  lambda x : x['pair'].leg1,LLType),
    NTupleSubObject("l1_l1",  lambda x : x['pair'].leg1.leg1,leptonTypeExtra),
    NTupleSubObject("l1_l2",  lambda x : x['pair'].leg1.leg2,leptonTypeExtra),
    NTupleSubObject("l2",  lambda x : x['pair'].leg2,metType),
    #NTupleVariable("CosdphiZMet",   lambda x : math.cos(x['pair'].deltaPhi()), float), 
    #NTupleVariable("CosZdeltaPhi",   lambda x : math.cos(x['pair'].leg1.deltaPhi()), float), 
    #NTupleVariable("dPTPara",   lambda x : (x['pair'].leg1.pt() + x['pair'].leg2.pt() * math.cos(x['pair'].deltaPhi())), float), 
    #NTupleVariable("dPTParaRel",   lambda x : (x['pair'].leg1.pt() + x['pair'].leg2.pt() * math.cos(x['pair'].deltaPhi()))/(x['pair'].leg1.pt()), float), 
    #NTupleVariable("dPTPerp",   lambda x : (x['pair'].leg2.pt() * math.sin(x['pair'].deltaPhi())), float), 
    #NTupleVariable("dPTPerpRel",   lambda x : (x['pair'].leg2.pt() * math.sin(x['pair'].deltaPhi()))/(x['pair'].leg1.pt()), float), 
    #NTupleVariable("metOvSqSET",   lambda x : (x['pair'].leg2.pt())/math.sqrt(x['pair'].leg2.sumEt()), float), 
])

PhotonJetType = NTupleObjectType("PhotonJetType", baseObjectTypes=[], variables = [
    NTupleVariable("deltaPhi",   lambda x : x['pair'].deltaPhi(), float),
    NTupleVariable("mt",   lambda x : x['pair'].mt(), float),
    NTupleSubObject("l1",  lambda x : x['pair'].leg1,photonType),
    NTupleSubObject("l2",  lambda x : x['pair'].leg2,metType),
    NTupleVariable("CosDeltaPhi",   lambda x : math.cos(x['pair'].deltaPhi()), float),
    NTupleVariable("dPTPara",   lambda x : (x['pair'].leg1.pt() + x['pair'].leg2.pt() * math.cos(x['pair'].deltaPhi())), float),
    NTupleVariable("dPTParaRel",   lambda x : (x['pair'].leg1.pt() + x['pair'].leg2.pt() * math.cos(x['pair'].deltaPhi()))/(x['pair'].leg1.pt()), float),
    NTupleVariable("dPTPerp",   lambda x : (x['pair'].leg2.pt() * math.sin(x['pair'].deltaPhi())), float),
    NTupleVariable("dPTPerpRel",   lambda x : (x['pair'].leg2.pt() * math.sin(x['pair'].deltaPhi()))/(x['pair'].leg1.pt()), float),
    NTupleVariable("metOvSqSET",   lambda x : (x['pair'].leg2.pt())/math.sqrt(x['pair'].leg2.sumEt()), float),
])


llpairType = NTupleObjectType("llpairType", baseObjectTypes=[], variables = [
    NTupleSubObject("Z",  lambda x : x,LLType),
    NTupleSubObject("l1",  lambda x : x.leg1,leptonType),
    NTupleSubObject("l2",  lambda x : x.leg2,leptonType),
])

metCompType = NTupleObjectType("metCompType", baseObjectTypes=[], variables = [
    NTupleVariable("Px",  lambda x : x[0]),
    NTupleVariable("Py",  lambda x : x[1]),
    NTupleVariable("Et",  lambda x : x[2]),
])

jet4metType = NTupleObjectType("jet4metType", baseObjectTypes=[], variables = [
    NTupleSubObject("rawP4forT1",  lambda x : x["rawP4forT1"], metCompType),
    NTupleSubObject("type1METCorr",  lambda x : x["type1METCorr"], metCompType),
    NTupleSubObject("corrP4forT1",  lambda x : x["corrP4forT1"], metCompType),
])

JetType = NTupleObjectType("xzzJetType", baseObjectTypes=[jetType], variables = [
    NTupleVariable("id",    lambda x : x.jetID("POG_PFID") , int, mcOnly=False,help="POG Loose jet ID"),
    NTupleVariable("area",   lambda x : x.jetArea(), help="Catchment area of jet"),
    NTupleVariable("rawFactor",   lambda x : x.rawFactor(), float, help="pt/rawfactor will give you the raw pt"),
    #NTupleVariable("corr_jer",  lambda x : getattr(x, 'corrJER', -99), float, mcOnly=True, help="JER corr factor"),
    NTupleVariable("btagCSV",   lambda x : x.btag('pfCombinedInclusiveSecondaryVertexV2BJetTags'), help="CSV-IVF v2 discriminator"),
    NTupleVariable("btagCMVA",  lambda x : x.btag('pfCombinedMVABJetTags'), help="CMVA discriminator"),
    NTupleVariable("mcPt",   lambda x : x.matchedGenJet.pt() if getattr(x,"matchedGenJet",None) else 0., mcOnly=True, help="p_{T} of associated gen jet"),
    NTupleVariable("mcEta",   lambda x : x.matchedGenJet.eta() if getattr(x,"matchedGenJet",None) else 0., mcOnly=True, help="eta of associated gen jet"),
    NTupleVariable("mcPhi",   lambda x : x.matchedGenJet.phi() if getattr(x,"matchedGenJet",None) else 0., mcOnly=True, help="phi of associated gen jet"),
    NTupleVariable("mcMass",   lambda x : x.matchedGenJet.mass() if getattr(x,"matchedGenJet",None) else 0., mcOnly=True, help="mass of associated gen jet"),
    NTupleVariable("mcFlavour", lambda x : x.partonFlavour(), int,     mcOnly=True, help="parton flavour (physics definition, i.e. including b's from shower)"),
    #NTupleVariable("btag",   lambda x : x.bTag(), float),
    #NTupleVariable("nConstituents",   lambda x : len(x.constituents), int),
    #NTupleVariable("looseID",   lambda x : x.looseID, int),
    #NTupleVariable("tightID",   lambda x : x.tightID, int),
    NTupleVariable("chargedHadronEnergyFraction",   lambda x : x.chargedHadronEnergyFraction(), float,  help="for Jet ID"),
    NTupleVariable("neutralHadronEnergyFraction",   lambda x : x.neutralHadronEnergyFraction(), float, help="for Jet ID"),
    NTupleVariable("neutralEmEnergyFraction",   lambda x : x.neutralEmEnergyFraction(), float, help="for Jet ID"),
    NTupleVariable("muonEnergyFraction",   lambda x : x.muonEnergyFraction(), float, help="for Jet ID"),
    NTupleVariable("chargedEmEnergyFraction",   lambda x : x.chargedEmEnergyFraction(), float, help="for Jet ID"),
    NTupleVariable("chargedHadronMultiplicity",   lambda x : x.chargedHadronMultiplicity(), float, help="for Jet ID"),
    NTupleVariable("chargedMultiplicity",   lambda x : x.chargedMultiplicity(), float, help="for Jet ID"),
    NTupleVariable("neutralMultiplicity",   lambda x : x.neutralMultiplicity(), float, help="for Jet ID"),
])

JetTypeExtra = NTupleObjectType("xzzJetTypeExtra", baseObjectTypes=[JetType], variables = [
    NTupleVariable("photonEnergyFraction",   lambda x : x.photonEnergyFraction(), float,),
    NTupleVariable("HFHadronEnergyFraction",   lambda x : x.HFHadronEnergyFraction(), float,),
    NTupleVariable("HFEMEnergyFraction",   lambda x : x.HFEMEnergyFraction(), float),
    NTupleVariable("electronEnergyFraction",   lambda x : x.electronEnergyFraction(), float),
])

