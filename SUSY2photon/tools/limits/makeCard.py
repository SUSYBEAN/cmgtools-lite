#!/bin/env python

import optparse
from ROOT import *
from XZZDataCardMaker import XZZDataCardMaker 
from Setup import *

#parser = optparse.OptionParser()
#parser.add_option("-m","--masses",dest="masses",default='limitPlot',help="Limit plot")


#tag='sig1pb'
tag='xzz2l2nu'

outdir='cards'

tag = outdir+'/'+tag

infile_name='tightzpt100met200_puWeight_metfilter_unblind_both_log_scale1_.root'

rfin = TFile(infile_name)

masses = [
         600,
         800,
          1000,
         1200,1400,
         1600,
         1800,
          2000,
         2500,
          3000,
         3500,
          4000,
#         4500
         ]

# common datacard inputs
channels=['ee','mm']
#luminosity=2.1691267
luminosity=2.318278305

# systematics
#lumi
lumi_unc=1.027
sys_lumi_values = {'ee': {'xzz':lumi_unc, 'zjets':lumi_unc, 'vvzreso':lumi_unc, 'vvnonreso':lumi_unc, 'tt':lumi_unc,},
                   'mm': {'xzz':lumi_unc, 'zjets':lumi_unc, 'vvzreso':lumi_unc, 'vvnonreso':lumi_unc, 'tt':lumi_unc,},
                  }

#error = Double(0.)

# get histograms

# data_obs
hObs = {}
hObs['ee'] = rfin.Get("h_data_obs_el")
hObs['mm'] = rfin.Get("h_data_obs_mu")

# ZJets
hZJets = {}
hZJets['ee'] = rfin.Get("h_ZJets_el") 
hZJets['mm'] = rfin.Get("h_ZJets_mu") 

# VVZReso
hVVZReso = {}
hVVZReso['ee'] = rfin.Get("h_VVZReso_el")
hVVZReso['mm'] = rfin.Get("h_VVZReso_mu")

# VVNonReso
hVVNonReso = {}
hVVNonReso['ee'] = rfin.Get("h_VVNonReso_el")
hVVNonReso['mm'] = rfin.Get("h_VVNonReso_mu")

# TT
hTT = {}
hTT['ee'] = rfin.Get("h_TT_el")
hTT['mm'] = rfin.Get("h_TT_mu")

# XZZ
hSigs = {}
for mass in masses:
   hsig = {}
   hsig['ee'] = rfin.Get("h_BulkGravToZZToZlepZinv_narrow_"+str(mass)+"_el")
   hsig['mm'] = rfin.Get("h_BulkGravToZZToZlepZinv_narrow_"+str(mass)+"_mu")
   hSigs[mass] = hsig


# get integral and error for list of channels with cuts from a hists list for a given process
def GetIntegralAndError3D(pname,channels,hists,cutx,cuty,cutz):

    error = Double(0.)
    yields = {}
    yields_err = {}
    for ch in channels:
        x_bin_min = hists[ch].GetXaxis().FindBin(cutx+0.1)
        x_bin_max = hists[ch].GetXaxis().GetNbins()
        y_bin_min = hists[ch].GetYaxis().FindBin(cuty+0.1)
        y_bin_max = hists[ch].GetYaxis().GetNbins()
        z_bin_min = hists[ch].GetZaxis().FindBin(cutz+0.1)
        z_bin_max = hists[ch].GetZaxis().GetNbins()
        yields[ch] = hists[ch].IntegralAndError(x_bin_min,x_bin_max,y_bin_min,y_bin_max,z_bin_min,z_bin_max,error)
        yields_err[ch] = {pname:1.0 if yields[ch]<=0 else 1.0+abs(error)/yields[ch]}
        if yields[ch]<=0.0 : yields[ch] = 1e-30 # protect against negative or zero yields
        print '  GetIntegralAndError::'+pname+'['+str(cutx)+','+str(cuty)+','+str(cutz)+']['+ch+']: yields={0:}, err={1:}, err/yields={2:}'.format(yields[ch], error, error/yields[ch])

    return yields,yields_err

# get integral and error for list of channels with cuts from a hists list for a given process
def GetIntegralAndError(pname,channels,hists,cutx):

    error = Double(0.)
    yields = {}
    yields_err = {}
    for ch in channels:
        x_bin_min = hists[ch].GetXaxis().FindBin(cutx+0.1)
        x_bin_max = hists[ch].GetXaxis().GetNbins()
        yields[ch] = hists[ch].IntegralAndError(x_bin_min,x_bin_max,error)
        yields_err[ch] = {pname:1.0 if yields[ch]<=0 else 1.0+abs(error)/yields[ch]}
        if yields[ch]<=0.0 : yields[ch] = 1e-30 # protect against negative or zero yields
        print '  GetIntegralAndError::'+pname+'['+str(cutx)+']['+ch+']: yields={0:}, err={1:}, err/yields={2:}'.format(yields[ch], error, error/yields[ch])

    return yields,yields_err

# make datacards for all mass points
for mass in masses:
    # loop mt cuts
    for mtcut in mtcuts[mass]:

        # card tag
        cardTag = '{tag}_m{mass:n}_mtcut{mtcut:n}'.format(tag=tag,mass=mass,mtcut=mtcut)
        #
        print '# Makeing card: '+cardTag+'.txt'

        # cardMaker
        cardMaker = XZZDataCardMaker(channels,luminosity,cardTag)

        # observation
        obs,obs_err = GetIntegralAndError('obs',channels,hObs,mtcut)
        cardMaker.addObservation(obs)

        # signal
        signal,signal_err = GetIntegralAndError('xzz',channels,hSigs[mass],mtcut)
        cardMaker.addContribution('xzz',0,signal)
        cardMaker.addSystematic('xzz_stat','lnN',signal_err)


        # zjets
        zjets,zjets_err = GetIntegralAndError('zjets',channels,hZJets,mtcut)
        cardMaker.addContribution('zjets',1,zjets)
        cardMaker.addSystematic('zjets_stat','lnN',zjets_err)

        # vvzreso
        vvzreso,vvzreso_err = GetIntegralAndError('vvzreso',channels,hVVZReso,mtcut)
        cardMaker.addContribution('vvzreso',2,vvzreso)
        cardMaker.addSystematic('vvzreso_stat','lnN',vvzreso_err)

        # vvnonreso
        vvnonreso,vvnonreso_err = GetIntegralAndError('vvnonreso',channels,hVVNonReso,mtcut)
        cardMaker.addContribution('vvnonreso',3,vvnonreso)
        cardMaker.addSystematic('vvnonreso_stat','lnN',vvnonreso_err)

        # tt
        tt,tt_err = GetIntegralAndError('tt',channels,hTT,mtcut)
        cardMaker.addContribution('tt',4,tt)
        cardMaker.addSystematic('tt_stat','lnN',tt_err)

        # other systematics
        # lumi
        cardMaker.addSystematic('lumi_unc','lnN',sys_lumi_values)

        # counts
        Nsig = 0
        Nbkg = 0
        for ch in channels: Nsig += signal[ch]
        for ch in channels: Nbkg += zjets[ch]
        for ch in channels: Nbkg += vvzreso[ch]
        for ch in channels: Nbkg += vvnonreso[ch]
        for ch in channels: Nbkg += tt[ch]
        print " yields [nSig, nBkg] : & "+str(Nsig)+" & "+str(Nbkg)

        # makeCard
        cardMaker.makeCard() 

    # end loop mtcuts
# end loop masses
      




