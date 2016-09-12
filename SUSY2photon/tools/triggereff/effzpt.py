from ROOT import *
gROOT.SetBatch()

def gethist(flnm,att,a,b,c,st=''):
    lst=[[0.0,TH1F(st+str(i),st+str(i),a,b,c)] for i in range(4)]
    f=TFile(flnm)
    t=f.Get('tree')
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        if t.llnunu_l1_mass[0]<70 or t.llnunu_l1_mass[0]>110 or t.llnunu_l1_pt[0]<200: continue
        if abs(t.llnunu_l1_l1_pdgId[0])==11 and t.llnunu_l1_l1_pt[0]>115:
            lst[0][0]+=1 if t.isData else t.puWeight
            lst[0][1].Fill(getattr(t,att)[0],1 if t.isData else t.puWeight)
            if t.llnunu_l1_l1_trigerob_pt[0]!=-100 or t.llnunu_l1_l2_trigerob_pt[0]!=-100:
                lst[1][0]+=1 if t.isData else t.puWeight
                lst[1][1].Fill(getattr(t,att)[0],1 if t.isData else t.puWeight)
        if abs(t.llnunu_l1_l1_pdgId[0])==13 and ((t.llnunu_l1_l1_pt[0]>50 and abs(t.llnunu_l1_l1_eta[0])<2.1) or (t.llnunu_l1_l2_pt[0]>50 and abs(t.llnunu_l1_l2_eta[0])<2.1)) and (t.llnunu_l1_l1_highPtID[0] or t.llnunu_l1_l2_highPtID[0]):
            lst[2][0]+=1 if t.isData else t.puWeight
            lst[2][1].Fill(getattr(t,att)[0],1 if t.isData else t.puWeight)
            if t.llnunu_l1_l1_trigerob_pt[0]!=-100 or t.llnunu_l1_l2_trigerob_pt[0]!=-100:
                lst[3][0]+=1 if t.isData else t.puWeight
                lst[3][1].Fill(getattr(t,att)[0],1 if t.isData else t.puWeight)
    f.Close()
    lst.append((lst[1][0]/lst[0][0],ercal(lst[1][0],lst[0][0]),lst[3][0]/lst[2][0],ercal(lst[3][0],lst[2][0])))
    print 'ele efficiency: {0}/{1}={2}% +- {3}%'.format(lst[1][0],lst[0][0],100*lst[4][0],lst[4][1]*100)
    print 'muon efficiency: {0}/{1}={2}% +- {3}%'.format(lst[3][0],lst[2][0],100*lst[4][2],lst[4][3]*100)
    return lst
def ercal(a0,b0):
    a1=b0-a0
    return (a1*a0/b0**3)**.5
def geteff(a,b,n=4,nm=''):
    if a.GetSize()!=b.GetSize():
        print 'ERROR: Different bin number!'
        return TH1F()
    e=TH1F(nm,nm,a.GetSize()-2,a.GetXaxis().GetXmin(),a.GetXaxis().GetXmax())
    for i in range(1,a.GetSize()-1): 
        a0=a.GetBinContent(i)
        b0=b.GetBinContent(i)
        if b0:
            e.SetBinContent(i,100.*a0/b0)
            a1=b0-a0
            e.SetBinError(i,(a1*a0/b0**3)**.5*100)
    e.SetMarkerColor(n)
    e.SetLineColor(n)
    return e
if __name__=='__main__':
    gStyle.SetOptStat(000000)
    gStyle.SetLegendBorderSize(0)
    c=TCanvas('c','c',1600,1600)
    leg=TLegend(.91,.8,1,.9)
    ag='/afs/cern.ch/work/y/yanchu/graviton/CMSSW_7_6_3_patch2/src/CMGTools/XZZ2l2nu/cfg/trgdt/JetHT_Run2015D_16Dec/vvTreeProducer/tree.root' # hadd all data tree.root files
    a=gethist(ag,'llnunu_l1_pt',30,200,800,'dt')
    ee=geteff(a[1][1],a[0][1],2,'e_data')
    em=geteff(a[3][1],a[2][1],2,'mu_data')
    ee.SetTitle('electron efficiency vs Pt_Z;Pt_Z [GeV];efficiency [%]')
    em.SetTitle('muon efficiency vs Pt_Z;Pt_Z [GeV];efficiency [%]')

    agmc='/afs/cern.ch/work/y/yanchu/graviton/CMSSW_7_6_3_patch2/src/CMGTools/XZZ2l2nu/cfg/trgmc/dy.root' #hadd DY samples
    amc=gethist(agmc,'llnunu_l1_pt',30,200,800,'dy')
    eemc=geteff(amc[1][1],amc[0][1],1,'e_DY')
    emmc=geteff(amc[3][1],amc[2][1],1,'mu_DY')

    agzmc='/afs/cern.ch/work/y/yanchu/graviton/CMSSW_7_6_3_patch2/src/CMGTools/XZZ2l2nu/cfg/trgmc/zz2l2nu.root' #hadd ZZ2l2nu samples
    azmc=gethist(agzmc,'llnunu_l1_pt',30,200,800,'zz')
    eezmc=geteff(azmc[1][1],azmc[0][1],7,'e_zz')
    emzmc=geteff(azmc[3][1],azmc[2][1],7,'mu_zz')
 
    print 'electron SF (DY): {0}+-{1}'.format(a[4][0]/amc[4][0],((a[4][1]/a[4][0])**2+(amc[4][1]/amc[4][0])**2)**.5*a[4][0]/amc[4][0])
    print 'muon SF (DY): {0}+-{1}'.format(a[4][2]/amc[4][2],((a[4][3]/a[4][2])**2+(amc[4][3]/amc[4][2])**2)**.5*a[4][2]/amc[4][2])
#    print 'electron SF (signal): {0}+-{1}'.format(a[4][0]/asmc[4][0],((a[4][1]/a[4][0])**2+(asmc[4][1]/asmc[4][0])**2)**.5*a[4][0]/asmc[4][0])
#    print 'muon SF (signal): {0}+-{1}'.format(a[4][2]/asmc[4][2],((a[4][3]/a[4][2])**2+(asmc[4][3]/asmc[4][2])**2)**.5*a[4][2]/asmc[4][2])
    print 'electron SF (ZZ): {0}+-{1}'.format(a[4][0]/azmc[4][0],((a[4][1]/a[4][0])**2+(azmc[4][1]/azmc[4][0])**2)**.5*a[4][0]/azmc[4][0])
    print 'muon SF (ZZ): {0}+-{1}'.format(a[4][2]/azmc[4][2],((a[4][3]/a[4][2])**2+(azmc[4][3]/azmc[4][2])**2)**.5*a[4][2]/azmc[4][2])
    ee.Draw('e')
    eemc.Draw('esame')
    eezmc.Draw('esame')
#    c.Print('mcelezpt.pdf')
#    c.Print('mcelezpt.png')
    leg.AddEntry(eemc,'DY sample')
    leg.AddEntry(eezmc,'ZZ2l2nu sample')
    leg.AddEntry(ee,'2015D data')
    leg.Draw('same')
    c.Print('elezpt.pdf')
    c.Print('elezpt.png')
    em.Draw('e')
    emmc.Draw('esame')
    emzmc.Draw('esame')
    leg.Draw('same')
    c.Print('muzpt.pdf')
    c.Print('muzpt.png')
    esc=TH1F('esc','esc',30,200,800)
    esc.SetTitle('electron scale factor vs Pt_Z;Pt_Z [GeV];SF')
    for i in range(1,esc.GetSize()-1):
        if ee.GetBinContent(i) and eemc.GetBinContent(i):
            esc.SetBinContent(i,ee.GetBinContent(i)/eemc.GetBinContent(i))
            esc.SetBinError(i,((ee.GetBinError(i)/ee.GetBinContent(i))**2+(eemc.GetBinError(i)/eemc.GetBinContent(i))**2)**.5*esc.GetBinContent(i))
    esc.Draw('e')
    c.Print('esfzpt.pdf')
    c.Print('esfzpt.png')
    msc=TH1F('msc','msc',30,200,800)
    msc.SetTitle('muon scale factor vs Pt_Z;Pt_Z [GeV];SF')
    for i in range(1,esc.GetSize()-1):
        if em.GetBinContent(i) and emmc.GetBinContent(i):
            msc.SetBinContent(i,em.GetBinContent(i)/emmc.GetBinContent(i))
            msc.SetBinError(i,((em.GetBinError(i)/em.GetBinContent(i))**2+(emmc.GetBinError(i)/emmc.GetBinContent(i))**2)**.5*msc.GetBinContent(i))
    msc.Draw('e')
    c.Print('msfzpt.pdf')
    c.Print('msfzpt.png')

