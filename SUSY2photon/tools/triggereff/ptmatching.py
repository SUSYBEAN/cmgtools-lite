from ROOT import *
gROOT.SetBatch()

def deltR(e1,e2,p1,p2):
    pi=3.1415926536
    de=e1-e2
    dp=p1-p2
    while dp>pi: dp-=2*pi
    while dp<-pi: dp+=2*pi
    return (de**2+dp**2)**.5

def gethist(flnm):
    lst=[[TH1F('','',200,-100,100),TH1F('','',4,0,4),TH1F('','',200,0,0.3)] for i in range(2)]
    f=TFile(flnm)
    t=f.Get('tree')
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        if t.llnunu_l1_mass[0]<70 or t.llnunu_l1_mass[0]>110 or t.llnunu_l1_pt[0]<200: continue
        if abs(t.llnunu_l1_l1_pdgId[0])==11 and t.llnunu_l1_l1_pt[0]>115:
            if t.llnunu_l1_l1_trigerob_pt[0]!=-100 :
                lst[0][0].Fill(t.llnunu_l1_l1_trigerob_pt[0]-t.llnunu_l1_l1_pt[0])
                lst[0][2].Fill(deltR(t.llnunu_l1_l1_trigerob_eta[0],t.llnunu_l1_l1_eta[0],t.llnunu_l1_l1_trigerob_phi[0],t.llnunu_l1_l1_phi[0]))
                if t.llnunu_l1_l2_trigerob_pt[0]!=-100:
                    lst[0][0].Fill(t.llnunu_l1_l2_trigerob_pt[0]-t.llnunu_l1_l2_pt[0])
                    lst[0][2].Fill(deltR(t.llnunu_l1_l2_trigerob_eta[0],t.llnunu_l1_l2_eta[0],t.llnunu_l1_l2_trigerob_phi[0],t.llnunu_l1_l2_phi[0]))
                    lst[0][1].Fill(1)
                else:lst[0][1].Fill(2)
            else:
                if t.llnunu_l1_l2_trigerob_pt[0]!=-100:
                    lst[0][0].Fill(t.llnunu_l1_l2_trigerob_pt[0]-t.llnunu_l1_l2_pt[0])
                    lst[0][2].Fill(deltR(t.llnunu_l1_l2_trigerob_eta[0],t.llnunu_l1_l2_eta[0],t.llnunu_l1_l2_trigerob_phi[0],t.llnunu_l1_l2_phi[0]))
                    lst[0][1].Fill(3)
                else:lst[0][1].Fill(0)

        if abs(t.llnunu_l1_l1_pdgId[0])==13 and ((t.llnunu_l1_l1_pt[0]>50 and abs(t.llnunu_l1_l1_eta[0])<2.1) or (t.llnunu_l1_l2_pt[0]>50 and abs(t.llnunu_l1_l2_eta[0])<2.1)) and (t.llnunu_l1_l1_highPtID[0] or t.llnunu_l1_l2_highPtID[0]):
            if t.llnunu_l1_l1_trigerob_pt[0]!=-100:
                lst[1][0].Fill(t.llnunu_l1_l1_trigerob_pt[0]-t.llnunu_l1_l1_pt[0])
                lst[1][2].Fill(deltR(t.llnunu_l1_l1_trigerob_eta[0],t.llnunu_l1_l1_eta[0],t.llnunu_l1_l1_trigerob_phi[0],t.llnunu_l1_l1_phi[0]))
                if t.llnunu_l1_l2_trigerob_pt[0]!=-100:
                    lst[1][0].Fill(t.llnunu_l1_l2_trigerob_pt[0]-t.llnunu_l1_l2_pt[0])
                    lst[1][2].Fill(deltR(t.llnunu_l1_l2_trigerob_eta[0],t.llnunu_l1_l2_eta[0],t.llnunu_l1_l2_trigerob_phi[0],t.llnunu_l1_l2_phi[0]))
                    lst[1][1].Fill(1)
                else:lst[1][1].Fill(2)
            else:
                if t.llnunu_l1_l2_trigerob_pt[0]!=-100:
                    lst[1][0].Fill(t.llnunu_l1_l2_trigerob_pt[0]-t.llnunu_l1_l2_pt[0])
                    lst[1][2].Fill(deltR(t.llnunu_l1_l2_trigerob_eta[0],t.llnunu_l1_l2_eta[0],t.llnunu_l1_l2_trigerob_phi[0],t.llnunu_l1_l2_phi[0]))
                    lst[1][1].Fill(3)
                else:lst[1][1].Fill(0)
    f.Close()
    print 'e selected:\t',lst[0][1].GetEntries()
    print 'e selected with match:\t',lst[0][1].GetEntries()-lst[0][1].GetBinContent(1)
    print "matched e1 and e2:\t",lst[0][1].GetBinContent(2),lst[0][1].GetBinContent(2)/(lst[0][1].GetEntries()-lst[0][1].GetBinContent(1))
    print "matched e1 but not e2:\t",lst[0][1].GetBinContent(3),lst[0][1].GetBinContent(3)/(lst[0][1].GetEntries()-lst[0][1].GetBinContent(1))
    print "matched e2 but not e1:\t",lst[0][1].GetBinContent(4),lst[0][1].GetBinContent(4)/(lst[0][1].GetEntries()-lst[0][1].GetBinContent(1))
    print 'matched neither e1 nor e2:\t',lst[0][1].GetBinContent(1)

    print 'm selected:\t',lst[1][1].GetEntries()
    print 'm selected with match:\t',lst[1][1].GetEntries()-lst[1][1].GetBinContent(1)
    print "matched m1 and m2:\t",lst[1][1].GetBinContent(2),lst[1][1].GetBinContent(2)/(lst[1][1].GetEntries()-lst[1][1].GetBinContent(1))
    print "matched m1 but not m2:\t",lst[1][1].GetBinContent(3),lst[1][1].GetBinContent(3)/(lst[1][1].GetEntries()-lst[1][1].GetBinContent(1))
    print "matched m2 but not m1:\t",lst[1][1].GetBinContent(4),lst[1][1].GetBinContent(4)/(lst[1][1].GetEntries()-lst[1][1].GetBinContent(1))
    print 'matched neither m1 nor m2:\t',lst[1][1].GetBinContent(1)

    return lst

if __name__=='__main__':
    a=gethist('/afs/cern.ch/work/y/yanchu/graviton/CMSSW_7_6_3_patch2/src/CMGTools/XZZ2l2nu/cfg/trgmc/dy.root') # hadd all DY tree.root files
#    a=gethist('/afs/cern.ch/work/y/yanchu/graviton/CMSSW_7_6_3_patch2/src/CMGTools/XZZ2l2nu/cfg/trgdt/JetHT_Run2015D_16Dec/vvTreeProducer/tree.root') # hadd all data tree.root files
    c=TCanvas()
    a[0][0].SetTitle('selected electrons trigger matching deltPt distribution;deltPt(GeV)')
    a[1][0].SetTitle('selected muons trigger matching deltPt distribution;deltPt(GeV)')
    a[0][2].SetTitle('selected electrons trigger matching deltR distribution;deltR')
    a[1][2].SetTitle('selected muons trigger matching deltR distribution;deltR')
    a[0][0].Draw()
    c.Print('edPtdy.png')
    a[1][0].Draw()
    c.Print('mdPtdy.png')
    a[0][2].Draw()
    c.Print('edRdy.png')
    a[1][2].Draw()
    c.Print('mdRdy.png')
