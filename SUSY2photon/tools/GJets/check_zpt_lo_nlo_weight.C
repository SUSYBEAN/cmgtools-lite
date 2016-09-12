{


TFile* file1 = TFile::Open("/home/heli/XZZ/80X_20160818_light_Skim/DYJetsToLL_M50_NoRecoil.root");
TFile* file2 = TFile::Open("/home/heli/XZZ/80X_20160818_light_Skim/DYJetsToLL_M50_MGMLM_Ext1_NoRecoil.root");

std::string tag = "NoRecoil";

std::string wt1 = "ZPtWeight";
std::string wt2 = "ZPtWeight";

//std::string wt1 = "(1)";
//std::string wt2 = "(0.434409-0.083122*TMath::Erf((genZ_pt-13.081120)/21.807790)+0.087146*TMath::Erf((genZ_pt-20.509608)/7.431639)--2.030941*TMath::Erf((genZ_pt-81.143629)/380.082763)+-0.778689*TMath::Erf((genZ_pt-170.436870)/195.281851)-0.657784*TMath::Erf((genZ_pt-347.533486)/521.639701))";
//std::string wt2 = "(0.791494-0.141873*TMath::Erf((genZ_pt-12.559276)/27.351067)+0.087455*TMath::Erf((genZ_pt-20.607234)/7.472463)-0.177693*TMath::Erf((genZ_pt-85.639301)/60.038844)+1.874484*TMath::Erf((genZ_pt-170.173370)/267.664379)-2.441469*TMath::Erf((genZ_pt-285.145333)/148.617628)+1.228828*TMath::Erf((genZ_pt-304.776950)/110.678623))";


//std::string wt1 = "(1)";
//std::string wt2 = "(1.037341-0.153197*TMath::Erf((genZ_pt-8.000000)/38.211660)+0.073979*TMath::Erf((genZ_pt-20.000000)/6.046473)-0.137818*TMath::Erf((genZ_pt-117.872834)/47.944872)+1.041567*TMath::Erf((genZ_pt-183.745762)/203.876399)-1.196846*TMath::Erf((genZ_pt-280.000000)/123.753850)+0.531018*TMath::Erf((genZ_pt-300.000000)/87.230449))";


//std::string wt1 = "(-10.139730-0.174889*TMath::Erf((genZ_pt-13.654379)/11.877161)+11.172971*TMath::Erf((genZ_pt--173.890605)/107.479445)-0.113491*TMath::Erf((genZ_pt-315.120674)/103.326396))";
//std::string wt2 = wt1+"*(1.037341-0.153197*TMath::Erf((genZ_pt-8.000000)/38.211660)+0.073979*TMath::Erf((genZ_pt-20.000000)/6.046473)-0.137818*TMath::Erf((genZ_pt-117.872834)/47.944872)+1.041567*TMath::Erf((genZ_pt-183.745762)/203.876399)-1.196846*TMath::Erf((genZ_pt-280.000000)/123.753850)+0.531018*TMath::Erf((genZ_pt-300.000000)/87.230449))";

//std::string wt1 = "(-0.389014-0.128125*TMath::Erf((genZ_pt-14.548919)/10.195989)+2.013448*TMath::Erf((genZ_pt--181.075481)/184.307672)-0.782235*TMath::Erf((genZ_pt--274.872052)/645.686104))";
//std::string wt2 = "(0.690238-0.068081*TMath::Erf((genZ_pt-10.254618)/5.074902)+3.460980*TMath::Erf((genZ_pt-151.407304)/235.391295)-3.168676*TMath::Erf((genZ_pt-184.684699)/214.128678))";
//std::string wt2 = "(-0.115221-0.074336*TMath::Erf((genZ_pt-10.262228)/5.500722)+2.239024*TMath::Erf((genZ_pt--128.514391)/279.350304)-1.211303*TMath::Erf((genZ_pt--7.540941)/381.585801))";

std::string sel1 = "(1)*(genWeight*"+wt1+")";
std::string sel2 = "(1)*(genWeight*"+wt2+")";

//std::string sel1 = "(1)*(puWeight68075*genWeight*"+wt1+")";
//std::string sel2 = "(1)*(puWeight68075*genWeight*"+wt2+")";

TTree* tree1 = (TTree*)file1->Get("tree");
TTree* tree2 = (TTree*)file2->Get("tree");

//Double_t ZPtBins[] = {0,1.25,2.5,3.75,5,6.25,7.5,8.75,10,11.25,12.5,15,17.5,20,25,30,35,40,45,50,60,70,80,90,100,110,130,150,170,190,220,250,300,500,1000};
//Double_t ZPtBins[] = {0, 1.25, 2.5, 3.75, 5, 6.25, 7.5, 8.75, 10, 11.25, 12.5, 15, 17.5, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 130, 150, 170, 190, 220, 250, 400, 1000};
Double_t ZPtBins[] = {0, 1.25, 2.5, 3.75, 5, 6.25, 7.5, 8.75, 10, 11.25, 12.5, 15, 17.5, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110, 130, 150, 170, 190, 220, 250, 400, 3000};
Int_t NZPtBins = sizeof(ZPtBins)/sizeof(ZPtBins[0]) - 1;

char name[1000];

TH1D* hzpt1 = new TH1D("hzpt1", "hzpg1", NZPtBins, ZPtBins);
TH1D* hzpt2 = new TH1D("hzpt2", "hzpg2", NZPtBins, ZPtBins);

//TH1D* hzpt1 = new TH1D("hzpt1", "hzpg1", 100, 0, 1500);
//TH1D* hzpt2 = new TH1D("hzpt2", "hzpg2", 100, 0, 1500);

hzpt1->Sumw2();
hzpt2->Sumw2();


tree1->Draw("llnunu_l1_pt>>hzpt1", sel1.c_str());
tree2->Draw("llnunu_l1_pt>>hzpt2", sel2.c_str());

//tree1->Draw("genZ_pt>>hzpt1", sel1.c_str());
//tree2->Draw("genZ_pt>>hzpt2", sel2.c_str());

hzpt1->SetLineColor(2);
hzpt2->SetLineColor(4);

hzpt1->Scale(1./hzpt1->Integral("width"));
hzpt2->Scale(1./hzpt2->Integral("width"));

hzpt1->Draw();
hzpt2->Draw("same");


TH1D* hzptr12 = (TH1D*)hzpt1->Clone("hzptr12");
hzptr12->Divide(hzpt2);

hzptr12->Draw();

/*
TH1D* hzeta1 = new TH1D("hzeta1", "hzeta1", 100, -10, 10);
TH1D* hzeta2 = new TH1D("hzeta2", "hzeta2", 100, -10, 10);

hzeta1->Sumw2();
hzeta2->Sumw2();


tree1->Draw("llnunu_l1_eta>>hzeta1", sel1.c_str());
tree2->Draw("llnunu_l1_eta>>hzeta2", sel2.c_str());

hzeta1->SetLineColor(2);
hzeta2->SetLineColor(4);


hzeta1->Scale(1./hzeta1->Integral());
hzeta2->Scale(1./hzeta2->Integral());

hzeta1->Draw();
hzeta2->Draw("same");



TH1D* hzmass1 = new TH1D("hzmass1", "hzmass1", 100, 50, 180);
TH1D* hzmass2 = new TH1D("hzmass2", "hzmass2", 100, 50, 180);

hzmass1->Sumw2();
hzmass2->Sumw2();


tree1->Draw("llnunu_l1_mass>>hzmass1", sel1.c_str());
tree2->Draw("llnunu_l1_mass>>hzmass2", sel2.c_str());

hzmass1->SetLineColor(2);
hzmass2->SetLineColor(4);


hzmass1->Scale(1./hzmass1->Integral());
hzmass2->Scale(1./hzmass2->Integral());

hzmass1->Draw();
hzmass2->Draw("same");



//TH1D* hzeta1_wt = new TH1D("hzeta1_wt", "hzeta1_wt", 100, -10, 10);
//TH1D* hzeta2_wt = new TH1D("hzeta2_wt", "hzeta2_wt", 100, -10, 10);

//hzeta1_wt->Sumw2();
//hzeta2_wt->Sumw2();
*/
sprintf(name, "check_zpt_lo_nlo_weight_%s.root", tag.c_str());
TFile* fout = new TFile(name, "recreate");

hzpt1->Write();
hzpt2->Write();
hzptr12->Write();
/*
hzeta1->Write();
hzeta2->Write();
hzmass1->Write();
hzmass2->Write();
*/
fout->Close();


}
