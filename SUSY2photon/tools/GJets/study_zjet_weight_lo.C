{
  TFile* fout = new TFile("dyjets_zpt_weight_lo_nlo_sel.root", "recreate");

  gROOT->ProcessLine(".x tdrstyle.C");
  
  TFile* fdyzpt = new TFile("UnfoldingOutputZPt.root");
  TH1D* hdyzptdt = (TH1D*)fdyzpt->Get("hUnfold");
  TH1D* hdyzptmc = (TH1D*)fdyzpt->Get("hTruth");
  hdyzptdt->SetName("hdyzptdt");
  hdyzptmc->SetName("hdyzptmc");
  TH1D* hdyzpt_dtmc_ratio = (TH1D*)hdyzptdt->Clone("hdyzpt_dtmc_ratio");
  hdyzpt_dtmc_ratio->Divide(hdyzptmc);
  TGraphErrors* gdyzpt_dtmc_ratio = new TGraphErrors(hdyzpt_dtmc_ratio);
  gdyzpt_dtmc_ratio->SetName("gdyzpt_dtmc_ratio");

  TF1* fczpt1 = new TF1("fczpt1", "([0]-[1]*TMath::Erf((x-[2])/[3])+[4]*TMath::Erf((x-[5])/[6]))",0,3000);

  fczpt1->SetParameters(1.02852,0.0949640,19.0422,10.4487,0.0758834,56.1146,41.1653);
 
//(1.02852 - 0.0949640*TMath::Erf((gen_ptll-19.0422)/10.4487) + 0.0758834*Tmath::Erf((gen_ptll- 56.1146)/41.1653)) 

  TH1D* hdyzpt_dtmc_ratio_smooth= (TH1D*)hdyzpt_dtmc_ratio->Clone("hdyzpt_dtmc_ratio_smooth");
  hdyzpt_dtmc_ratio_smooth->Smooth();
  hdyzpt_dtmc_ratio_smooth->SetLineColor(4);
  hdyzpt_dtmc_ratio_smooth->SetMarkerColor(4);

  TF1* fczpt2 = new TF1("fczpt2", "[0]-[1]*TMath::Erf((x-[2])/[3])+[4]*TMath::Erf((x-[5])/[6])-[7]*TMath::Erf((x-[8])/[9])",0,3000);;
  fczpt2->SetParameters(1.02852,0.0949640,19.0422,10.4487,0.0758834,56.1146,41.1653, 0.1, 100, 100);
  //fczpt2->SetParameters(-2.978752, 0.131293, 14.500684, 10.332572, 4.097619, -252.264980, 196.546957, 0.258976, 142.987571, 394.031326);

//Par   0                   p0 =   -2.978752 
//Par   1                   p1 =    0.131293 
//Par   2                   p2 =   14.500684 
//Par   3                   p3 =   10.332572 
//Par   4                   p4 =    4.097619 
//Par   5                   p5 =  -252.264980 
//Par   6                   p6 =  196.546957 
//Par   7                   p7 =    0.258976 
//Par   8                   p8 =  142.987571 
//Par   9                   p9 =  394.031326 

  Int_t nfit=3;
  for (int i=0; i<nfit; i++){
    hdyzpt_dtmc_ratio->Fit(fczpt2);
  }

  // print function
  std::cout << "  fczpt2->SetParameters(" ;
  for (int i=0; i<fczpt2->GetNpar(); i++){
    if (i==0) std::cout << fczpt2->GetParameter(i);
    else std::cout << "," << fczpt2->GetParameter(i);
  }
  std::cout << ");" << std::endl;

  char name[3000];
  //sprintf(name, "  double weight = %f-%f*TMath::Erf((x-%f)/%f)+%f*TMath::Erf((x-%f)/%f)-%f*TMath::Erf((x-%f)/%f)",
  sprintf(name, "  double weight = (%f-%f*TMath::Erf((genZ_pt-%f)/%f)+%f*TMath::Erf((genZ_pt-%f)/%f)-%f*TMath::Erf((genZ_pt-%f)/%f))", 
     fczpt2->GetParameter(0),
     fczpt2->GetParameter(1),
     fczpt2->GetParameter(2),
     fczpt2->GetParameter(3),
     fczpt2->GetParameter(4),
     fczpt2->GetParameter(5),
     fczpt2->GetParameter(6),
     fczpt2->GetParameter(7),
     fczpt2->GetParameter(8),
     fczpt2->GetParameter(9)
  );
  std::cout << name << std::endl;
  
  //TFile* flotonlo = TFile::Open("study_zpt_new2.root");
  //TFile* flotonlo = TFile::Open("study_zpt_new3.root");
  TFile* flotonlo = TFile::Open("study_zpt_old.root");
  TH1D* hzptnlo = (TH1D*)flotonlo->Get("hzpt1");
  TH1D* hzptlo = (TH1D*)flotonlo->Get("hzpt2");
  TH1D* hzptrnlolo = (TH1D*)flotonlo->Get("hzptr12");
  hzptnlo->SetName("hzptnlo");
  hzptlo->SetName("hzptlo");
  hzptrnlolo->SetName("hzptrnlolo");

  TH1D* hzpt3knlo = (TH1D*)flotonlo->Get("hzpt3k1");
  TH1D* hzpt3klo = (TH1D*)flotonlo->Get("hzpt3k2");
  //TH1D* hzpt3krnlolo = (TH1D*)flotonlo->Get("hzpt3kr12");
  TH1D* hzpt3krnlolo = (TH1D*)flotonlo->Get("hzptr12");
  hzpt3knlo->SetName("hzpt3knlo");
  hzpt3klo->SetName("hzpt3klo");
  hzpt3krnlolo->SetName("hzpt3krnlolo");

  TH1D* hdyzpt_dtmc_lo_ratio = (TH1D*)hdyzpt_dtmc_ratio->Clone("hdyzpt_dtmc_lo_ratio");
  hdyzpt_dtmc_lo_ratio->Multiply(hzptrnlolo);
  hdyzpt_dtmc_lo_ratio->Draw();
 
  TF1* fczpt3 = new TF1("fczpt3", "[0]-[1]*TMath::Erf((x-[2])/[3])+[4]*TMath::Erf((x-[5])/[6])-[7]*TMath::Erf((x-[8])/[9])",0,3000);
  fczpt3->SetParameters(-0.388696,0.128165,14.5515,10.1973,2.0109,-181.141,184.501,0.783707,-275.078,652.422); 
  //hdyzpt_dtmc_lo_ratio->Fit(fczpt3);
  
  //Int_t nfit=2;
  nfit=1;
  for (int i=0; i<nfit; i++){
    hdyzpt_dtmc_lo_ratio->Fit(fczpt3);
  }

  // print function
  std::cout << "  fczpt3->SetParameters(" ;
  for (int i=0; i<fczpt3->GetNpar(); i++){
    if (i==0) std::cout << fczpt3->GetParameter(i);
    else std::cout << "," << fczpt3->GetParameter(i);
  }
  std::cout << ");" << std::endl;

  //char name[3000];
  //sprintf(name, "  double weight = %f-%f*TMath::Erf((x-%f)/%f)+%f*TMath::Erf((x-%f)/%f)-%f*TMath::Erf((x-%f)/%f)",
  sprintf(name, "  double weight = (%f-%f*TMath::Erf((genZ_pt-%f)/%f)+%f*TMath::Erf((genZ_pt-%f)/%f)-%f*TMath::Erf((genZ_pt-%f)/%f))", 
     fczpt3->GetParameter(0),
     fczpt3->GetParameter(1),
     fczpt3->GetParameter(2),
     fczpt3->GetParameter(3),
     fczpt3->GetParameter(4),
     fczpt3->GetParameter(5),
     fczpt3->GetParameter(6),
     fczpt3->GetParameter(7),
     fczpt3->GetParameter(8),
     fczpt3->GetParameter(9)
  );
  std::cout << name << std::endl;


  TH1D* hdyzpt_dtmc_lo_ratio_smooth = (TH1D*)hdyzpt_dtmc_lo_ratio->Clone("hdyzpt_dtmc_lo_ratio_smooth");
  hdyzpt_dtmc_lo_ratio_smooth->Smooth();


  // fit mc nlo/lo
  //TF1* fczpt4 = new TF1("fczpt4", "[0]-[1]*TMath::Erf((x-[2])/[3])+[4]*TMath::Erf((x-[5])/[6])-[7]*TMath::Erf((x-[8])/[9])",0,3000);
  //fczpt4->SetParameters(-0.388696,0.128165,14.5515,10.1973,2.0109,-181.141,184.501,0.783707,-275.078,652.422);

  hzpt3krnlolo->SetLineColor(1);
  hzpt3krnlolo->SetMarkerColor(1);
  hzpt3krnlolo->SetMarkerStyle(20);
  hzpt3krnlolo->GetXaxis()->SetTitle("gen. p_{T}(Z) (GeV)");
  hzpt3krnlolo->GetYaxis()->SetTitle("NLO/LO ratio");

  hzpt3krnlolo->GetXaxis()->SetTitleOffset(1.25);
  hzpt3krnlolo->GetYaxis()->SetTitleOffset(1.25);
  hzpt3krnlolo->GetXaxis()->SetTitleSize(0.05);
  hzpt3krnlolo->GetYaxis()->SetTitleSize(0.05);
  hzpt3krnlolo->GetXaxis()->SetLabelSize(0.04);
  hzpt3krnlolo->GetYaxis()->SetLabelSize(0.04); 
 
  TF1* fczpt4 = new TF1("fczpt4", "[0]-[1]*TMath::Erf((x-[2])/[3])+[4]*TMath::Erf((x-[5])/[6])-[7]*TMath::Erf((x-[8])/[9])+[10]*TMath::Erf((x-[11])/[12])-[13]*TMath::Erf((x-[14])/[15])+[16]*TMath::Erf((x-[17])/[18])",0,3000);
  fczpt4->SetLineWidth(2);
  fczpt4->SetNpx(10000);

  //Double_t params[] = {1.02852,   0.1,10,10,  1,20,50, 1,60,50, 1,100,50,  1,300,100 };
  //Double_t params[] = {0.960208,0.157272,9.33256,34.3082,0.0757899,20,6.1308,0.17918,125.932,63.9047,1.55286,180.682,228.41,1.0,300,177.811,0.1,400,10};
  //Double_t params[] = {0.975835,0.152064,8.0001,36.3054,0.0746964,20,6.08109,0.0902824,118.121,43.2299,1.12978,174.853,241.396,1.20754,300.073,173.306,0.488323,300,85.7321};
  //Double_t params[] = {0.264091,0.0844949,12.9514,21.9499,0.0870793,20.5089,7.43006,-3.20637,91.0083,481.45,-0.869759,179.36,210.235,1.38243,281.541,551.735,0,0,0};
  Double_t params[] = {0.264091,0.0844949,12.9514,21.9499,0.0870793,20.5089,7.43006,-3.20637,91.0083,481.45,-0.869759,179.36,210.235,0.3,300,50,0,0,0};
  fczpt4->SetParameters(params);

  fczpt4->SetParLimits(2,8,16);
  fczpt4->SetParLimits(5,20,40);
  fczpt4->SetParLimits(8,80,140);
  fczpt4->SetParLimits(11,170,250);
  fczpt4->SetParLimits(14,280,350);
  //fczpt4->SetParLimits(17,300,600);

  fczpt4->FixParameter(16,0);
  fczpt4->FixParameter(17,0);
  fczpt4->FixParameter(18,0);
  
  //Int_t nfit=2;
  nfit=100;
  for (int i=0; i<nfit; i++){
    hzpt3krnlolo->Fit(fczpt4);
  }

  // print function
  std::cout << "  fczpt4->SetParameters(" ;
  for (int i=0; i<fczpt4->GetNpar(); i++){
    if (i==0) std::cout << fczpt4->GetParameter(i);
    else std::cout << "," << fczpt4->GetParameter(i);
  }
  std::cout << ");" << std::endl;

  //char name[3000];
  sprintf(name, "  double weight = (%f-%f*TMath::Erf((genZ_pt-%f)/%f)+%f*TMath::Erf((genZ_pt-%f)/%f)-%f*TMath::Erf((genZ_pt-%f)/%f)+%f*TMath::Erf((genZ_pt-%f)/%f)-%f*TMath::Erf((genZ_pt-%f)/%f)+%f*TMath::Erf((genZ_pt-%f)/%f))",
     fczpt4->GetParameter(0),
     fczpt4->GetParameter(1),
     fczpt4->GetParameter(2),
     fczpt4->GetParameter(3),
     fczpt4->GetParameter(4),
     fczpt4->GetParameter(5),
     fczpt4->GetParameter(6),
     fczpt4->GetParameter(7),
     fczpt4->GetParameter(8),
     fczpt4->GetParameter(9),
     fczpt4->GetParameter(10),
     fczpt4->GetParameter(11),
     fczpt4->GetParameter(12),
     fczpt4->GetParameter(13),
     fczpt4->GetParameter(14),
     fczpt4->GetParameter(15),
     fczpt4->GetParameter(16),
     fczpt4->GetParameter(17),
     fczpt4->GetParameter(18)
  );
  std::cout << name << std::endl;

  fout->cd();
  hdyzptdt->Write();
  hdyzptmc->Write();
  hdyzpt_dtmc_ratio->Write(); 
  hdyzpt_dtmc_ratio_smooth->Write(); 
  gdyzpt_dtmc_ratio->Write(); 
  fczpt1->Write(); 
  fczpt2->Write("fcdyzpt_dtmc_ratio"); 
 
  hzptnlo->Write();
  hzptlo->Write();
  hzptrnlolo->Write();
  hzpt3knlo->Write();
  hzpt3klo->Write();
  hzpt3krnlolo->Write("hdyzpt_mc_nlo_lo_ratio");

  hdyzpt_dtmc_lo_ratio->Write();
  hdyzpt_dtmc_lo_ratio_smooth->Write();
  fczpt3->Write("fcdyzpt_dtmc_lo_ratio"); 
  fczpt4->Write("fcdyzpt_mc_nlo_lo_ratio"); 
}
