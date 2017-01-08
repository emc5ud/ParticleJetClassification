#include <TH1F.h>
#include <stdio.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <fstream>
#include <string>
#include <TProfile2D.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>
#include <stdio.h>
//#include <iterator>
//#include <algorithm>
#include <stdlib.h>

using namespace std;

vector<string> split(const char *str, char c);

void jet_structure(){
  ifstream infile("../j_plus_c_output.txt");
  string line;
  int numJets = 0;
  TProfile2D *  jetShapeLight = new TProfile2D("hprof2dl", "Light Jet Profile; #Delta phi; #Delta eta",30,-.5,.5,40,-.3,.3,0,2000000);
  TProfile2D *  jetShapeTop   = new TProfile2D("hprof2dt", "Top Jet Profile; #Delta phi; #Delta eta",30,-.5,.5,40,-.3,.3,0,2000000);
  TProfile2D *  jetShapeGluon = new TProfile2D("hprof2dg", "Gluon Jet Profile;#Delta phi;#Delta eta",30,-.5,.5,40,-.3,.3,0,2000000);
  jetShapeLight ->  GetXaxis() -> CenterTitle();
  jetShapeTop ->    GetXaxis() -> CenterTitle();
  jetShapeGluon -> GetXaxis() -> CenterTitle();
  jetShapeLight ->  GetYaxis() -> CenterTitle();
  jetShapeTop ->    GetYaxis() -> CenterTitle();
  jetShapeGluon -> GetYaxis() -> CenterTitle();
  jetShapeLight ->  GetXaxis() -> SetTitleSize(.05);
  jetShapeTop ->    GetXaxis() -> SetTitleSize(.05);
  jetShapeGluon -> GetXaxis() -> SetTitleSize(.05);
  jetShapeLight ->  GetYaxis() -> SetTitleSize(.05);
  jetShapeTop ->    GetYaxis() -> SetTitleSize(.05);
  jetShapeGluon -> GetYaxis() -> SetTitleSize(.05);
  //TH2F *  jetShape= new TH2F("hprof2d", "Jet Constituent Profile;phi;eta",30,-.5,.5,40,-.5,.5);
  

  vector<double> pt, eta, phi;
  double j_eta = 0, j_phi = 0;
  int source = -1;
  while(getline(infile, line)){
    if(line[0] == 'J' && numJets > 0){
      for(unsigned int i = 0; i < pt.size(); i++){
	if(abs(source) == 21)
	  jetShapeGluon->Fill(phi[i]-j_phi,eta[i]-j_eta,pt[i]);
	else if(abs(source) == 6)
	  jetShapeTop->Fill(phi[i]-j_phi,eta[i]-j_eta,pt[i]);
	else if(abs(source) != 4)
	  jetShapeLight->Fill(phi[i]-j_phi,eta[i]-j_eta,pt[i]);
      }
      vector<string> toks = split(line.c_str(), ' ');
      source = atoi(toks[2].c_str());
      j_eta = atof(toks[9].c_str());
      j_phi = atof(toks[7].c_str());
      eta.clear();
      pt.clear();
      phi.clear();
      numJets+=1;
    }
    else if(numJets == 0 && line[0] == 'J'){
      numJets+=1;
      vector<string> toks = split(line.c_str(), ' ');
      source = atoi(toks[2].c_str());
      j_eta = atof(toks[9].c_str());
      j_phi = atof(toks[7].c_str());
    }
    if(line[0] == 'C'){
      vector<string> toks = split(line.c_str(), ' ');
      pt.push_back(atof(toks[4].c_str()));
      eta.push_back(atof(toks[7].c_str()));
      phi.push_back(atof(toks[8].c_str()));
    }
  }
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas();
  c1->SetRightMargin(.13);
  jetShapeGluon->Draw("colz");
  TCanvas *c2 = new TCanvas();
  c2->SetRightMargin(.13);
  jetShapeTop->Draw("colz");
  TCanvas *c3 = new TCanvas();
  c3->SetRightMargin(.13);
  jetShapeLight->Draw("colz");
    
  c3->SaveAs("lightjetPZ.png");
  c2->SaveAs("topjetPZ.png");
  c1->SaveAs("gluonjetPZ.png");

}

vector<string> split(const char *str, char c)
{
    vector<string> result;

    do
    {
        const char *begin = str;

        while(*str != c && *str)
            str++;

        result.push_back(string(begin, str));
    } while (0 != *str++);

    return result;
}
