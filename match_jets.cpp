#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <string>
//#include <iterator>
//#include <algorithm>
#include <stdlib.h>
#include "fastjet/ClusterSequence.hh"


using namespace fastjet;
using namespace std;

// simple particle object
struct particle {
  int pid;
  int status;
  vector<double> fourVect; //assumed to have 4 elements
  int nextVert;
  particle(int id, int s, vector<double> v, int vid) : pid(id), status(s),
    fourVect(v), nextVert(vid)
  {
  }
  void print(){
    cout << "ID: " << pid << "\t4Vec: (" << fourVect[0] << ", " << fourVect[1] << ", " <<
      fourVect[2] << ", " << fourVect[3] << ")" <<  endl;
  }
};


// simple vertex object
struct vertex {
  std::vector <particle> particles;
  void addPart(particle p){
    particles.push_back(p);
  }
  int vid;
  int nPart;
  vertex(int id, int n) : vid(id), nPart(n)
  {
  }
  void print(){
    cout << "Vertex ID: " << vid << endl;
    for(int i = 0; i < particles.size(); i++)
      
      particles[i].print();
  }
};


double r_squared(const PseudoJet & j, const particle & parton);
vector<string> split( std::string const& original, char separator );
void readEvent(ifstream &infile, ofstream &outfile);
void writeOutJet(const PseudoJet & jet, int jid ,int pid, ofstream & outfile);
void banner(ofstream &outputFile);


int main(int argc, char *argv[]){

  if(argc != 3){
    cout<<"usage: "<< argv[0] <<" <input filename> " << "<output filename>\n";
    return 1;
  }
  ifstream infile(argv[1]);
  ofstream outfile(argv[2]);
  int events = 0;
  

  string line;
  banner(outfile);
  while(getline(infile, line)){
    if(line[0] == 'C'){
      outfile << "E " << events << endl; 
      readEvent(infile, outfile); // reads to end of event
      events++;
      //if(events == 1000) break;
    }
    //cout << tokens[0] << endl;
  }
  cout << events << endl;
}

void readEvent(ifstream &infile, ofstream &outfile){
  vector<particle> finalState;
  vector<particle> partons;
  vector<PseudoJet> particles;
  string line;
  while(getline(infile,line)){
    if(line[0] == 'P'){ 
      vector<string> toks = split(line, ' ');
      vector<double> fourVect;
      fourVect.push_back(atof(toks[3].c_str()));
      fourVect.push_back(atof(toks[4].c_str()));
      fourVect.push_back(atof(toks[5].c_str()));
      fourVect.push_back(atof(toks[6].c_str()));

      //storing pid, status, 4vector, and id of next vertex (not used)
      particle p(atoi(toks[2].c_str()), atoi(toks[8].c_str()),
		 fourVect, abs(atoi(toks[11].c_str())));
      if(p.status == 1){
	finalState.push_back(p);
	particles.push_back(fourVect);
      }
      //the first four particles have the good stuff
      if((p.status == 23 || p.status == 22) && partons.size() < 4){ 
	partons.push_back(p);
      }
    }
    if(line[0] == 'E'){
      //cout << finalState.size() << endl;
      // signifies end of event
      // now we create our jets 
      double R = 0.7;
      int num_const=0;
      JetDefinition jet_def(antikt_algorithm, R);
      // run the clustering, extract the jets    
      ClusterSequence cs(particles, jet_def);
      vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
      //matching jets
      for(int j_idx = 0; j_idx < jets.size(); j_idx++){
	PseudoJet jet = jets[j_idx];
	int nMatch = 0;
	int matchID; //pid of matching particle
	for(int p = 0; p < 4 && p < partons.size(); p++){
	  if(r_squared(jet, partons[p]) < .3){
	    matchID = partons[p].pid;
	    writeOutJet(jet,j_idx, partons[p].pid, outfile);
	  }
	  
	}
      }
      
      return;
    }
  }
}



double r_squared(const PseudoJet & j, const particle & parton){
  PseudoJet p(parton.fourVect);
  return (pow(j.rap()-p.rap(),2) + pow(j.phi()-p.phi(),2));
}

vector<string> split( std::string const& original, char separator )
{
    std::vector<std::string> results;
    std::string::const_iterator start = original.begin();
    std::string::const_iterator end = original.end();
    std::string::const_iterator next = std::find( start, end, separator );
    while ( next != end ) {
        results.push_back( std::string( start, next ) );
        start = next + 1;
        next = std::find( start, end, separator );
    }
    results.push_back( std::string( start, next ) );
    return results;
}


//Writing out jet and constituent particle info to file
void writeOutJet(const PseudoJet & jet, int jid ,int pid, ofstream & outfile){
  string s = " ";
  outfile << "J " << jid << s << pid << s
	  << jet.E() << s << jet.px() << s << jet.py() << s << jet.pz()
	  << s << jet.phi() << s << jet.rap() << s << jet.eta() << s
	  << jet.pt2() << s << jet.pt() <<  s
	  << jet.m2() << s  << jet.m() << s
	  << jet.mt2() << s << jet.mt()  << endl;
  
  vector<PseudoJet> cs = jet.constituents();
  for (unsigned j = 0; j < cs.size(); j++) {
    outfile << "C "
	    << cs[j].E() << s << cs[j].px() << s << cs[j].py() << s <<cs[j].pz()
	    << s << cs[j].pt() << s << cs[j].rap() << s << cs[j].eta() << s
	    << cs[j].phi() << s  << endl;
  }
  return;
}

void banner(ofstream &outputFile){
  outputFile << "##############################################" << endl;
  outputFile << "# This is the Output from Fastjet Clustering #" << endl;
  outputFile << "# The File contains all of the information   #" << endl;
  outputFile << "# from the fastjet jet class. You can ammend #" << endl;
  outputFile << "# the amount of information dumped if you    #" << endl;
  outputFile << "# want less. Fastjet Clustering output was   #" << endl;
  outputFile << "# matched to initial partons by trajectory.  #" << endl;
  outputFile << "#                                            #" << endl;
  outputFile << "#            - Chris Wiens. cdw9bf           #" << endl;
  outputFile << "#            - Eric Culbertson. emc5ud       #" << endl;
  outputFile << "##############################################" << endl;
  outputFile << "#                 FORMAT                     #" << endl;
  outputFile << "# Rows with J are jets, C are particles, E   #" << endl;
  outputFile << "# are new events with event number following #" << endl;
  outputFile << "##############################################" << endl;
  outputFile << "# Jet Column Headers {only for rows with J at the beginning}" << endl;
  outputFile << "# number sourcePID E px py pz phi rap eta pt2 pt m2 m mt2 mt" << endl;
  outputFile << "# Particle Column Headers {only for rows with C}" << endl;
  outputFile << "# E px py pz pt rap eta phi" << endl;
  outputFile << "#" << endl;
  outputFile << "# For columns with P, xx xx xx xx are the header values" << endl;

  return;
}
