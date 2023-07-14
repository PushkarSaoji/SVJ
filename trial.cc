#include "sift.h"
#include "Pythia8/Pythia.h"
#include "fstream"
using namespace Pythia8;

int main(int argc, char *argv[]) {
  // Number of events, generated and listed ones.
  int nEvent = 1000;
  int seed = argc>1? std::stoi(argv[1]) : 0;
  
  // Generator. LHC process and output selection. Initialization.
  Pythia pythia;
  pythia.readFile("QCD.dat");
  
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + to_string(seed));
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Next:numberShowInfo = 0");
  
  pythia.init();
  
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    vector<pair<Particle,vector<int>>> jet_list = sift(pythia.event);
    if (jet_list.size() < 1) continue;
  
    jet_list.erase(std::remove_if(jet_list.begin(), jet_list.end(), [](const std::pair<Pythia8::Particle, std::vector<int>>& k1) {
    return k1.first.pT() < 10;
}), jet_list.end());
    
    
    sort(jet_list.begin(),jet_list.end(),sortFunction); 
    
    for (const auto& k3: jet_list){
      cout << k3.first.pT() << "  ";
    }
    cout << "\n" << endl;
    
    if (jet_list.size() < 2) continue;
  
    Event event_j=pythia.process;
    for (int i=0; i< 2 ; i++) {                        // for the 2 hardest jets in an event
      event_j.reset();
      if ( jet_list[i].second.size() < 2) continue;
      for (int j:jet_list[i].second ) {    // for all constituents in a jet
        event_j.append(pythia.event[j]);
      }      
      vector<pair<Particle,vector<int>>> subjet_list = sub_sift(event_j);
      
      sort(subjet_list.begin(),subjet_list.end(),sortFunction);
      
      cout << "jet :" << i << " ,  n_subjets : " << subjet_list.size() << endl;
      for (const auto& k3: subjet_list){
        cout << k3.first.pT() << "  ";
      }
      cout << "\n" << endl;
    
    }
  }  
  return 0;
}
  

