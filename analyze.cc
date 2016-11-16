#include<iostream>
#include<fstream>
#include <utility> // std::pair
#include <vector> // std::pair
#include "menulib.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"


struct L1Seed{
     size_t index;
     size_t hits;
     std::string name;
     
     L1Seed(size_t i, size_t h, std::string n)
    : index(i), hits(h), name(n) {}
};

using namespace std;
int main(){
     TFile *myFile = TFile::Open("L1Ntuple.root");
     if (!myFile || myFile->IsZombie()) {
          return 1;
     }

     const size_t max_triggers = 512;
     // form a list of trigger IDs
     ofstream report("report.txt");
     std::vector<L1Seed> L1Seeds;
     for(size_t i=0; i<max_triggers; i++){
          if (getNameFromId(i) != "") {
               L1Seeds.emplace_back(i, 0, getNameFromId(i));
               //report << L1Seeds.size() -1 << ":" << getNameFromId(i);
               
          }
     }
     cout << "Algos in menu: " << L1Seeds.size() << endl;

   
     TTreeReader myReader("l1UpgradeTree/L1UpgradeTree", myFile);
     TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat> data(myReader, "L1Upgrade");
   

     AlgorithmFunction trigger;
     size_t events = 0;
     while (myReader.Next()) {
          ++events;
          for (size_t i=0; i<L1Seeds.size() ; i++){
               trigger=getFuncFromId(L1Seeds[i].index);
               if ( (*trigger)(data) )
                    ++L1Seeds[i].hits; // = L1Seeds[i].hits + 1;
          }
     }
     cout << "Analyzed " << events << " events" << endl;
     //float effPercent;
     for (size_t i=0; i<L1Seeds.size() ; i++){
          //effPercent = 100.0 * (float)L1Seeds[i].hits / (float)events;
          report << L1Seeds[i].index << " : " << L1Seeds[i].name << " : " << L1Seeds[i].hits << endl;
     }
   
     report.close();
     cout << "report written to report.txt\n";
     return 0;
}