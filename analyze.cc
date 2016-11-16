#include<iostream>
#include "menulib.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

using namespace std;
int main(){
     TFile *myFile = TFile::Open("L1Ntuple.root");
   if (!myFile || myFile->IsZombie()) {
      return 1;
   }
   
   TTreeReader myReader("l1UpgradeTree", myFile);
   TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat> data(myReader, "L1UpgradeTree");
   
   while (myReader.Next()) {
    bool a = L1_AlwaysTrue(data); 
    cout << a << endl;
   }
   
   return 0;
}