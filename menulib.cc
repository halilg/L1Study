/* automatically generated from L1Menu_Collisions2016_v6r8 with menu2lib.py */
/* https://gitlab.cern.ch/cms-l1t-utm/scripts */

#include <algorithm>
#include <map>
#include <string>

#include "menulib.hh"

// temporary hack
namespace L1Analysis
{
  int kMinBiasHFP0 = 11;
  int kMinBiasHFM0 = 12;
  int kMinBiasHFP1 = 13;
  int kMinBiasHFM1 = 14;
  int kMissingEtHF = 8;
  int kTotalEtEm = 16;
}


// utility methods
void
getCombination(int N,
               int K,
               std::vector<std::vector<int> >& combination)
{
  std::string bitmask(K, 1);
  bitmask.resize(N, 0);

  do
  {
    std::vector<int> set;
    for (int ii = 0; ii < N; ++ii)
    {
      if (bitmask[ii]) set.push_back(ii);
    }
    combination.push_back(set);
  }
  while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}


void
getPermutation(int N,
               std::vector<std::vector<int> >& permutation)
{
  std::vector<int> indicies(N);
  for (int ii = 0; ii < N; ii++) indicies.at(ii) = ii;

  do
  {
    std::vector<int> set;
    for (int ii = 0; ii < N; ++ii)
    {
      set.push_back(indicies.at(ii));
    }
    permutation.push_back(set);
  }
  while (std::next_permutation(indicies.begin(), indicies.end()));
}





// generate conditions
      


        
                            

bool
CaloCaloCorrelation_14433217633607694784
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    bool etaWindow1;
                              // EG18: ET >= 36 at BX = 0
      if (not (data->egIEt.at(ii) >= 36)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(ii)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
          if (not etaWindow1) continue;
    for (size_t jj = 0; jj < data->tauBx.size(); jj++)
    {
      if (not (data->tauBx.at(jj) == 0)) continue;
                                // TAU24: ET >= 48 at BX = 0
      if (not (data->tauIEt.at(jj) >= 48)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(jj)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(jj)) and (data->tauIEta.at(jj) <= 48));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
          // 0.2 <= DeltaEta <= 10.0
      int iEta = data->egIEta.at(ii);
    unsigned int deltaIEta = abs(iEta - data->tauIEta.at(jj));
  unsigned int deltaEta = LUT_DETA_EG_TAU[deltaIEta];

    minimum = (long long)(0.2 * POW10[3]);
    maximum = (long long)(10.0 * POW10[3]);
    if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


        


        
                            

bool
CaloCaloCorrelation_14500771630165735872
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    bool etaWindow1;
                              // EG20: ET >= 40 at BX = 0
      if (not (data->egIEt.at(ii) >= 40)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(ii)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
          if (not etaWindow1) continue;
    for (size_t jj = 0; jj < data->tauBx.size(); jj++)
    {
      if (not (data->tauBx.at(jj) == 0)) continue;
                                // TAU25: ET >= 50 at BX = 0
      if (not (data->tauIEt.at(jj) >= 50)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(jj)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(jj)) and (data->tauIEta.at(jj) <= 48));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
          // 0.2 <= DeltaEta <= 10.0
      int iEta = data->egIEta.at(ii);
    unsigned int deltaIEta = abs(iEta - data->tauIEta.at(jj));
  unsigned int deltaEta = LUT_DETA_EG_TAU[deltaIEta];

    minimum = (long long)(0.2 * POW10[3]);
    maximum = (long long)(10.0 * POW10[3]);
    if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


        


        
                            

bool
CaloCaloCorrelation_14501897532220062144
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    bool etaWindow1;
                              // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(ii) >= 44)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(ii)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
          if (not etaWindow1) continue;
    for (size_t jj = 0; jj < data->tauBx.size(); jj++)
    {
      if (not (data->tauBx.at(jj) == 0)) continue;
                                // TAU26: ET >= 52 at BX = 0
      if (not (data->tauIEt.at(jj) >= 52)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(jj)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(jj)) and (data->tauIEta.at(jj) <= 48));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
          // 0.2 <= DeltaEta <= 10.0
      int iEta = data->egIEta.at(ii);
    unsigned int deltaIEta = abs(iEta - data->tauIEta.at(jj));
  unsigned int deltaEta = LUT_DETA_EG_TAU[deltaIEta];

    minimum = (long long)(0.2 * POW10[3]);
    maximum = (long long)(10.0 * POW10[3]);
    if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


        


        
                            

bool
CaloCaloCorrelation_2085349615869404806
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    bool etaWindow1;
                              // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(ii) >= 44)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(ii)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
          if (not etaWindow1) continue;
    for (size_t jj = 0; jj < data->tauBx.size(); jj++)
    {
      if (not (data->tauBx.at(jj) == 0)) continue;
                                // TAU20: ET >= 40 at BX = 0
      if (not (data->tauIEt.at(jj) >= 40)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(jj)) and (data->tauIEta.at(jj) <= 48));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
          // 0.2 <= DeltaEta <= 10.0
      int iEta = data->egIEta.at(ii);
    unsigned int deltaIEta = abs(iEta - data->tauIEta.at(jj));
  unsigned int deltaEta = LUT_DETA_EG_TAU[deltaIEta];

    minimum = (long long)(0.2 * POW10[3]);
    maximum = (long long)(10.0 * POW10[3]);
    if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


        


          

  bool
CaloEsumCorrelation_16768129600233686289
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
                              // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(ii) >= 120)) continue;

          

    for (size_t jj = 0; jj < data->sumBx.size(); jj++)
    {
      if (not (data->sumType.at(jj) == L1Analysis::kMissingEt)) continue;
      if (not (data->sumBx.at(jj) == 0)) continue;
                          // ETM75: ET >= 75.0 at BX = 0
      if (not (data->sumEt.at(jj) >= 75.0)) continue;
          
                  int iPhi = data->sumIPhi.at(jj);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(ii));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  unsigned int deltaPhi = LUT_DPHI_JET_ETM[deltaIPhi];

    long long minimum = (long long)(0.4 * POW10[3]);
    long long maximum = (long long)(3.15 * POW10[3]);
    if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;
    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        



      

bool
CaloMuonCorrelation_10674670645420326056
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    
                              // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(ii) >= 20)) continue;

          
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(jj) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
          // 1.0 <= DeltaPhi <= 3.15
      int iPhi = data->egIPhi.at(ii);
      iPhi = LUT_PHI_EG2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhi.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  unsigned int deltaPhi = LUT_DPHI_EG_MU[deltaIPhi];

    minimum = (long long)(1.0 * POW10[3]);
    maximum = (long long)(3.15 * POW10[3]);
    if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

        



      

bool
CaloMuonCorrelation_10791898730651162912
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    
                              // JET32: ET >= 64 at BX = 0
      if (not (data->jetIEt.at(ii) >= 64)) continue;

          
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(jj) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
          // 0.0 <= DeltaPhi <= 0.4
      int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhi.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];

    minimum = (long long)(0.0 * POW10[3]);
    maximum = (long long)(0.4 * POW10[3]);
    if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

        



      

bool
CaloMuonCorrelation_15993852978349077723
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    bool etaWindow1;
                              // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(ii) >= 240)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 68));
            
          if (not etaWindow1) continue;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
          // 0.0 <= DeltaPhi <= 0.4
      int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhi.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];

    minimum = (long long)(0.0 * POW10[3]);
    maximum = (long long)(0.4 * POW10[3]);
    if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;

          // 0.0 <= DeltaEta <= 0.4
      int iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    unsigned int deltaIEta = abs(iEta - data->muonIEta.at(jj));
  unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];

    minimum = (long long)(0.0 * POW10[3]);
    maximum = (long long)(0.4 * POW10[3]);
    if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

        



      

bool
CaloMuonCorrelation_16240387826857744385
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    bool etaWindow1;
                              // JET16: ET >= 32 at BX = 0
      if (not (data->jetIEt.at(ii) >= 32)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 68));
            
          if (not etaWindow1) continue;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
          // 0.0 <= DeltaPhi <= 0.4
      int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhi.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];

    minimum = (long long)(0.0 * POW10[3]);
    maximum = (long long)(0.4 * POW10[3]);
    if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;

          // 0.0 <= DeltaEta <= 0.4
      int iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    unsigned int deltaIEta = abs(iEta - data->muonIEta.at(jj));
  unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];

    minimum = (long long)(0.0 * POW10[3]);
    maximum = (long long)(0.4 * POW10[3]);
    if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

        



      

bool
CaloMuonCorrelation_16240389188362377217
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    bool etaWindow1;
                              // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(ii) >= 120)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 68));
            
          if (not etaWindow1) continue;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
          // 0.0 <= DeltaPhi <= 0.4
      int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhi.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];

    minimum = (long long)(0.0 * POW10[3]);
    maximum = (long long)(0.4 * POW10[3]);
    if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;

          // 0.0 <= DeltaEta <= 0.4
      int iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    unsigned int deltaIEta = abs(iEta - data->muonIEta.at(jj));
  unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];

    minimum = (long long)(0.0 * POW10[3]);
    maximum = (long long)(0.4 * POW10[3]);
    if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

        
bool
DoubleEG_14367260113818400607
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleEG_14367282104050956127
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG15: ET >= 30 at BX = 0
      if (not (data->egIEt.at(idx) >= 30)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleEG_14367290900143979231
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG17: ET >= 34 at BX = 0
      if (not (data->egIEt.at(idx) >= 34)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG17: ET >= 34 at BX = 0
      if (not (data->egIEt.at(idx) >= 34)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleEG_14367295298190490335
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG18: ET >= 36 at BX = 0
      if (not (data->egIEt.at(idx) >= 36)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG17: ET >= 34 at BX = 0
      if (not (data->egIEt.at(idx) >= 34)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleEG_14367823063771822943
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG20: ET >= 40 at BX = 0
      if (not (data->egIEt.at(idx) >= 40)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG18: ET >= 36 at BX = 0
      if (not (data->egIEt.at(idx) >= 36)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleEG_14367831859864844127
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(idx) >= 44)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleEG_14367831859864844383
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(idx) >= 44)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleEG_14367831859864844767
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(idx) >= 44)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG15: ET >= 30 at BX = 0
      if (not (data->egIEt.at(idx) >= 30)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleEG_14367836257911355231
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG23: ET >= 46 at BX = 0
      if (not (data->egIEt.at(idx) >= 46)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleEG_14367840655957867231
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(idx) >= 48)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG17: ET >= 34 at BX = 0
      if (not (data->egIEt.at(idx) >= 34)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleEG_14367845054004377695
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG25: ET >= 50 at BX = 0
      if (not (data->egIEt.at(idx) >= 50)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleEG_8902241742241126126
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG6: ET >= 12 at BX = 0
      if (not (data->egIEt.at(idx) >= 12)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG6: ET >= 12 at BX = 0
      if (not (data->egIEt.at(idx) >= 12)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_15894421920862285922
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx) >= 200)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx) >= 200)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_15903572090988376162
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET112: ET >= 224 at BX = 0
      if (not (data->jetIEt.at(idx) >= 224)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET112: ET >= 224 at BX = 0
      if (not (data->jetIEt.at(idx) >= 224)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_15912440717418279010
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_5010010172296896555
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET8: ET >= 16 at BX = 0
      if (not (data->jetIEt.at(idx) >= 16)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET8: ET >= 16 at BX = 0
      if (not (data->jetIEt.at(idx) >= 16)) continue;

                        // 3.0015 <= eta <= 5.0
              etaWindow1 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8281320341886584868
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET12: ET >= 24 at BX = 0
      if (not (data->jetIEt.at(idx) >= 24)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET12: ET >= 24 at BX = 0
      if (not (data->jetIEt.at(idx) >= 24)) continue;

                        // 3.0015 <= eta <= 5.0
              etaWindow1 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8281320350476519461
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET16: ET >= 32 at BX = 0
      if (not (data->jetIEt.at(idx) >= 32)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET16: ET >= 32 at BX = 0
      if (not (data->jetIEt.at(idx) >= 32)) continue;

                        // 3.0015 <= eta <= 5.0
              etaWindow1 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659083538331519699
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx) >= 60)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx) >= 60)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659084672202885843
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET32: ET >= 64 at BX = 0
      if (not (data->jetIEt.at(idx) >= 64)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET32: ET >= 64 at BX = 0
      if (not (data->jetIEt.at(idx) >= 64)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659086373009935059
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659156106098952915
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659160641584417491
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET48: ET >= 96 at BX = 0
      if (not (data->jetIEt.at(idx) >= 96)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET48: ET >= 96 at BX = 0
      if (not (data->jetIEt.at(idx) >= 96)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659228673866386131
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET50: ET >= 100 at BX = 0
      if (not (data->jetIEt.at(idx) >= 100)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET50: ET >= 100 at BX = 0
      if (not (data->jetIEt.at(idx) >= 100)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659301241633819347
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659301379072772819
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET68: ET >= 136 at BX = 0
      if (not (data->jetIEt.at(idx) >= 136)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET48: ET >= 96 at BX = 0
      if (not (data->jetIEt.at(idx) >= 96)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659370613945584339
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET72: ET >= 144 at BX = 0
      if (not (data->jetIEt.at(idx) >= 144)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET56: ET >= 112 at BX = 0
      if (not (data->jetIEt.at(idx) >= 112)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659374977632357075
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET76: ET >= 152 at BX = 0
      if (not (data->jetIEt.at(idx) >= 152)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET64: ET >= 128 at BX = 0
      if (not (data->jetIEt.at(idx) >= 128)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659439917537872595
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET84: ET >= 168 at BX = 0
      if (not (data->jetIEt.at(idx) >= 168)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET48: ET >= 96 at BX = 0
      if (not (data->jetIEt.at(idx) >= 96)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659444281224645331
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET88: ET >= 176 at BX = 0
      if (not (data->jetIEt.at(idx) >= 176)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET56: ET >= 112 at BX = 0
      if (not (data->jetIEt.at(idx) >= 112)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659444315584383699
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET84: ET >= 168 at BX = 0
      if (not (data->jetIEt.at(idx) >= 168)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET68: ET >= 136 at BX = 0
      if (not (data->jetIEt.at(idx) >= 136)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659446377168685779
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(idx) >= 160)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(idx) >= 160)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659448610551679699
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET88: ET >= 176 at BX = 0
      if (not (data->jetIEt.at(idx) >= 176)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET72: ET >= 144 at BX = 0
      if (not (data->jetIEt.at(idx) >= 144)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659513516097456851
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET92: ET >= 184 at BX = 0
      if (not (data->jetIEt.at(idx) >= 184)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET64: ET >= 128 at BX = 0
      if (not (data->jetIEt.at(idx) >= 128)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleJET_8659515749480450771
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET92: ET >= 184 at BX = 0
      if (not (data->jetIEt.at(idx) >= 184)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET76: ET >= 152 at BX = 0
      if (not (data->jetIEt.at(idx) >= 152)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleMU_14585777620730815295
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleMU_14585778097856672575
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleMU_14585796862184301375
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU7: ET >= 15 at BX = 0
      if (not (data->muonIEt.at(idx) >= 15)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU7: ET >= 15 at BX = 0
      if (not (data->muonIEt.at(idx) >= 15)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleMU_16961145543694661636
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU10: ET >= 21 at BX = 0
      if (not (data->muonIEt.at(idx) >= 21)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleMU_16961154507842811908
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU11: ET >= 23 at BX = 0
      if (not (data->muonIEt.at(idx) >= 23)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU4: ET >= 9 at BX = 0
      if (not (data->muonIEt.at(idx) >= 9)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleMU_16961157256621881348
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(idx) >= 25)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleMU_16961158905889323012
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU15: ET >= 31 at BX = 0
      if (not (data->muonIEt.at(idx) >= 31)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleMU_16961160005400950788
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU13: ET >= 27 at BX = 0
      if (not (data->muonIEt.at(idx) >= 27)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU6: ET >= 13 at BX = 0
      if (not (data->muonIEt.at(idx) >= 13)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleMU_16961163853691648004
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(idx) >= 25)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU8: ET >= 17 at BX = 0
      if (not (data->muonIEt.at(idx) >= 17)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleMU_18206240164090448142
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU10: ET >= 21 at BX = 0
      if (not (data->muonIEt.at(idx) >= 21)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3p5: ET >= 8 at BX = 0
      if (not (data->muonIEt.at(idx) >= 8)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleTAU_10196652277112847102
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU36: ET >= 72 at BX = 0
      if (not (data->tauIEt.at(idx) >= 72)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU36: ET >= 72 at BX = 0
      if (not (data->tauIEt.at(idx) >= 72)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleTAU_14808338227894500078
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU28: ET >= 56 at BX = 0
      if (not (data->tauIEt.at(idx) >= 56)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU28: ET >= 56 at BX = 0
      if (not (data->tauIEt.at(idx) >= 56)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleTAU_14808338292319009533
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU30: ET >= 60 at BX = 0
      if (not (data->tauIEt.at(idx) >= 60)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU30: ET >= 60 at BX = 0
      if (not (data->tauIEt.at(idx) >= 60)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleTAU_15233202657361500387
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU50: ET >= 100 at BX = 0
      if (not (data->tauIEt.at(idx) >= 100)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU50: ET >= 100 at BX = 0
      if (not (data->tauIEt.at(idx) >= 100)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleTAU_3279123247861152510
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU33: ET >= 66 at BX = 0
      if (not (data->tauIEt.at(idx) >= 66)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU33: ET >= 66 at BX = 0
      if (not (data->tauIEt.at(idx) >= 66)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleTAU_5584966257611717374
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU34: ET >= 68 at BX = 0
      if (not (data->tauIEt.at(idx) >= 68)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU34: ET >= 68 at BX = 0
      if (not (data->tauIEt.at(idx) >= 68)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleTAU_7890809267362282238
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU35: ET >= 70 at BX = 0
      if (not (data->tauIEt.at(idx) >= 70)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU35: ET >= 70 at BX = 0
      if (not (data->tauIEt.at(idx) >= 70)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
DoubleTAU_973280238110587646
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU32: ET >= 64 at BX = 0
      if (not (data->tauIEt.at(idx) >= 64)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU32: ET >= 64 at BX = 0
      if (not (data->tauIEt.at(idx) >= 64)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        


bool
MuonMuonCorrelation_16040223250608453060
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -1.6040625 <= eta <= 1.6040625
              etaWindow1 = ((-147 <= data->muonIEta.at(idx0)) and (data->muonIEta.at(idx0) <= 147));
            
          if (not etaWindow1) continue;
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -1.6040625 <= eta <= 1.6040625
              etaWindow1 = ((-147 <= data->muonIEta.at(idx1)) and (data->muonIEta.at(idx1) <= 147));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
          int iEta = data->muonIEta.at(idx0);
    unsigned int deltaIEta = abs(iEta - data->muonIEta.at(idx1));
  unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];

    minimum = (long long)(0.0 * POW10[3]);
    maximum = (long long)(1.8 * POW10[3]);
    if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        


bool
MuonMuonCorrelation_6226381454046753505
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU10: ET >= 21 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 21)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;
          int iPhi = data->muonIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhi.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];

    minimum = (long long)(1.0 * POW10[3]);
    maximum = (long long)(3.15 * POW10[3]);
    if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        


bool
MuonMuonCorrelation_7972376774213455602
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -1.4083125 <= eta <= 1.4083125
              etaWindow1 = ((-129 <= data->muonIEta.at(idx0)) and (data->muonIEta.at(idx0) <= 129));
            
          if (not etaWindow1) continue;
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -1.4083125 <= eta <= 1.4083125
              etaWindow1 = ((-129 <= data->muonIEta.at(idx1)) and (data->muonIEta.at(idx1) <= 129));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
        const std::string OS = "os";
    const std::string SS = "ss";
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
            int iEta = data->muonIEta.at(idx0);
    unsigned int deltaIEta = abs(iEta - data->muonIEta.at(idx1));
  unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];

    minimum = (long long)(0.0 * POW10[3]);
    maximum = (long long)(1.8 * POW10[3]);
    if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        


bool
MuonMuonCorrelation_8008405571232419574
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -1.6040625 <= eta <= 1.6040625
              etaWindow1 = ((-147 <= data->muonIEta.at(idx0)) and (data->muonIEta.at(idx0) <= 147));
            
          if (not etaWindow1) continue;
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -1.6040625 <= eta <= 1.6040625
              etaWindow1 = ((-147 <= data->muonIEta.at(idx1)) and (data->muonIEta.at(idx1) <= 147));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
        const std::string OS = "os";
    const std::string SS = "ss";
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
            int iEta = data->muonIEta.at(idx0);
    unsigned int deltaIEta = abs(iEta - data->muonIEta.at(idx1));
  unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];

    minimum = (long long)(0.0 * POW10[3]);
    maximum = (long long)(1.8 * POW10[3]);
    if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        


bool
MuonMuonCorrelation_8772456668275224612
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU10: ET >= 21 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 21)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;
          int iEta = data->muonIEta.at(idx0);
    unsigned int deltaIEta = abs(iEta - data->muonIEta.at(idx1));
  unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];

    minimum = (long long)(0.0 * POW10[3]);
    maximum = (long long)(1.8 * POW10[3]);
    if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
QuadJET_2680186536839014580
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET36: ET >= 72 at BX = 0
      if (not (data->jetIEt.at(idx) >= 72)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET36: ET >= 72 at BX = 0
      if (not (data->jetIEt.at(idx) >= 72)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET36: ET >= 72 at BX = 0
      if (not (data->jetIEt.at(idx) >= 72)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET36: ET >= 72 at BX = 0
      if (not (data->jetIEt.at(idx) >= 72)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
QuadJET_2751081844007168180
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
QuadJET_2825463805626214580
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET50: ET >= 100 at BX = 0
      if (not (data->jetIEt.at(idx) >= 100)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET50: ET >= 100 at BX = 0
      if (not (data->jetIEt.at(idx) >= 100)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET50: ET >= 100 at BX = 0
      if (not (data->jetIEt.at(idx) >= 100)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET50: ET >= 100 at BX = 0
      if (not (data->jetIEt.at(idx) >= 100)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
QuadJET_2899845767245260980
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
QuadMU_509409160461874775
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(3)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_1139634
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG2: ET >= 4 at BX = 0
      if (not (data->egIEt.at(idx) >= 4)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_1139637
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG5: ET >= 10 at BX = 0
      if (not (data->egIEt.at(idx) >= 10)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_1139639
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG7: ET >= 14 at BX = 0
      if (not (data->egIEt.at(idx) >= 14)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_12507428088042853440
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG6: ET >= 12 at BX = 0
      if (not (data->egIEt.at(idx) >= 12)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_12507579852048143424
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_12507579852056532032
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG18: ET >= 36 at BX = 0
      if (not (data->egIEt.at(idx) >= 36)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_12507579852182361152
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG20: ET >= 40 at BX = 0
      if (not (data->egIEt.at(idx) >= 40)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_12507579852184458304
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(idx) >= 44)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_12507579852186555456
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(idx) >= 48)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_12507579852188652608
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(idx) >= 52)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_12507579852190749760
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_12507579852316578880
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG30: ET >= 60 at BX = 0
      if (not (data->egIEt.at(idx) >= 60)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_12507579852318676032
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG32: ET >= 64 at BX = 0
      if (not (data->egIEt.at(idx) >= 64)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_12507579852320773184
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG34: ET >= 68 at BX = 0
      if (not (data->egIEt.at(idx) >= 68)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_12507579852322870336
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG36: ET >= 72 at BX = 0
      if (not (data->egIEt.at(idx) >= 72)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_12507579852324967488
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG38: ET >= 76 at BX = 0
      if (not (data->egIEt.at(idx) >= 76)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_14262501742662192051
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG25: ET >= 50 at BX = 0
      if (not (data->egIEt.at(idx) >= 50)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_14262501742930627507
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG27: ET >= 54 at BX = 0
      if (not (data->egIEt.at(idx) >= 54)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873072
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873076
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG14: ET >= 28 at BX = 0
      if (not (data->egIEt.at(idx) >= 28)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873077
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG15: ET >= 30 at BX = 0
      if (not (data->egIEt.at(idx) >= 30)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873079
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG17: ET >= 34 at BX = 0
      if (not (data->egIEt.at(idx) >= 34)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873080
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG18: ET >= 36 at BX = 0
      if (not (data->egIEt.at(idx) >= 36)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873200
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG20: ET >= 40 at BX = 0
      if (not (data->egIEt.at(idx) >= 40)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873203
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG23: ET >= 46 at BX = 0
      if (not (data->egIEt.at(idx) >= 46)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873204
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(idx) >= 48)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873206
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(idx) >= 52)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873208
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873328
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG30: ET >= 60 at BX = 0
      if (not (data->egIEt.at(idx) >= 60)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873330
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG32: ET >= 64 at BX = 0
      if (not (data->egIEt.at(idx) >= 64)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873332
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG34: ET >= 68 at BX = 0
      if (not (data->egIEt.at(idx) >= 68)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873334
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG36: ET >= 72 at BX = 0
      if (not (data->egIEt.at(idx) >= 72)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873336
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG38: ET >= 76 at BX = 0
      if (not (data->egIEt.at(idx) >= 76)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873456
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG40: ET >= 80 at BX = 0
      if (not (data->egIEt.at(idx) >= 80)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_145873461
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG45: ET >= 90 at BX = 0
      if (not (data->egIEt.at(idx) >= 90)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_6872811427746276593
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG18: ET >= 36 at BX = 0
      if (not (data->egIEt.at(idx) >= 36)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_6872943369141609713
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG20: ET >= 40 at BX = 0
      if (not (data->egIEt.at(idx) >= 40)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_6872945568164865265
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(idx) >= 44)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_6872947767188120817
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(idx) >= 48)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_6872949966211376369
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(idx) >= 52)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_6872952165234631921
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_6873084106629965041
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG30: ET >= 60 at BX = 0
      if (not (data->egIEt.at(idx) >= 60)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_6873086305653220593
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG32: ET >= 64 at BX = 0
      if (not (data->egIEt.at(idx) >= 64)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleEG_6873088504676476145
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG34: ET >= 68 at BX = 0
      if (not (data->egIEt.at(idx) >= 68)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
  bool
SingleETM_18699475376
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM30: ET >= 30.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 30.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleETM_18699475504
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM40: ET >= 40.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 40.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleETM_18699475632
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM50: ET >= 50.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 50.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleETM_18699475637
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM55: ET >= 55.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 55.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleETM_18699475760
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM60: ET >= 60.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 60.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleETM_18699475888
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM70: ET >= 70.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 70.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleETM_18699475893
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM75: ET >= 75.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 75.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleETM_18699476016
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM80: ET >= 80.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 80.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleETM_18699476021
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM85: ET >= 85.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 85.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleETM_18699476144
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM90: ET >= 90.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 90.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleETM_18699476149
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM95: ET >= 95.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 95.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleETM_2393532815408
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM100: ET >= 100.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 100.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleETM_2393532815664
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM120: ET >= 120.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 120.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleETT_18699589941
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETT25: ET >= 25.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 25.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleETT_18699590192
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETT40: ET >= 40.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 40.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleETT_18699590320
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETT50: ET >= 50.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 50.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleETT_18699590448
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETT60: ET >= 60.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 60.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleETT_18699590576
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETT70: ET >= 70.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 70.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleEXT_10104897634845317422
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_10371607390599051624
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_14414193171404190569
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_14715923867298343304
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_15455824636181887404
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_15455824636181887405
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_15455824636181887660
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_15455824636181887661
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_17118203077108929635
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_17561531836164454591
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_17561531836164454592
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_17561531836164454847
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_17561531836164454848
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_17833638493488257651
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_2629888000553438421
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_6395198100430131034
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_6873400283626490434
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_6912739140295604792
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_6926915327998939228
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_7332905005558692114
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
bool
SingleEXT_7332905005558692115
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

        
  bool
SingleHTM_19504782000
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTM50: ET >= 50.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 50.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleHTM_19504782128
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTM60: ET >= 60.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 60.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleHTM_19504782256
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTM70: ET >= 70.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 70.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleHTM_19504782384
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTM80: ET >= 80.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 80.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleHTM_2496612030512
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTM100: ET >= 100.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 100.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleHTM_2496612030768
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTM120: ET >= 120.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 120.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleHTM_2496612030896
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTM130: ET >= 130.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 130.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleHTM_2496612031024
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTM140: ET >= 140.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 140.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
  bool
SingleHTM_2496612031152
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTM150: ET >= 150.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 150.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleHTT_2496626710832
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT120: ET >= 120.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 120.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleHTT_2496626710837
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT125: ET >= 125.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 125.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleHTT_2496626711216
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT150: ET >= 150.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 150.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleHTT_2496626711344
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT160: ET >= 160.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 160.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleHTT_2496626726960
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT200: ET >= 200.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 200.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleHTT_2496626727216
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT220: ET >= 220.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 220.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleHTT_2496626727472
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT240: ET >= 240.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 240.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleHTT_2496626727605
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT255: ET >= 255.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 255.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleHTT_2496626727728
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT260: ET >= 260.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 260.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleHTT_2496626727856
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT270: ET >= 270.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 270.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleHTT_2496626727984
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT280: ET >= 280.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 280.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleHTT_2496626743344
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT300: ET >= 300.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 300.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleHTT_2496626743600
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT320: ET >= 320.0 at BX = 0
      if (not (data->sumEt.at(ii) >= 320.0)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleJET_15014918520304220377
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_156330552
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET8: ET >= 16 at BX = 0
      if (not (data->jetIEt.at(idx) >= 16)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_20010309810
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET12: ET >= 24 at BX = 0
      if (not (data->jetIEt.at(idx) >= 24)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_20010309814
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET16: ET >= 32 at BX = 0
      if (not (data->jetIEt.at(idx) >= 32)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_20010309936
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET20: ET >= 40 at BX = 0
      if (not (data->jetIEt.at(idx) >= 40)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_20010310069
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_20010310448
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_20010310832
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET90: ET >= 180 at BX = 0
      if (not (data->jetIEt.at(idx) >= 180)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_2561319655728
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_2561319655984
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET140: ET >= 280 at BX = 0
      if (not (data->jetIEt.at(idx) >= 280)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_2561319656112
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET150: ET >= 300 at BX = 0
      if (not (data->jetIEt.at(idx) >= 300)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_2561319656240
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET160: ET >= 320 at BX = 0
      if (not (data->jetIEt.at(idx) >= 320)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_2561319656368
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET170: ET >= 340 at BX = 0
      if (not (data->jetIEt.at(idx) >= 340)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_2561319656496
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET180: ET >= 360 at BX = 0
      if (not (data->jetIEt.at(idx) >= 360)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_2561319671856
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET200: ET >= 400 at BX = 0
      if (not (data->jetIEt.at(idx) >= 400)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5967545293332986055
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET16: ET >= 32 at BX = 0
      if (not (data->jetIEt.at(idx) >= 32)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5967545309707548871
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET20: ET >= 40 at BX = 0
      if (not (data->jetIEt.at(idx) >= 40)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5967545327155853511
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET32: ET >= 64 at BX = 0
      if (not (data->jetIEt.at(idx) >= 64)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5967545327558506695
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5967545344067287239
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5967545345141029063
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET48: ET >= 96 at BX = 0
      if (not (data->jetIEt.at(idx) >= 96)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5967545361247156423
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET50: ET >= 100 at BX = 0
      if (not (data->jetIEt.at(idx) >= 100)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5967545378427025607
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5974006375341702652
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET48: ET >= 96 at BX = 0
      if (not (data->jetIEt.at(idx) >= 96)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5974075644574252540
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET56: ET >= 112 at BX = 0
      if (not (data->jetIEt.at(idx) >= 112)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5974144913806802428
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET64: ET >= 128 at BX = 0
      if (not (data->jetIEt.at(idx) >= 128)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5974147112830057980
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET68: ET >= 136 at BX = 0
      if (not (data->jetIEt.at(idx) >= 136)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5974214183039352316
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET72: ET >= 144 at BX = 0
      if (not (data->jetIEt.at(idx) >= 144)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5974216382062607868
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET76: ET >= 152 at BX = 0
      if (not (data->jetIEt.at(idx) >= 152)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5974285651295157756
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET84: ET >= 168 at BX = 0
      if (not (data->jetIEt.at(idx) >= 168)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5974287850318413308
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET88: ET >= 176 at BX = 0
      if (not (data->jetIEt.at(idx) >= 176)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleJET_5974354920527707644
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET92: ET >= 184 at BX = 0
      if (not (data->jetIEt.at(idx) >= 184)) continue;

                        // -5.0 <= eta <= -3.0015
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
                        // 3.0015 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMBT0HFM_43640316738250417
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMinBiasHFM0)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                    

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleMBT0HFP_43640316738250801
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMinBiasHFP0)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                    

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleMBT1HFM_43640317006685873
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMinBiasHFM1)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                    

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleMBT1HFP_43640317006686257
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMinBiasHFP1)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                    

    pass = true;
    break;
  }

  return pass;
}

        
bool
SingleMU_11649248473972557216
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU30: ET >= 61 at BX = 0
      if (not (data->muonIEt.at(idx) >= 61)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEta.at(idx)) and (data->muonIEta.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_1272496
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_14769293018627052229
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_14769293071236239813
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_14769293105595978181
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_14769293122775847365
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU6: ET >= 13 at BX = 0
      if (not (data->muonIEt.at(idx) >= 13)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_14769293139955716549
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU7: ET >= 15 at BX = 0
      if (not (data->muonIEt.at(idx) >= 15)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_14769293157135585733
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU8: ET >= 17 at BX = 0
      if (not (data->muonIEt.at(idx) >= 17)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_16260934492399787300
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEta.at(idx)) and (data->muonIEta.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_17545683021081726533
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU10: ET >= 21 at BX = 0
      if (not (data->muonIEt.at(idx) >= 21)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_17545683059493081541
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(idx) >= 25)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_17545683093852819909
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU14: ET >= 29 at BX = 0
      if (not (data->muonIEt.at(idx) >= 29)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_17545683128212558277
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU16: ET >= 33 at BX = 0
      if (not (data->muonIEt.at(idx) >= 33)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_17545683162572296645
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU18: ET >= 37 at BX = 0
      if (not (data->muonIEt.at(idx) >= 37)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_17545685224156598725
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU20: ET >= 41 at BX = 0
      if (not (data->muonIEt.at(idx) >= 41)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_17545685258516337093
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU22: ET >= 45 at BX = 0
      if (not (data->muonIEt.at(idx) >= 45)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_17545685275696206277
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU23: ET >= 47 at BX = 0
      if (not (data->muonIEt.at(idx) >= 47)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_17545685292876075461
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU24: ET >= 49 at BX = 0
      if (not (data->muonIEt.at(idx) >= 49)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_17545685310055944645
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU25: ET >= 51 at BX = 0
      if (not (data->muonIEt.at(idx) >= 51)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_17545685361595552197
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU28: ET >= 57 at BX = 0
      if (not (data->muonIEt.at(idx) >= 57)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_17545687423179854277
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU30: ET >= 61 at BX = 0
      if (not (data->muonIEt.at(idx) >= 61)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_17545687457539592645
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU32: ET >= 65 at BX = 0
      if (not (data->muonIEt.at(idx) >= 65)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_17545687560618807749
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU38: ET >= 77 at BX = 0
      if (not (data->muonIEt.at(idx) >= 77)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_7037562455545169312
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU10: ET >= 21 at BX = 0
      if (not (data->muonIEt.at(idx) >= 21)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEta.at(idx)) and (data->muonIEta.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_7109620049583097248
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU14: ET >= 29 at BX = 0
      if (not (data->muonIEt.at(idx) >= 29)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEta.at(idx)) and (data->muonIEta.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_7145648846602061216
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU16: ET >= 33 at BX = 0
      if (not (data->muonIEt.at(idx) >= 33)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEta.at(idx)) and (data->muonIEta.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_7181677643621025184
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU18: ET >= 37 at BX = 0
      if (not (data->muonIEt.at(idx) >= 37)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEta.at(idx)) and (data->muonIEta.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_9343405464758863264
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU20: ET >= 41 at BX = 0
      if (not (data->muonIEt.at(idx) >= 41)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEta.at(idx)) and (data->muonIEta.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_9379434261777827232
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU22: ET >= 45 at BX = 0
      if (not (data->muonIEt.at(idx) >= 45)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEta.at(idx)) and (data->muonIEta.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleMU_9433477457306273184
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU25: ET >= 51 at BX = 0
      if (not (data->muonIEt.at(idx) >= 51)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEta.at(idx)) and (data->muonIEta.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleTAU_16608830939767017288
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU20: ET >= 40 at BX = 0
      if (not (data->tauIEt.at(idx) >= 40)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleTAU_16608831008486494024
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU24: ET >= 48 at BX = 0
      if (not (data->tauIEt.at(idx) >= 48)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleTAU_16608844133906550600
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU80: ET >= 160 at BX = 0
      if (not (data->tauIEt.at(idx) >= 160)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleTAU_218368042610145022
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU26: ET >= 52 at BX = 0
      if (not (data->tauIEt.at(idx) >= 52)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleTAU_22686292658
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU52: ET >= 104 at BX = 0
      if (not (data->tauIEt.at(idx) >= 104)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleTAU_3484211327656040900
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU100: ET >= 200 at BX = 0
      if (not (data->tauIEt.at(idx) >= 200)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
SingleTAU_3484215725702552004
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU120: ET >= 240 at BX = 0
      if (not (data->tauIEt.at(idx) >= 240)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
TripleEG_4430569450691365292
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG14: ET >= 28 at BX = 0
      if (not (data->egIEt.at(idx) >= 28)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // EG8: ET >= 16 at BX = 0
      if (not (data->egIEt.at(idx) >= 16)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
TripleEG_4430569691209534124
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG18: ET >= 36 at BX = 0
      if (not (data->egIEt.at(idx) >= 36)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG17: ET >= 34 at BX = 0
      if (not (data->egIEt.at(idx) >= 34)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // EG8: ET >= 16 at BX = 0
      if (not (data->egIEt.at(idx) >= 16)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
TripleJET_7918843258706703477
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET32: ET >= 64 at BX = 0
      if (not (data->jetIEt.at(idx) >= 64)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET32: ET >= 64 at BX = 0
      if (not (data->jetIEt.at(idx) >= 64)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET32: ET >= 64 at BX = 0
      if (not (data->jetIEt.at(idx) >= 64)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
TripleJET_7918897736071885941
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
TripleJET_7919042235951592565
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET50: ET >= 100 at BX = 0
      if (not (data->jetIEt.at(idx) >= 100)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
TripleJET_7921131308044366965
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
TripleJET_7921276581018186869
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET48: ET >= 96 at BX = 0
      if (not (data->jetIEt.at(idx) >= 96)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET48: ET >= 96 at BX = 0
      if (not (data->jetIEt.at(idx) >= 96)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET48: ET >= 96 at BX = 0
      if (not (data->jetIEt.at(idx) >= 96)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
TripleJET_7923455675625485429
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET50: ET >= 100 at BX = 0
      if (not (data->jetIEt.at(idx) >= 100)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET50: ET >= 100 at BX = 0
      if (not (data->jetIEt.at(idx) >= 100)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET50: ET >= 100 at BX = 0
      if (not (data->jetIEt.at(idx) >= 100)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
TripleJET_7930354149017105525
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET84: ET >= 168 at BX = 0
      if (not (data->jetIEt.at(idx) >= 168)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET68: ET >= 136 at BX = 0
      if (not (data->jetIEt.at(idx) >= 136)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET48: ET >= 96 at BX = 0
      if (not (data->jetIEt.at(idx) >= 96)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
TripleJET_7930493752634094709
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET88: ET >= 176 at BX = 0
      if (not (data->jetIEt.at(idx) >= 176)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET72: ET >= 144 at BX = 0
      if (not (data->jetIEt.at(idx) >= 144)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET56: ET >= 112 at BX = 0
      if (not (data->jetIEt.at(idx) >= 112)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
TripleJET_7932644363018286197
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET92: ET >= 184 at BX = 0
      if (not (data->jetIEt.at(idx) >= 184)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET76: ET >= 152 at BX = 0
      if (not (data->jetIEt.at(idx) >= 152)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET64: ET >= 128 at BX = 0
      if (not (data->jetIEt.at(idx) >= 128)) continue;

                        // -3.0015 <= eta <= 3.0015
              etaWindow1 = ((-69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 68));
            
          if (not etaWindow1) continue;
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
TripleMU_3324682852515662879
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
TripleMU_3324683539710430239
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

        
bool
TripleMU_3324692885559266335
(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

  

// generate algorithms
bool
L1_AlwaysTrue(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEXT_10371607390599051624(data) or (not SingleEXT_10371607390599051624(data));
}
bool
L1_BRIL_TRIG0_AND(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEXT_17833638493488257651(data);
}
bool
L1_BRIL_TRIG0_FstBunchInTrain(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEXT_6873400283626490434(data);
}
bool
L1_BRIL_TRIG0_OR(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEXT_6912739140295604792(data);
}
bool
L1_BRIL_TRIG0_delayedAND(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEXT_14715923867298343304(data);
}
bool
L1_BeamGasB1(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEXT_7332905005558692114(data);
}
bool
L1_BeamGasB2(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEXT_7332905005558692115(data);
}
bool
L1_BeamGasMinus(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEXT_14414193171404190569(data);
}
bool
L1_BeamGasPlus(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEXT_17118203077108929635(data);
}
bool
L1_BptxMinus(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEXT_6926915327998939228(data);
}
bool
L1_BptxOR(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEXT_10104897634845317422(data);
}
bool
L1_BptxPlus(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEXT_6395198100430131034(data);
}
bool
L1_BptxXOR(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return (SingleEXT_6395198100430131034(data) and (not SingleEXT_6926915327998939228(data))) or (SingleEXT_6926915327998939228(data) and (not SingleEXT_6395198100430131034(data)));
}
bool
L1_DoubleEG6_HTT255(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleEG_8902241742241126126(data) and SingleHTT_2496626727605(data);
}
bool
L1_DoubleEG_15_10(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleEG_14367282104050956127(data);
}
bool
L1_DoubleEG_18_17(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleEG_14367295298190490335(data);
}
bool
L1_DoubleEG_20_18(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleEG_14367823063771822943(data);
}
bool
L1_DoubleEG_22_10(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleEG_14367831859864844127(data);
}
bool
L1_DoubleEG_22_12(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleEG_14367831859864844383(data);
}
bool
L1_DoubleEG_22_15(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleEG_14367831859864844767(data);
}
bool
L1_DoubleEG_23_10(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleEG_14367836257911355231(data);
}
bool
L1_DoubleEG_24_17(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleEG_14367840655957867231(data);
}
bool
L1_DoubleEG_25_12(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleEG_14367845054004377695(data);
}
bool
L1_DoubleIsoTau28er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleTAU_14808338227894500078(data);
}
bool
L1_DoubleIsoTau30er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleTAU_14808338292319009533(data);
}
bool
L1_DoubleIsoTau32er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleTAU_973280238110587646(data);
}
bool
L1_DoubleIsoTau33er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleTAU_3279123247861152510(data);
}
bool
L1_DoubleIsoTau34er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleTAU_5584966257611717374(data);
}
bool
L1_DoubleIsoTau35er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleTAU_7890809267362282238(data);
}
bool
L1_DoubleIsoTau36er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleTAU_10196652277112847102(data);
}
bool
L1_DoubleJet12_ForwardBackward(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleJET_8281320341886584868(data);
}
bool
L1_DoubleJet16_ForwardBackward(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleJET_8281320350476519461(data);
}
bool
L1_DoubleJet8_ForwardBackward(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleJET_5010010172296896555(data);
}
bool
L1_DoubleJetC100(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleJET_15894421920862285922(data);
}
bool
L1_DoubleJetC112(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleJET_15903572090988376162(data);
}
bool
L1_DoubleJetC120(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleJET_15912440717418279010(data);
}
bool
L1_DoubleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleJET_8659156106098952915(data);
}
bool
L1_DoubleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleJET_8659228673866386131(data);
}
bool
L1_DoubleJetC60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleJET_8659301241633819347(data);
}
bool
L1_DoubleJetC60_ETM60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleJET_8659301241633819347(data) and SingleETM_18699475760(data);
}
bool
L1_DoubleJetC80(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleJET_8659446377168685779(data);
}
bool
L1_DoubleMu0(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleMU_14585777620730815295(data);
}
bool
L1_DoubleMu0_ETM40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleMU_14585777620730815295(data) and SingleETM_18699475504(data);
}
bool
L1_DoubleMu0_ETM55(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleMU_14585777620730815295(data) and SingleETM_18699475637(data);
}
bool
L1_DoubleMu0er1p4_dEta_Max1p8_OS(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return MuonMuonCorrelation_7972376774213455602(data);
}
bool
L1_DoubleMu0er1p6_dEta_Max1p8(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return MuonMuonCorrelation_16040223250608453060(data);
}
bool
L1_DoubleMu0er1p6_dEta_Max1p8_OS(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return MuonMuonCorrelation_8008405571232419574(data);
}
bool
L1_DoubleMu7_EG14(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleMU_14585796862184301375(data) and SingleEG_145873076(data);
}
bool
L1_DoubleMu7_EG7(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleMU_14585796862184301375(data) and SingleEG_1139639(data);
}
bool
L1_DoubleMuOpen(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleMU_14585778097856672575(data);
}
bool
L1_DoubleMu_10_0_dEta_Max1p8(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return MuonMuonCorrelation_8772456668275224612(data);
}
bool
L1_DoubleMu_10_3p5(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleMU_18206240164090448142(data);
}
bool
L1_DoubleMu_10_Open(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleMU_16961145543694661636(data);
}
bool
L1_DoubleMu_11_4(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleMU_16961154507842811908(data);
}
bool
L1_DoubleMu_12_5(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleMU_16961157256621881348(data);
}
bool
L1_DoubleMu_12_8(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleMU_16961163853691648004(data);
}
bool
L1_DoubleMu_13_6(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleMU_16961160005400950788(data);
}
bool
L1_DoubleMu_15_5(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleMU_16961158905889323012(data);
}
bool
L1_DoubleTau50er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return DoubleTAU_15233202657361500387(data);
}
bool
L1_EG25er_HTT125(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_14262501742662192051(data) and SingleHTT_2496626710837(data);
}
bool
L1_EG27er_HTT200(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_14262501742930627507(data) and SingleHTT_2496626726960(data);
}
bool
L1_ETM100(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETM_2393532815408(data);
}
bool
L1_ETM120(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETM_2393532815664(data);
}
bool
L1_ETM30(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETM_18699475376(data);
}
bool
L1_ETM40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETM_18699475504(data);
}
bool
L1_ETM50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETM_18699475632(data);
}
bool
L1_ETM60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETM_18699475760(data);
}
bool
L1_ETM70(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETM_18699475888(data);
}
bool
L1_ETM75(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETM_18699475893(data);
}
bool
L1_ETM75_Jet60_dPhi_Min0p4(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return CaloEsumCorrelation_16768129600233686289(data);
}
bool
L1_ETM80(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETM_18699476016(data);
}
bool
L1_ETM85(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETM_18699476021(data);
}
bool
L1_ETM90(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETM_18699476144(data);
}
bool
L1_ETM95(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETM_18699476149(data);
}
bool
L1_ETT25(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETT_18699589941(data);
}
bool
L1_ETT40_BptxAND(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETT_18699590192(data) and SingleEXT_10371607390599051624(data);
}
bool
L1_ETT50_BptxAND(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETT_18699590320(data) and SingleEXT_10371607390599051624(data);
}
bool
L1_ETT60_BptxAND(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETT_18699590448(data) and SingleEXT_10371607390599051624(data);
}
bool
L1_ETT70_BptxAND(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleETT_18699590576(data) and SingleEXT_10371607390599051624(data);
}
bool
L1_FirstBunchAfterTrain(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEXT_15455824636181887661(data) and SingleEXT_15455824636181887660(data) and (not SingleEXT_10104897634845317422(data)) and (not SingleEXT_17561531836164454591(data)) and (not SingleEXT_17561531836164454592(data));
}
bool
L1_FirstBunchInTrain(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return (not SingleEXT_17561531836164454848(data)) and (not SingleEXT_17561531836164454847(data)) and SingleEXT_10371607390599051624(data) and SingleEXT_15455824636181887404(data) and SingleEXT_15455824636181887405(data);
}
bool
L1_HTM100(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTM_2496612030512(data);
}
bool
L1_HTM120(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTM_2496612030768(data);
}
bool
L1_HTM130(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTM_2496612030896(data);
}
bool
L1_HTM140(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTM_2496612031024(data);
}
bool
L1_HTM150(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTM_2496612031152(data);
}
bool
L1_HTM50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTM_19504782000(data);
}
bool
L1_HTM60_HTT260(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTM_19504782128(data) and SingleHTT_2496626727728(data);
}
bool
L1_HTM70(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTM_19504782256(data);
}
bool
L1_HTM80(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTM_19504782384(data);
}
bool
L1_HTM80_HTT220(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTM_19504782384(data) and SingleHTT_2496626727216(data);
}
bool
L1_HTT120(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTT_2496626710832(data);
}
bool
L1_HTT160(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTT_2496626711344(data);
}
bool
L1_HTT200(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTT_2496626726960(data);
}
bool
L1_HTT220(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTT_2496626727216(data);
}
bool
L1_HTT240(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTT_2496626727472(data);
}
bool
L1_HTT255(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTT_2496626727605(data);
}
bool
L1_HTT270(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTT_2496626727856(data);
}
bool
L1_HTT280(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTT_2496626727984(data);
}
bool
L1_HTT300(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTT_2496626743344(data);
}
bool
L1_HTT320(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleHTT_2496626743600(data);
}
bool
L1_IsoEG18er_IsoTau24er_dEta_Min0p2(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return CaloCaloCorrelation_14433217633607694784(data);
}
bool
L1_IsoEG20er_IsoTau25er_dEta_Min0p2(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return CaloCaloCorrelation_14500771630165735872(data);
}
bool
L1_IsoEG22er_IsoTau26er_dEta_Min0p2(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return CaloCaloCorrelation_14501897532220062144(data);
}
bool
L1_IsoEG22er_Tau20er_dEta_Min0p2(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return CaloCaloCorrelation_2085349615869404806(data);
}
bool
L1_IsoEG24_DoubleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and DoubleJET_8659084672202885843(data);
}
bool
L1_IsoEG24_DoubleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and DoubleJET_8659086373009935059(data);
}
bool
L1_IsoEG24_DoubleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and DoubleJET_8659156106098952915(data);
}
bool
L1_IsoEG24_DoubleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and DoubleJET_8659160641584417491(data);
}
bool
L1_IsoEG24_DoubleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and DoubleJET_8659228673866386131(data);
}
bool
L1_IsoEG24_JetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and SingleJET_5967545327155853511(data);
}
bool
L1_IsoEG24_JetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and SingleJET_5967545327558506695(data);
}
bool
L1_IsoEG24_JetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and SingleJET_5967545344067287239(data);
}
bool
L1_IsoEG24_JetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and SingleJET_5967545345141029063(data);
}
bool
L1_IsoEG24_JetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and SingleJET_5967545361247156423(data);
}
bool
L1_IsoEG24_TripleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and TripleJET_7918843258706703477(data);
}
bool
L1_IsoEG24_TripleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and TripleJET_7918897736071885941(data);
}
bool
L1_IsoEG24_TripleJetC35_50_60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and TripleJET_7919042235951592565(data);
}
bool
L1_IsoEG24_TripleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and TripleJET_7921131308044366965(data);
}
bool
L1_IsoEG24_TripleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and TripleJET_7921276581018186869(data);
}
bool
L1_IsoEG24_TripleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data) and TripleJET_7923455675625485429(data);
}
bool
L1_IsoEG26_DoubleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and DoubleJET_8659084672202885843(data);
}
bool
L1_IsoEG26_DoubleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and DoubleJET_8659086373009935059(data);
}
bool
L1_IsoEG26_DoubleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and DoubleJET_8659156106098952915(data);
}
bool
L1_IsoEG26_DoubleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and DoubleJET_8659160641584417491(data);
}
bool
L1_IsoEG26_DoubleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and DoubleJET_8659228673866386131(data);
}
bool
L1_IsoEG26_JetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and SingleJET_5967545327155853511(data);
}
bool
L1_IsoEG26_JetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and SingleJET_5967545327558506695(data);
}
bool
L1_IsoEG26_JetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and SingleJET_5967545344067287239(data);
}
bool
L1_IsoEG26_JetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and SingleJET_5967545345141029063(data);
}
bool
L1_IsoEG26_JetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and SingleJET_5967545361247156423(data);
}
bool
L1_IsoEG26_TripleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and TripleJET_7918843258706703477(data);
}
bool
L1_IsoEG26_TripleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and TripleJET_7918897736071885941(data);
}
bool
L1_IsoEG26_TripleJetC35_50_60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and TripleJET_7919042235951592565(data);
}
bool
L1_IsoEG26_TripleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and TripleJET_7921131308044366965(data);
}
bool
L1_IsoEG26_TripleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and TripleJET_7921276581018186869(data);
}
bool
L1_IsoEG26_TripleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data) and TripleJET_7923455675625485429(data);
}
bool
L1_IsoEG28_DoubleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and DoubleJET_8659084672202885843(data);
}
bool
L1_IsoEG28_DoubleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and DoubleJET_8659086373009935059(data);
}
bool
L1_IsoEG28_DoubleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and DoubleJET_8659156106098952915(data);
}
bool
L1_IsoEG28_DoubleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and DoubleJET_8659160641584417491(data);
}
bool
L1_IsoEG28_DoubleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and DoubleJET_8659228673866386131(data);
}
bool
L1_IsoEG28_JetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and SingleJET_5967545327155853511(data);
}
bool
L1_IsoEG28_JetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and SingleJET_5967545327558506695(data);
}
bool
L1_IsoEG28_JetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and SingleJET_5967545344067287239(data);
}
bool
L1_IsoEG28_JetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and SingleJET_5967545345141029063(data);
}
bool
L1_IsoEG28_JetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and SingleJET_5967545361247156423(data);
}
bool
L1_IsoEG28_TripleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and TripleJET_7918843258706703477(data);
}
bool
L1_IsoEG28_TripleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and TripleJET_7918897736071885941(data);
}
bool
L1_IsoEG28_TripleJetC35_50_60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and TripleJET_7919042235951592565(data);
}
bool
L1_IsoEG28_TripleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and TripleJET_7921131308044366965(data);
}
bool
L1_IsoEG28_TripleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and TripleJET_7921276581018186869(data);
}
bool
L1_IsoEG28_TripleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data) and TripleJET_7923455675625485429(data);
}
bool
L1_IsoEG32_DoubleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and DoubleJET_8659084672202885843(data);
}
bool
L1_IsoEG32_DoubleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and DoubleJET_8659086373009935059(data);
}
bool
L1_IsoEG32_DoubleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and DoubleJET_8659156106098952915(data);
}
bool
L1_IsoEG32_DoubleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and DoubleJET_8659160641584417491(data);
}
bool
L1_IsoEG32_DoubleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and DoubleJET_8659228673866386131(data);
}
bool
L1_IsoEG32_JetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and SingleJET_5967545327155853511(data);
}
bool
L1_IsoEG32_JetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and SingleJET_5967545327558506695(data);
}
bool
L1_IsoEG32_JetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and SingleJET_5967545344067287239(data);
}
bool
L1_IsoEG32_JetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and SingleJET_5967545345141029063(data);
}
bool
L1_IsoEG32_JetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and SingleJET_5967545361247156423(data);
}
bool
L1_IsoEG32_TripleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and TripleJET_7918843258706703477(data);
}
bool
L1_IsoEG32_TripleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and TripleJET_7918897736071885941(data);
}
bool
L1_IsoEG32_TripleJetC35_50_60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and TripleJET_7919042235951592565(data);
}
bool
L1_IsoEG32_TripleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and TripleJET_7921131308044366965(data);
}
bool
L1_IsoEG32_TripleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and TripleJET_7921276581018186869(data);
}
bool
L1_IsoEG32_TripleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data) and TripleJET_7923455675625485429(data);
}
bool
L1_IsoEG34_DoubleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and DoubleJET_8659084672202885843(data);
}
bool
L1_IsoEG34_DoubleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and DoubleJET_8659086373009935059(data);
}
bool
L1_IsoEG34_DoubleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and DoubleJET_8659156106098952915(data);
}
bool
L1_IsoEG34_DoubleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and DoubleJET_8659160641584417491(data);
}
bool
L1_IsoEG34_DoubleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and DoubleJET_8659228673866386131(data);
}
bool
L1_IsoEG34_JetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and SingleJET_5967545327155853511(data);
}
bool
L1_IsoEG34_JetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and SingleJET_5967545327558506695(data);
}
bool
L1_IsoEG34_JetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and SingleJET_5967545344067287239(data);
}
bool
L1_IsoEG34_JetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and SingleJET_5967545345141029063(data);
}
bool
L1_IsoEG34_JetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and SingleJET_5967545361247156423(data);
}
bool
L1_IsoEG34_TripleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and TripleJET_7918843258706703477(data);
}
bool
L1_IsoEG34_TripleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and TripleJET_7918897736071885941(data);
}
bool
L1_IsoEG34_TripleJetC35_50_60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and TripleJET_7919042235951592565(data);
}
bool
L1_IsoEG34_TripleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and TripleJET_7921131308044366965(data);
}
bool
L1_IsoEG34_TripleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and TripleJET_7921276581018186869(data);
}
bool
L1_IsoEG34_TripleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data) and TripleJET_7923455675625485429(data);
}
bool
L1_IsoEG38_DoubleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and DoubleJET_8659084672202885843(data);
}
bool
L1_IsoEG38_DoubleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and DoubleJET_8659086373009935059(data);
}
bool
L1_IsoEG38_DoubleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and DoubleJET_8659156106098952915(data);
}
bool
L1_IsoEG38_DoubleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and DoubleJET_8659160641584417491(data);
}
bool
L1_IsoEG38_DoubleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and DoubleJET_8659228673866386131(data);
}
bool
L1_IsoEG38_JetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and SingleJET_5967545327155853511(data);
}
bool
L1_IsoEG38_JetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and SingleJET_5967545327558506695(data);
}
bool
L1_IsoEG38_JetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and SingleJET_5967545344067287239(data);
}
bool
L1_IsoEG38_JetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and SingleJET_5967545345141029063(data);
}
bool
L1_IsoEG38_JetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and SingleJET_5967545361247156423(data);
}
bool
L1_IsoEG38_TripleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and TripleJET_7918843258706703477(data);
}
bool
L1_IsoEG38_TripleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and TripleJET_7918897736071885941(data);
}
bool
L1_IsoEG38_TripleJetC35_50_60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and TripleJET_7919042235951592565(data);
}
bool
L1_IsoEG38_TripleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and TripleJET_7921131308044366965(data);
}
bool
L1_IsoEG38_TripleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and TripleJET_7921276581018186869(data);
}
bool
L1_IsoEG38_TripleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852324967488(data) and TripleJET_7923455675625485429(data);
}
bool
L1_IsolatedBunch(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return (not SingleEXT_17561531836164454848(data)) and (not SingleEXT_17561531836164454847(data)) and SingleEXT_10371607390599051624(data) and (not SingleEXT_17561531836164454591(data)) and (not SingleEXT_17561531836164454592(data));
}
bool
L1_Jet32_DoubleMu_10_0_dPhi_Jet_Mu0_Max0p4_dPhi_Mu_Mu_Min1p0(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return CaloMuonCorrelation_10791898730651162912(data) and MuonMuonCorrelation_6226381454046753505(data);
}
bool
L1_Jet32_Mu0_EG10_dPhi_Jet_Mu_Max0p4_dPhi_Mu_EG_Min1p0(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return CaloMuonCorrelation_10791898730651162912(data) and CaloMuonCorrelation_10674670645420326056(data);
}
bool
L1_MU20_EG15(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685224156598725(data) and SingleEG_145873077(data);
}
bool
L1_MinimumBiasHF0_AND(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMBT0HFP_43640316738250801(data) and SingleMBT0HFM_43640316738250417(data);
}
bool
L1_MinimumBiasHF0_AND_BptxAND(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return (SingleMBT0HFP_43640316738250801(data) and SingleMBT0HFM_43640316738250417(data)) and SingleEXT_10371607390599051624(data);
}
bool
L1_MinimumBiasHF0_OR(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMBT0HFP_43640316738250801(data) or SingleMBT0HFM_43640316738250417(data);
}
bool
L1_MinimumBiasHF0_OR_BptxAND(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return (SingleMBT0HFP_43640316738250801(data) or SingleMBT0HFM_43640316738250417(data)) and SingleEXT_10371607390599051624(data);
}
bool
L1_MinimumBiasHF1_AND(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMBT1HFP_43640317006686257(data) and SingleMBT1HFM_43640317006685873(data);
}
bool
L1_MinimumBiasHF1_AND_BptxAND(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return (SingleMBT1HFP_43640317006686257(data) and SingleMBT1HFM_43640317006685873(data)) and SingleEXT_10371607390599051624(data);
}
bool
L1_MinimumBiasHF1_OR(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMBT1HFP_43640317006686257(data) or SingleMBT1HFM_43640317006685873(data);
}
bool
L1_MinimumBiasHF1_OR_BptxAND(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return (SingleMBT1HFP_43640317006686257(data) or SingleMBT1HFM_43640317006685873(data)) and SingleEXT_10371607390599051624(data);
}
bool
L1_Mu0er_ETM40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_16260934492399787300(data) and SingleETM_18699475504(data);
}
bool
L1_Mu0er_ETM55(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_16260934492399787300(data) and SingleETM_18699475637(data);
}
bool
L1_Mu10er_ETM30(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_7037562455545169312(data) and SingleETM_18699475376(data);
}
bool
L1_Mu10er_ETM50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_7037562455545169312(data) and SingleETM_18699475632(data);
}
bool
L1_Mu12_EG10(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545683059493081541(data) and SingleEG_145873072(data);
}
bool
L1_Mu14er_ETM30(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_7109620049583097248(data) and SingleETM_18699475376(data);
}
bool
L1_Mu16er_Tau20er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_7145648846602061216(data) and SingleTAU_16608830939767017288(data);
}
bool
L1_Mu16er_Tau24er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_7145648846602061216(data) and SingleTAU_16608831008486494024(data);
}
bool
L1_Mu18er_IsoTau26er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_7181677643621025184(data) and SingleTAU_218368042610145022(data);
}
bool
L1_Mu18er_Tau20er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_7181677643621025184(data) and SingleTAU_16608830939767017288(data);
}
bool
L1_Mu18er_Tau24er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_7181677643621025184(data) and SingleTAU_16608831008486494024(data);
}
bool
L1_Mu20_EG10(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685224156598725(data) and SingleEG_145873072(data);
}
bool
L1_Mu20_EG17(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685224156598725(data) and SingleEG_145873079(data);
}
bool
L1_Mu20_IsoEG6(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685224156598725(data) and SingleEG_12507428088042853440(data);
}
bool
L1_Mu20er_IsoTau26er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_9343405464758863264(data) and SingleTAU_218368042610145022(data);
}
bool
L1_Mu22er_IsoTau26er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_9379434261777827232(data) and SingleTAU_218368042610145022(data);
}
bool
L1_Mu23_EG10(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685275696206277(data) and SingleEG_145873072(data);
}
bool
L1_Mu23_IsoEG10(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685275696206277(data) and SingleEG_12507579852048143424(data);
}
bool
L1_Mu25er_IsoTau26er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_9433477457306273184(data) and SingleTAU_218368042610145022(data);
}
bool
L1_Mu3_JetC120(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293071236239813(data) and SingleJET_15014918520304220377(data);
}
bool
L1_Mu3_JetC120_dEta_Max0p4_dPhi_Max0p4(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return CaloMuonCorrelation_15993852978349077723(data);
}
bool
L1_Mu3_JetC16(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293071236239813(data) and SingleJET_5967545293332986055(data);
}
bool
L1_Mu3_JetC16_dEta_Max0p4_dPhi_Max0p4(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return CaloMuonCorrelation_16240387826857744385(data);
}
bool
L1_Mu3_JetC60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293071236239813(data) and SingleJET_5967545378427025607(data);
}
bool
L1_Mu3_JetC60_dEta_Max0p4_dPhi_Max0p4(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return CaloMuonCorrelation_16240389188362377217(data);
}
bool
L1_Mu5_EG15(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293105595978181(data) and SingleEG_145873077(data);
}
bool
L1_Mu5_EG20(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293105595978181(data) and SingleEG_145873200(data);
}
bool
L1_Mu5_EG23(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293105595978181(data) and SingleEG_145873203(data);
}
bool
L1_Mu5_IsoEG18(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293105595978181(data) and SingleEG_12507579852056532032(data);
}
bool
L1_Mu5_IsoEG20(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293105595978181(data) and SingleEG_12507579852182361152(data);
}
bool
L1_Mu6_DoubleEG10(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293122775847365(data) and DoubleEG_14367260113818400607(data);
}
bool
L1_Mu6_DoubleEG17(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293122775847365(data) and DoubleEG_14367290900143979231(data);
}
bool
L1_Mu6_HTT200(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293122775847365(data) and SingleHTT_2496626726960(data);
}
bool
L1_Mu8_HTT150(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293157135585733(data) and SingleHTT_2496626711216(data);
}
bool
L1_NotBptxOR(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return not SingleEXT_10104897634845317422(data);
}
bool
L1_QuadJetC36_Tau52(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return QuadJET_2680186536839014580(data) and SingleTAU_22686292658(data);
}
bool
L1_QuadJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return QuadJET_2751081844007168180(data);
}
bool
L1_QuadJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return QuadJET_2825463805626214580(data);
}
bool
L1_QuadJetC60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return QuadJET_2899845767245260980(data);
}
bool
L1_QuadMu0(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return QuadMU_509409160461874775(data);
}
bool
L1_SingleEG10(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_145873072(data);
}
bool
L1_SingleEG15(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_145873077(data);
}
bool
L1_SingleEG18(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_145873080(data);
}
bool
L1_SingleEG24(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_145873204(data);
}
bool
L1_SingleEG26(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_145873206(data);
}
bool
L1_SingleEG28(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_145873208(data);
}
bool
L1_SingleEG2_BptxAND(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_1139634(data) and SingleEXT_10371607390599051624(data);
}
bool
L1_SingleEG30(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_145873328(data);
}
bool
L1_SingleEG32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_145873330(data);
}
bool
L1_SingleEG34(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_145873332(data);
}
bool
L1_SingleEG36(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_145873334(data);
}
bool
L1_SingleEG38(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_145873336(data);
}
bool
L1_SingleEG40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_145873456(data);
}
bool
L1_SingleEG45(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_145873461(data);
}
bool
L1_SingleEG5(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_1139637(data);
}
bool
L1_SingleIsoEG18(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852056532032(data);
}
bool
L1_SingleIsoEG18er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_6872811427746276593(data);
}
bool
L1_SingleIsoEG20(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852182361152(data);
}
bool
L1_SingleIsoEG20er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_6872943369141609713(data);
}
bool
L1_SingleIsoEG22(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852184458304(data);
}
bool
L1_SingleIsoEG22er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_6872945568164865265(data);
}
bool
L1_SingleIsoEG24(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852186555456(data);
}
bool
L1_SingleIsoEG24er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_6872947767188120817(data);
}
bool
L1_SingleIsoEG26(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852188652608(data);
}
bool
L1_SingleIsoEG26er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_6872949966211376369(data);
}
bool
L1_SingleIsoEG28(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852190749760(data);
}
bool
L1_SingleIsoEG28er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_6872952165234631921(data);
}
bool
L1_SingleIsoEG30(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852316578880(data);
}
bool
L1_SingleIsoEG30er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_6873084106629965041(data);
}
bool
L1_SingleIsoEG32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852318676032(data);
}
bool
L1_SingleIsoEG32er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_6873086305653220593(data);
}
bool
L1_SingleIsoEG34(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852320773184(data);
}
bool
L1_SingleIsoEG34er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_6873088504676476145(data);
}
bool
L1_SingleIsoEG36(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEG_12507579852322870336(data);
}
bool
L1_SingleJet120(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_2561319655728(data);
}
bool
L1_SingleJet12_BptxAND(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_20010309810(data) and SingleEXT_10371607390599051624(data);
}
bool
L1_SingleJet140(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_2561319655984(data);
}
bool
L1_SingleJet150(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_2561319656112(data);
}
bool
L1_SingleJet16(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_20010309814(data);
}
bool
L1_SingleJet160(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_2561319656240(data);
}
bool
L1_SingleJet170(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_2561319656368(data);
}
bool
L1_SingleJet180(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_2561319656496(data);
}
bool
L1_SingleJet20(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_20010309936(data);
}
bool
L1_SingleJet200(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_2561319671856(data);
}
bool
L1_SingleJet35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_20010310069(data);
}
bool
L1_SingleJet60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_20010310448(data);
}
bool
L1_SingleJet8_BptxAND(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_156330552(data) and SingleEXT_10371607390599051624(data);
}
bool
L1_SingleJet90(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_20010310832(data);
}
bool
L1_SingleJetC20_NotBptxOR(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_5967545309707548871(data) and (not SingleEXT_10104897634845317422(data));
}
bool
L1_SingleJetC20_NotBptxOR_3BX(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_5967545309707548871(data) and (not SingleEXT_17561531836164454847(data)) and (not SingleEXT_10104897634845317422(data)) and (not SingleEXT_17561531836164454591(data));
}
bool
L1_SingleJetC32_NotBptxOR(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_5967545327155853511(data) and (not SingleEXT_10104897634845317422(data));
}
bool
L1_SingleJetC32_NotBptxOR_3BX(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_5967545327155853511(data) and (not SingleEXT_17561531836164454847(data)) and (not SingleEXT_10104897634845317422(data)) and (not SingleEXT_17561531836164454591(data));
}
bool
L1_SingleJetC40_NotBptxOR_3BX(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleJET_5967545344067287239(data) and (not SingleEXT_17561531836164454847(data)) and (not SingleEXT_10104897634845317422(data)) and (not SingleEXT_17561531836164454591(data));
}
bool
L1_SingleMu10_LowQ(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545683021081726533(data);
}
bool
L1_SingleMu12(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545683059493081541(data);
}
bool
L1_SingleMu14(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545683093852819909(data);
}
bool
L1_SingleMu14er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_7109620049583097248(data);
}
bool
L1_SingleMu16(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545683128212558277(data);
}
bool
L1_SingleMu16er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_7145648846602061216(data);
}
bool
L1_SingleMu18(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545683162572296645(data);
}
bool
L1_SingleMu18er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_7181677643621025184(data);
}
bool
L1_SingleMu20(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685224156598725(data);
}
bool
L1_SingleMu20er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_9343405464758863264(data);
}
bool
L1_SingleMu22(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685258516337093(data);
}
bool
L1_SingleMu22er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_9379434261777827232(data);
}
bool
L1_SingleMu24_DoubleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and DoubleJET_8659084672202885843(data);
}
bool
L1_SingleMu24_DoubleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and DoubleJET_8659086373009935059(data);
}
bool
L1_SingleMu24_DoubleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and DoubleJET_8659156106098952915(data);
}
bool
L1_SingleMu24_DoubleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and DoubleJET_8659160641584417491(data);
}
bool
L1_SingleMu24_DoubleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and DoubleJET_8659228673866386131(data);
}
bool
L1_SingleMu24_JetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and SingleJET_5967545327155853511(data);
}
bool
L1_SingleMu24_JetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and SingleJET_5967545327558506695(data);
}
bool
L1_SingleMu24_JetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and SingleJET_5967545344067287239(data);
}
bool
L1_SingleMu24_JetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and SingleJET_5967545345141029063(data);
}
bool
L1_SingleMu24_JetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and SingleJET_5967545361247156423(data);
}
bool
L1_SingleMu24_TripleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and TripleJET_7918843258706703477(data);
}
bool
L1_SingleMu24_TripleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and TripleJET_7918897736071885941(data);
}
bool
L1_SingleMu24_TripleJetC35_50_60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and TripleJET_7919042235951592565(data);
}
bool
L1_SingleMu24_TripleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and TripleJET_7921131308044366965(data);
}
bool
L1_SingleMu24_TripleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and TripleJET_7921276581018186869(data);
}
bool
L1_SingleMu24_TripleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685292876075461(data) and TripleJET_7923455675625485429(data);
}
bool
L1_SingleMu25(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685310055944645(data);
}
bool
L1_SingleMu25_DoubleJetC30_(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685310055944645(data) and DoubleJET_8659083538331519699(data);
}
bool
L1_SingleMu25er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_9433477457306273184(data);
}
bool
L1_SingleMu28_DoubleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and DoubleJET_8659084672202885843(data);
}
bool
L1_SingleMu28_DoubleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and DoubleJET_8659086373009935059(data);
}
bool
L1_SingleMu28_DoubleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and DoubleJET_8659156106098952915(data);
}
bool
L1_SingleMu28_DoubleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and DoubleJET_8659160641584417491(data);
}
bool
L1_SingleMu28_DoubleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and DoubleJET_8659228673866386131(data);
}
bool
L1_SingleMu28_JetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and SingleJET_5967545327155853511(data);
}
bool
L1_SingleMu28_JetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and SingleJET_5967545327558506695(data);
}
bool
L1_SingleMu28_JetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and SingleJET_5967545344067287239(data);
}
bool
L1_SingleMu28_JetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and SingleJET_5967545345141029063(data);
}
bool
L1_SingleMu28_JetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and SingleJET_5967545361247156423(data);
}
bool
L1_SingleMu28_TripleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and TripleJET_7918843258706703477(data);
}
bool
L1_SingleMu28_TripleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and TripleJET_7918897736071885941(data);
}
bool
L1_SingleMu28_TripleJetC35_50_60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and TripleJET_7919042235951592565(data);
}
bool
L1_SingleMu28_TripleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and TripleJET_7921131308044366965(data);
}
bool
L1_SingleMu28_TripleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and TripleJET_7921276581018186869(data);
}
bool
L1_SingleMu28_TripleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545685361595552197(data) and TripleJET_7923455675625485429(data);
}
bool
L1_SingleMu3(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293071236239813(data);
}
bool
L1_SingleMu30(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687423179854277(data);
}
bool
L1_SingleMu30er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_11649248473972557216(data);
}
bool
L1_SingleMu32_DoubleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and DoubleJET_8659084672202885843(data);
}
bool
L1_SingleMu32_DoubleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and DoubleJET_8659086373009935059(data);
}
bool
L1_SingleMu32_DoubleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and DoubleJET_8659156106098952915(data);
}
bool
L1_SingleMu32_DoubleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and DoubleJET_8659160641584417491(data);
}
bool
L1_SingleMu32_DoubleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and DoubleJET_8659228673866386131(data);
}
bool
L1_SingleMu32_JetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and SingleJET_5967545327155853511(data);
}
bool
L1_SingleMu32_JetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and SingleJET_5967545327558506695(data);
}
bool
L1_SingleMu32_JetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and SingleJET_5967545344067287239(data);
}
bool
L1_SingleMu32_JetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and SingleJET_5967545345141029063(data);
}
bool
L1_SingleMu32_JetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and SingleJET_5967545361247156423(data);
}
bool
L1_SingleMu32_TripleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and TripleJET_7918843258706703477(data);
}
bool
L1_SingleMu32_TripleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and TripleJET_7918897736071885941(data);
}
bool
L1_SingleMu32_TripleJetC35_50_60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and TripleJET_7919042235951592565(data);
}
bool
L1_SingleMu32_TripleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and TripleJET_7921131308044366965(data);
}
bool
L1_SingleMu32_TripleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and TripleJET_7921276581018186869(data);
}
bool
L1_SingleMu32_TripleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687457539592645(data) and TripleJET_7923455675625485429(data);
}
bool
L1_SingleMu38_DoubleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and DoubleJET_8659084672202885843(data);
}
bool
L1_SingleMu38_DoubleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and DoubleJET_8659086373009935059(data);
}
bool
L1_SingleMu38_DoubleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and DoubleJET_8659156106098952915(data);
}
bool
L1_SingleMu38_DoubleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and DoubleJET_8659160641584417491(data);
}
bool
L1_SingleMu38_DoubleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and DoubleJET_8659228673866386131(data);
}
bool
L1_SingleMu38_JetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and SingleJET_5967545327155853511(data);
}
bool
L1_SingleMu38_JetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and SingleJET_5967545327558506695(data);
}
bool
L1_SingleMu38_JetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and SingleJET_5967545344067287239(data);
}
bool
L1_SingleMu38_JetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and SingleJET_5967545345141029063(data);
}
bool
L1_SingleMu38_JetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and SingleJET_5967545361247156423(data);
}
bool
L1_SingleMu38_TripleJetC32(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and TripleJET_7918843258706703477(data);
}
bool
L1_SingleMu38_TripleJetC35(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and TripleJET_7918897736071885941(data);
}
bool
L1_SingleMu38_TripleJetC35_50_60(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and TripleJET_7919042235951592565(data);
}
bool
L1_SingleMu38_TripleJetC40(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and TripleJET_7921131308044366965(data);
}
bool
L1_SingleMu38_TripleJetC48(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and TripleJET_7921276581018186869(data);
}
bool
L1_SingleMu38_TripleJetC50(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_17545687560618807749(data) and TripleJET_7923455675625485429(data);
}
bool
L1_SingleMu5(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293105595978181(data);
}
bool
L1_SingleMu7(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293139955716549(data);
}
bool
L1_SingleMuCosmics(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_1272496(data);
}
bool
L1_SingleMuOpen(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293018627052229(data);
}
bool
L1_SingleMuOpen_NotBptxOR(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293018627052229(data) and (not SingleEXT_10104897634845317422(data));
}
bool
L1_SingleMuOpen_NotBptxOR_3BX(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleMU_14769293018627052229(data) and (not SingleEXT_17561531836164454847(data)) and (not SingleEXT_10104897634845317422(data)) and (not SingleEXT_17561531836164454591(data));
}
bool
L1_SingleTau100er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleTAU_3484211327656040900(data);
}
bool
L1_SingleTau120er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleTAU_3484215725702552004(data);
}
bool
L1_SingleTau80er(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleTAU_16608844133906550600(data);
}
bool
L1_TripleEG_14_10_8(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return TripleEG_4430569450691365292(data);
}
bool
L1_TripleEG_18_17_8(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return TripleEG_4430569691209534124(data);
}
bool
L1_TripleJet_84_68_48_VBF(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return TripleJET_7930354149017105525(data) or (DoubleJET_8659444315584383699(data) and SingleJET_5974006375341702652(data)) or (DoubleJET_8659439917537872595(data) and SingleJET_5974147112830057980(data)) or (DoubleJET_8659301379072772819(data) and SingleJET_5974285651295157756(data));
}
bool
L1_TripleJet_88_72_56_VBF(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return TripleJET_7930493752634094709(data) or (DoubleJET_8659448610551679699(data) and SingleJET_5974075644574252540(data)) or (DoubleJET_8659444281224645331(data) and SingleJET_5974214183039352316(data)) or (DoubleJET_8659370613945584339(data) and SingleJET_5974287850318413308(data));
}
bool
L1_TripleJet_92_76_64_VBF(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return TripleJET_7932644363018286197(data) or (DoubleJET_8659515749480450771(data) and SingleJET_5974144913806802428(data)) or (DoubleJET_8659513516097456851(data) and SingleJET_5974216382062607868(data)) or (DoubleJET_8659374977632357075(data) and SingleJET_5974354920527707644(data));
}
bool
L1_TripleMu0(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return TripleMU_3324682852515662879(data);
}
bool
L1_TripleMu_5_0_0(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return TripleMU_3324683539710430239(data);
}
bool
L1_TripleMu_5_5_3(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return TripleMU_3324692885559266335(data);
}
bool
L1_ZeroBias(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEXT_10371607390599051624(data);
}
bool
L1_ZeroBias_FirstCollidingBunch(TTreeReaderValue<L1Analysis::L1AnalysisL1UpgradeDataFormat>& data)
{
  return SingleEXT_2629888000553438421(data);
}


std::string getNameFromId(const int index)
{
  static const std::pair<int, std::string> id2name[] = {
          std::make_pair(206, "L1_AlwaysTrue"),          std::make_pair(240, "L1_BRIL_TRIG0_AND"),          std::make_pair(243, "L1_BRIL_TRIG0_FstBunchInTrain"),          std::make_pair(242, "L1_BRIL_TRIG0_OR"),          std::make_pair(241, "L1_BRIL_TRIG0_delayedAND"),          std::make_pair(223, "L1_BeamGasB1"),          std::make_pair(224, "L1_BeamGasB2"),          std::make_pair(222, "L1_BeamGasMinus"),          std::make_pair(221, "L1_BeamGasPlus"),          std::make_pair(208, "L1_BptxMinus"),          std::make_pair(209, "L1_BptxOR"),          std::make_pair(207, "L1_BptxPlus"),          std::make_pair(220, "L1_BptxXOR"),          std::make_pair(176, "L1_DoubleEG6_HTT255"),          std::make_pair(71, "L1_DoubleEG_15_10"),          std::make_pair(72, "L1_DoubleEG_18_17"),          std::make_pair(73, "L1_DoubleEG_20_18"),          std::make_pair(75, "L1_DoubleEG_22_10"),          std::make_pair(196, "L1_DoubleEG_22_12"),          std::make_pair(197, "L1_DoubleEG_22_15"),          std::make_pair(76, "L1_DoubleEG_23_10"),          std::make_pair(77, "L1_DoubleEG_24_17"),          std::make_pair(277, "L1_DoubleEG_25_12"),          std::make_pair(109, "L1_DoubleIsoTau28er"),          std::make_pair(110, "L1_DoubleIsoTau30er"),          std::make_pair(111, "L1_DoubleIsoTau32er"),          std::make_pair(264, "L1_DoubleIsoTau33er"),          std::make_pair(265, "L1_DoubleIsoTau34er"),          std::make_pair(266, "L1_DoubleIsoTau35er"),          std::make_pair(278, "L1_DoubleIsoTau36er"),          std::make_pair(213, "L1_DoubleJet12_ForwardBackward"),          std::make_pair(214, "L1_DoubleJet16_ForwardBackward"),          std::make_pair(212, "L1_DoubleJet8_ForwardBackward"),          std::make_pair(96, "L1_DoubleJetC100"),          std::make_pair(97, "L1_DoubleJetC112"),          std::make_pair(98, "L1_DoubleJetC120"),          std::make_pair(92, "L1_DoubleJetC40"),          std::make_pair(93, "L1_DoubleJetC50"),          std::make_pair(94, "L1_DoubleJetC60"),          std::make_pair(181, "L1_DoubleJetC60_ETM60"),          std::make_pair(95, "L1_DoubleJetC80"),          std::make_pair(24, "L1_DoubleMu0"),          std::make_pair(256, "L1_DoubleMu0_ETM40"),          std::make_pair(257, "L1_DoubleMu0_ETM55"),          std::make_pair(267, "L1_DoubleMu0er1p4_dEta_Max1p8_OS"),          std::make_pair(32, "L1_DoubleMu0er1p6_dEta_Max1p8"),          std::make_pair(33, "L1_DoubleMu0er1p6_dEta_Max1p8_OS"),          std::make_pair(166, "L1_DoubleMu7_EG14"),          std::make_pair(167, "L1_DoubleMu7_EG7"),          std::make_pair(23, "L1_DoubleMuOpen"),          std::make_pair(35, "L1_DoubleMu_10_0_dEta_Max1p8"),          std::make_pair(26, "L1_DoubleMu_10_3p5"),          std::make_pair(25, "L1_DoubleMu_10_Open"),          std::make_pair(27, "L1_DoubleMu_11_4"),          std::make_pair(28, "L1_DoubleMu_12_5"),          std::make_pair(31, "L1_DoubleMu_12_8"),          std::make_pair(29, "L1_DoubleMu_13_6"),          std::make_pair(30, "L1_DoubleMu_15_5"),          std::make_pair(114, "L1_DoubleTau50er"),          std::make_pair(175, "L1_EG25er_HTT125"),          std::make_pair(174, "L1_EG27er_HTT200"),          std::make_pair(142, "L1_ETM100"),          std::make_pair(143, "L1_ETM120"),          std::make_pair(137, "L1_ETM30"),          std::make_pair(138, "L1_ETM40"),          std::make_pair(139, "L1_ETM50"),          std::make_pair(140, "L1_ETM60"),          std::make_pair(141, "L1_ETM70"),          std::make_pair(272, "L1_ETM75"),          std::make_pair(271, "L1_ETM75_Jet60_dPhi_Min0p4"),          std::make_pair(125, "L1_ETM80"),          std::make_pair(273, "L1_ETM85"),          std::make_pair(126, "L1_ETM90"),          std::make_pair(274, "L1_ETM95"),          std::make_pair(135, "L1_ETT25"),          std::make_pair(136, "L1_ETT40_BptxAND"),          std::make_pair(244, "L1_ETT50_BptxAND"),          std::make_pair(245, "L1_ETT60_BptxAND"),          std::make_pair(193, "L1_ETT70_BptxAND"),          std::make_pair(282, "L1_FirstBunchAfterTrain"),          std::make_pair(281, "L1_FirstBunchInTrain"),          std::make_pair(130, "L1_HTM100"),          std::make_pair(131, "L1_HTM120"),          std::make_pair(132, "L1_HTM130"),          std::make_pair(133, "L1_HTM140"),          std::make_pair(134, "L1_HTM150"),          std::make_pair(127, "L1_HTM50"),          std::make_pair(188, "L1_HTM60_HTT260"),          std::make_pair(128, "L1_HTM70"),          std::make_pair(129, "L1_HTM80"),          std::make_pair(187, "L1_HTM80_HTT220"),          std::make_pair(115, "L1_HTT120"),          std::make_pair(116, "L1_HTT160"),          std::make_pair(117, "L1_HTT200"),          std::make_pair(118, "L1_HTT220"),          std::make_pair(119, "L1_HTT240"),          std::make_pair(120, "L1_HTT255"),          std::make_pair(121, "L1_HTT270"),          std::make_pair(122, "L1_HTT280"),          std::make_pair(123, "L1_HTT300"),          std::make_pair(124, "L1_HTT320"),          std::make_pair(269, "L1_IsoEG18er_IsoTau24er_dEta_Min0p2"),          std::make_pair(270, "L1_IsoEG20er_IsoTau25er_dEta_Min0p2"),          std::make_pair(268, "L1_IsoEG22er_IsoTau26er_dEta_Min0p2"),          std::make_pair(199, "L1_IsoEG22er_Tau20er_dEta_Min0p2"),          std::make_pair(330, "L1_IsoEG24_DoubleJetC32"),          std::make_pair(336, "L1_IsoEG24_DoubleJetC35"),          std::make_pair(342, "L1_IsoEG24_DoubleJetC40"),          std::make_pair(348, "L1_IsoEG24_DoubleJetC48"),          std::make_pair(354, "L1_IsoEG24_DoubleJetC50"),          std::make_pair(300, "L1_IsoEG24_JetC32"),          std::make_pair(306, "L1_IsoEG24_JetC35"),          std::make_pair(312, "L1_IsoEG24_JetC40"),          std::make_pair(318, "L1_IsoEG24_JetC48"),          std::make_pair(324, "L1_IsoEG24_JetC50"),          std::make_pair(360, "L1_IsoEG24_TripleJetC32"),          std::make_pair(366, "L1_IsoEG24_TripleJetC35"),          std::make_pair(390, "L1_IsoEG24_TripleJetC35_50_60"),          std::make_pair(372, "L1_IsoEG24_TripleJetC40"),          std::make_pair(378, "L1_IsoEG24_TripleJetC48"),          std::make_pair(384, "L1_IsoEG24_TripleJetC50"),          std::make_pair(331, "L1_IsoEG26_DoubleJetC32"),          std::make_pair(337, "L1_IsoEG26_DoubleJetC35"),          std::make_pair(343, "L1_IsoEG26_DoubleJetC40"),          std::make_pair(349, "L1_IsoEG26_DoubleJetC48"),          std::make_pair(355, "L1_IsoEG26_DoubleJetC50"),          std::make_pair(301, "L1_IsoEG26_JetC32"),          std::make_pair(307, "L1_IsoEG26_JetC35"),          std::make_pair(313, "L1_IsoEG26_JetC40"),          std::make_pair(319, "L1_IsoEG26_JetC48"),          std::make_pair(325, "L1_IsoEG26_JetC50"),          std::make_pair(361, "L1_IsoEG26_TripleJetC32"),          std::make_pair(367, "L1_IsoEG26_TripleJetC35"),          std::make_pair(391, "L1_IsoEG26_TripleJetC35_50_60"),          std::make_pair(373, "L1_IsoEG26_TripleJetC40"),          std::make_pair(379, "L1_IsoEG26_TripleJetC48"),          std::make_pair(385, "L1_IsoEG26_TripleJetC50"),          std::make_pair(332, "L1_IsoEG28_DoubleJetC32"),          std::make_pair(338, "L1_IsoEG28_DoubleJetC35"),          std::make_pair(344, "L1_IsoEG28_DoubleJetC40"),          std::make_pair(350, "L1_IsoEG28_DoubleJetC48"),          std::make_pair(356, "L1_IsoEG28_DoubleJetC50"),          std::make_pair(302, "L1_IsoEG28_JetC32"),          std::make_pair(308, "L1_IsoEG28_JetC35"),          std::make_pair(314, "L1_IsoEG28_JetC40"),          std::make_pair(320, "L1_IsoEG28_JetC48"),          std::make_pair(326, "L1_IsoEG28_JetC50"),          std::make_pair(362, "L1_IsoEG28_TripleJetC32"),          std::make_pair(368, "L1_IsoEG28_TripleJetC35"),          std::make_pair(392, "L1_IsoEG28_TripleJetC35_50_60"),          std::make_pair(374, "L1_IsoEG28_TripleJetC40"),          std::make_pair(380, "L1_IsoEG28_TripleJetC48"),          std::make_pair(386, "L1_IsoEG28_TripleJetC50"),          std::make_pair(333, "L1_IsoEG32_DoubleJetC32"),          std::make_pair(339, "L1_IsoEG32_DoubleJetC35"),          std::make_pair(345, "L1_IsoEG32_DoubleJetC40"),          std::make_pair(351, "L1_IsoEG32_DoubleJetC48"),          std::make_pair(357, "L1_IsoEG32_DoubleJetC50"),          std::make_pair(303, "L1_IsoEG32_JetC32"),          std::make_pair(309, "L1_IsoEG32_JetC35"),          std::make_pair(315, "L1_IsoEG32_JetC40"),          std::make_pair(321, "L1_IsoEG32_JetC48"),          std::make_pair(327, "L1_IsoEG32_JetC50"),          std::make_pair(363, "L1_IsoEG32_TripleJetC32"),          std::make_pair(369, "L1_IsoEG32_TripleJetC35"),          std::make_pair(393, "L1_IsoEG32_TripleJetC35_50_60"),          std::make_pair(375, "L1_IsoEG32_TripleJetC40"),          std::make_pair(381, "L1_IsoEG32_TripleJetC48"),          std::make_pair(387, "L1_IsoEG32_TripleJetC50"),          std::make_pair(334, "L1_IsoEG34_DoubleJetC32"),          std::make_pair(340, "L1_IsoEG34_DoubleJetC35"),          std::make_pair(346, "L1_IsoEG34_DoubleJetC40"),          std::make_pair(352, "L1_IsoEG34_DoubleJetC48"),          std::make_pair(358, "L1_IsoEG34_DoubleJetC50"),          std::make_pair(304, "L1_IsoEG34_JetC32"),          std::make_pair(310, "L1_IsoEG34_JetC35"),          std::make_pair(316, "L1_IsoEG34_JetC40"),          std::make_pair(322, "L1_IsoEG34_JetC48"),          std::make_pair(328, "L1_IsoEG34_JetC50"),          std::make_pair(364, "L1_IsoEG34_TripleJetC32"),          std::make_pair(370, "L1_IsoEG34_TripleJetC35"),          std::make_pair(394, "L1_IsoEG34_TripleJetC35_50_60"),          std::make_pair(376, "L1_IsoEG34_TripleJetC40"),          std::make_pair(382, "L1_IsoEG34_TripleJetC48"),          std::make_pair(388, "L1_IsoEG34_TripleJetC50"),          std::make_pair(335, "L1_IsoEG38_DoubleJetC32"),          std::make_pair(341, "L1_IsoEG38_DoubleJetC35"),          std::make_pair(347, "L1_IsoEG38_DoubleJetC40"),          std::make_pair(353, "L1_IsoEG38_DoubleJetC48"),          std::make_pair(359, "L1_IsoEG38_DoubleJetC50"),          std::make_pair(305, "L1_IsoEG38_JetC32"),          std::make_pair(311, "L1_IsoEG38_JetC35"),          std::make_pair(317, "L1_IsoEG38_JetC40"),          std::make_pair(323, "L1_IsoEG38_JetC48"),          std::make_pair(329, "L1_IsoEG38_JetC50"),          std::make_pair(365, "L1_IsoEG38_TripleJetC32"),          std::make_pair(371, "L1_IsoEG38_TripleJetC35"),          std::make_pair(395, "L1_IsoEG38_TripleJetC35_50_60"),          std::make_pair(377, "L1_IsoEG38_TripleJetC40"),          std::make_pair(383, "L1_IsoEG38_TripleJetC48"),          std::make_pair(389, "L1_IsoEG38_TripleJetC50"),          std::make_pair(219, "L1_IsolatedBunch"),          std::make_pair(179, "L1_Jet32_DoubleMu_10_0_dPhi_Jet_Mu0_Max0p4_dPhi_Mu_Mu_Min1p0"),          std::make_pair(180, "L1_Jet32_Mu0_EG10_dPhi_Jet_Mu_Max0p4_dPhi_Mu_EG_Min1p0"),          std::make_pair(198, "L1_MU20_EG15"),          std::make_pair(251, "L1_MinimumBiasHF0_AND"),          std::make_pair(247, "L1_MinimumBiasHF0_AND_BptxAND"),          std::make_pair(250, "L1_MinimumBiasHF0_OR"),          std::make_pair(246, "L1_MinimumBiasHF0_OR_BptxAND"),          std::make_pair(253, "L1_MinimumBiasHF1_AND"),          std::make_pair(249, "L1_MinimumBiasHF1_AND_BptxAND"),          std::make_pair(252, "L1_MinimumBiasHF1_OR"),          std::make_pair(248, "L1_MinimumBiasHF1_OR_BptxAND"),          std::make_pair(182, "L1_Mu0er_ETM40"),          std::make_pair(183, "L1_Mu0er_ETM55"),          std::make_pair(184, "L1_Mu10er_ETM30"),          std::make_pair(185, "L1_Mu10er_ETM50"),          std::make_pair(149, "L1_Mu12_EG10"),          std::make_pair(186, "L1_Mu14er_ETM30"),          std::make_pair(154, "L1_Mu16er_Tau20er"),          std::make_pair(155, "L1_Mu16er_Tau24er"),          std::make_pair(158, "L1_Mu18er_IsoTau26er"),          std::make_pair(156, "L1_Mu18er_Tau20er"),          std::make_pair(157, "L1_Mu18er_Tau24er"),          std::make_pair(150, "L1_Mu20_EG10"),          std::make_pair(151, "L1_Mu20_EG17"),          std::make_pair(275, "L1_Mu20_IsoEG6"),          std::make_pair(159, "L1_Mu20er_IsoTau26er"),          std::make_pair(279, "L1_Mu22er_IsoTau26er"),          std::make_pair(153, "L1_Mu23_EG10"),          std::make_pair(152, "L1_Mu23_IsoEG10"),          std::make_pair(280, "L1_Mu25er_IsoTau26er"),          std::make_pair(217, "L1_Mu3_JetC120"),          std::make_pair(210, "L1_Mu3_JetC120_dEta_Max0p4_dPhi_Max0p4"),          std::make_pair(215, "L1_Mu3_JetC16"),          std::make_pair(170, "L1_Mu3_JetC16_dEta_Max0p4_dPhi_Max0p4"),          std::make_pair(216, "L1_Mu3_JetC60"),          std::make_pair(171, "L1_Mu3_JetC60_dEta_Max0p4_dPhi_Max0p4"),          std::make_pair(144, "L1_Mu5_EG15"),          std::make_pair(145, "L1_Mu5_EG20"),          std::make_pair(146, "L1_Mu5_EG23"),          std::make_pair(147, "L1_Mu5_IsoEG18"),          std::make_pair(148, "L1_Mu5_IsoEG20"),          std::make_pair(169, "L1_Mu6_DoubleEG10"),          std::make_pair(168, "L1_Mu6_DoubleEG17"),          std::make_pair(172, "L1_Mu6_HTT200"),          std::make_pair(173, "L1_Mu8_HTT150"),          std::make_pair(254, "L1_NotBptxOR"),          std::make_pair(177, "L1_QuadJetC36_Tau52"),          std::make_pair(102, "L1_QuadJetC40"),          std::make_pair(103, "L1_QuadJetC50"),          std::make_pair(104, "L1_QuadJetC60"),          std::make_pair(38, "L1_QuadMu0"),          std::make_pair(40, "L1_SingleEG10"),          std::make_pair(41, "L1_SingleEG15"),          std::make_pair(42, "L1_SingleEG18"),          std::make_pair(43, "L1_SingleEG24"),          std::make_pair(44, "L1_SingleEG26"),          std::make_pair(45, "L1_SingleEG28"),          std::make_pair(192, "L1_SingleEG2_BptxAND"),          std::make_pair(46, "L1_SingleEG30"),          std::make_pair(258, "L1_SingleEG32"),          std::make_pair(48, "L1_SingleEG34"),          std::make_pair(259, "L1_SingleEG36"),          std::make_pair(260, "L1_SingleEG38"),          std::make_pair(50, "L1_SingleEG40"),          std::make_pair(52, "L1_SingleEG45"),          std::make_pair(39, "L1_SingleEG5"),          std::make_pair(53, "L1_SingleIsoEG18"),          std::make_pair(62, "L1_SingleIsoEG18er"),          std::make_pair(54, "L1_SingleIsoEG20"),          std::make_pair(63, "L1_SingleIsoEG20er"),          std::make_pair(55, "L1_SingleIsoEG22"),          std::make_pair(64, "L1_SingleIsoEG22er"),          std::make_pair(56, "L1_SingleIsoEG24"),          std::make_pair(65, "L1_SingleIsoEG24er"),          std::make_pair(57, "L1_SingleIsoEG26"),          std::make_pair(66, "L1_SingleIsoEG26er"),          std::make_pair(59, "L1_SingleIsoEG28"),          std::make_pair(68, "L1_SingleIsoEG28er"),          std::make_pair(60, "L1_SingleIsoEG30"),          std::make_pair(69, "L1_SingleIsoEG30er"),          std::make_pair(261, "L1_SingleIsoEG32"),          std::make_pair(263, "L1_SingleIsoEG32er"),          std::make_pair(61, "L1_SingleIsoEG34"),          std::make_pair(70, "L1_SingleIsoEG34er"),          std::make_pair(262, "L1_SingleIsoEG36"),          std::make_pair(85, "L1_SingleJet120"),          std::make_pair(195, "L1_SingleJet12_BptxAND"),          std::make_pair(86, "L1_SingleJet140"),          std::make_pair(87, "L1_SingleJet150"),          std::make_pair(80, "L1_SingleJet16"),          std::make_pair(88, "L1_SingleJet160"),          std::make_pair(89, "L1_SingleJet170"),          std::make_pair(90, "L1_SingleJet180"),          std::make_pair(81, "L1_SingleJet20"),          std::make_pair(91, "L1_SingleJet200"),          std::make_pair(82, "L1_SingleJet35"),          std::make_pair(83, "L1_SingleJet60"),          std::make_pair(194, "L1_SingleJet8_BptxAND"),          std::make_pair(84, "L1_SingleJet90"),          std::make_pair(191, "L1_SingleJetC20_NotBptxOR"),          std::make_pair(201, "L1_SingleJetC20_NotBptxOR_3BX"),          std::make_pair(190, "L1_SingleJetC32_NotBptxOR"),          std::make_pair(202, "L1_SingleJetC32_NotBptxOR_3BX"),          std::make_pair(203, "L1_SingleJetC40_NotBptxOR_3BX"),          std::make_pair(14, "L1_SingleMu10_LowQ"),          std::make_pair(6, "L1_SingleMu12"),          std::make_pair(7, "L1_SingleMu14"),          std::make_pair(16, "L1_SingleMu14er"),          std::make_pair(8, "L1_SingleMu16"),          std::make_pair(17, "L1_SingleMu16er"),          std::make_pair(9, "L1_SingleMu18"),          std::make_pair(18, "L1_SingleMu18er"),          std::make_pair(10, "L1_SingleMu20"),          std::make_pair(19, "L1_SingleMu20er"),          std::make_pair(11, "L1_SingleMu22"),          std::make_pair(20, "L1_SingleMu22er"),          std::make_pair(416, "L1_SingleMu24_DoubleJetC32"),          std::make_pair(420, "L1_SingleMu24_DoubleJetC35"),          std::make_pair(424, "L1_SingleMu24_DoubleJetC40"),          std::make_pair(428, "L1_SingleMu24_DoubleJetC48"),          std::make_pair(432, "L1_SingleMu24_DoubleJetC50"),          std::make_pair(396, "L1_SingleMu24_JetC32"),          std::make_pair(400, "L1_SingleMu24_JetC35"),          std::make_pair(404, "L1_SingleMu24_JetC40"),          std::make_pair(408, "L1_SingleMu24_JetC48"),          std::make_pair(412, "L1_SingleMu24_JetC50"),          std::make_pair(436, "L1_SingleMu24_TripleJetC32"),          std::make_pair(440, "L1_SingleMu24_TripleJetC35"),          std::make_pair(456, "L1_SingleMu24_TripleJetC35_50_60"),          std::make_pair(444, "L1_SingleMu24_TripleJetC40"),          std::make_pair(448, "L1_SingleMu24_TripleJetC48"),          std::make_pair(452, "L1_SingleMu24_TripleJetC50"),          std::make_pair(12, "L1_SingleMu25"),          std::make_pair(15, "L1_SingleMu25_DoubleJetC30_"),          std::make_pair(21, "L1_SingleMu25er"),          std::make_pair(417, "L1_SingleMu28_DoubleJetC32"),          std::make_pair(421, "L1_SingleMu28_DoubleJetC35"),          std::make_pair(425, "L1_SingleMu28_DoubleJetC40"),          std::make_pair(429, "L1_SingleMu28_DoubleJetC48"),          std::make_pair(433, "L1_SingleMu28_DoubleJetC50"),          std::make_pair(397, "L1_SingleMu28_JetC32"),          std::make_pair(401, "L1_SingleMu28_JetC35"),          std::make_pair(405, "L1_SingleMu28_JetC40"),          std::make_pair(409, "L1_SingleMu28_JetC48"),          std::make_pair(413, "L1_SingleMu28_JetC50"),          std::make_pair(437, "L1_SingleMu28_TripleJetC32"),          std::make_pair(441, "L1_SingleMu28_TripleJetC35"),          std::make_pair(457, "L1_SingleMu28_TripleJetC35_50_60"),          std::make_pair(445, "L1_SingleMu28_TripleJetC40"),          std::make_pair(449, "L1_SingleMu28_TripleJetC48"),          std::make_pair(453, "L1_SingleMu28_TripleJetC50"),          std::make_pair(3, "L1_SingleMu3"),          std::make_pair(13, "L1_SingleMu30"),          std::make_pair(22, "L1_SingleMu30er"),          std::make_pair(418, "L1_SingleMu32_DoubleJetC32"),          std::make_pair(422, "L1_SingleMu32_DoubleJetC35"),          std::make_pair(426, "L1_SingleMu32_DoubleJetC40"),          std::make_pair(430, "L1_SingleMu32_DoubleJetC48"),          std::make_pair(434, "L1_SingleMu32_DoubleJetC50"),          std::make_pair(398, "L1_SingleMu32_JetC32"),          std::make_pair(402, "L1_SingleMu32_JetC35"),          std::make_pair(406, "L1_SingleMu32_JetC40"),          std::make_pair(410, "L1_SingleMu32_JetC48"),          std::make_pair(414, "L1_SingleMu32_JetC50"),          std::make_pair(438, "L1_SingleMu32_TripleJetC32"),          std::make_pair(442, "L1_SingleMu32_TripleJetC35"),          std::make_pair(458, "L1_SingleMu32_TripleJetC35_50_60"),          std::make_pair(446, "L1_SingleMu32_TripleJetC40"),          std::make_pair(450, "L1_SingleMu32_TripleJetC48"),          std::make_pair(454, "L1_SingleMu32_TripleJetC50"),          std::make_pair(419, "L1_SingleMu38_DoubleJetC32"),          std::make_pair(423, "L1_SingleMu38_DoubleJetC35"),          std::make_pair(427, "L1_SingleMu38_DoubleJetC40"),          std::make_pair(431, "L1_SingleMu38_DoubleJetC48"),          std::make_pair(435, "L1_SingleMu38_DoubleJetC50"),          std::make_pair(399, "L1_SingleMu38_JetC32"),          std::make_pair(403, "L1_SingleMu38_JetC35"),          std::make_pair(407, "L1_SingleMu38_JetC40"),          std::make_pair(411, "L1_SingleMu38_JetC48"),          std::make_pair(415, "L1_SingleMu38_JetC50"),          std::make_pair(439, "L1_SingleMu38_TripleJetC32"),          std::make_pair(443, "L1_SingleMu38_TripleJetC35"),          std::make_pair(459, "L1_SingleMu38_TripleJetC35_50_60"),          std::make_pair(447, "L1_SingleMu38_TripleJetC40"),          std::make_pair(451, "L1_SingleMu38_TripleJetC48"),          std::make_pair(455, "L1_SingleMu38_TripleJetC50"),          std::make_pair(4, "L1_SingleMu5"),          std::make_pair(5, "L1_SingleMu7"),          std::make_pair(218, "L1_SingleMuCosmics"),          std::make_pair(2, "L1_SingleMuOpen"),          std::make_pair(189, "L1_SingleMuOpen_NotBptxOR"),          std::make_pair(200, "L1_SingleMuOpen_NotBptxOR_3BX"),          std::make_pair(106, "L1_SingleTau100er"),          std::make_pair(107, "L1_SingleTau120er"),          std::make_pair(105, "L1_SingleTau80er"),          std::make_pair(78, "L1_TripleEG_14_10_8"),          std::make_pair(79, "L1_TripleEG_18_17_8"),          std::make_pair(99, "L1_TripleJet_84_68_48_VBF"),          std::make_pair(100, "L1_TripleJet_88_72_56_VBF"),          std::make_pair(101, "L1_TripleJet_92_76_64_VBF"),          std::make_pair(36, "L1_TripleMu0"),          std::make_pair(276, "L1_TripleMu_5_0_0"),          std::make_pair(37, "L1_TripleMu_5_5_3"),          std::make_pair(0, "L1_ZeroBias"),          std::make_pair(211, "L1_ZeroBias_FirstCollidingBunch")      };

  static const std::map<int, std::string> Id2Name(id2name, id2name + sizeof(id2name) / sizeof(id2name[0]));
  const std::map<int, std::string>::const_iterator rc = Id2Name.find(index);
  std::string name;
  if (rc != Id2Name.end()) name = rc->second;
  return name;
}


int getIdFromName(const std::string& name)
{
  static const std::pair<std::string, int> name2id[] = {
          std::make_pair("L1_AlwaysTrue", 206),          std::make_pair("L1_BRIL_TRIG0_AND", 240),          std::make_pair("L1_BRIL_TRIG0_FstBunchInTrain", 243),          std::make_pair("L1_BRIL_TRIG0_OR", 242),          std::make_pair("L1_BRIL_TRIG0_delayedAND", 241),          std::make_pair("L1_BeamGasB1", 223),          std::make_pair("L1_BeamGasB2", 224),          std::make_pair("L1_BeamGasMinus", 222),          std::make_pair("L1_BeamGasPlus", 221),          std::make_pair("L1_BptxMinus", 208),          std::make_pair("L1_BptxOR", 209),          std::make_pair("L1_BptxPlus", 207),          std::make_pair("L1_BptxXOR", 220),          std::make_pair("L1_DoubleEG6_HTT255", 176),          std::make_pair("L1_DoubleEG_15_10", 71),          std::make_pair("L1_DoubleEG_18_17", 72),          std::make_pair("L1_DoubleEG_20_18", 73),          std::make_pair("L1_DoubleEG_22_10", 75),          std::make_pair("L1_DoubleEG_22_12", 196),          std::make_pair("L1_DoubleEG_22_15", 197),          std::make_pair("L1_DoubleEG_23_10", 76),          std::make_pair("L1_DoubleEG_24_17", 77),          std::make_pair("L1_DoubleEG_25_12", 277),          std::make_pair("L1_DoubleIsoTau28er", 109),          std::make_pair("L1_DoubleIsoTau30er", 110),          std::make_pair("L1_DoubleIsoTau32er", 111),          std::make_pair("L1_DoubleIsoTau33er", 264),          std::make_pair("L1_DoubleIsoTau34er", 265),          std::make_pair("L1_DoubleIsoTau35er", 266),          std::make_pair("L1_DoubleIsoTau36er", 278),          std::make_pair("L1_DoubleJet12_ForwardBackward", 213),          std::make_pair("L1_DoubleJet16_ForwardBackward", 214),          std::make_pair("L1_DoubleJet8_ForwardBackward", 212),          std::make_pair("L1_DoubleJetC100", 96),          std::make_pair("L1_DoubleJetC112", 97),          std::make_pair("L1_DoubleJetC120", 98),          std::make_pair("L1_DoubleJetC40", 92),          std::make_pair("L1_DoubleJetC50", 93),          std::make_pair("L1_DoubleJetC60", 94),          std::make_pair("L1_DoubleJetC60_ETM60", 181),          std::make_pair("L1_DoubleJetC80", 95),          std::make_pair("L1_DoubleMu0", 24),          std::make_pair("L1_DoubleMu0_ETM40", 256),          std::make_pair("L1_DoubleMu0_ETM55", 257),          std::make_pair("L1_DoubleMu0er1p4_dEta_Max1p8_OS", 267),          std::make_pair("L1_DoubleMu0er1p6_dEta_Max1p8", 32),          std::make_pair("L1_DoubleMu0er1p6_dEta_Max1p8_OS", 33),          std::make_pair("L1_DoubleMu7_EG14", 166),          std::make_pair("L1_DoubleMu7_EG7", 167),          std::make_pair("L1_DoubleMuOpen", 23),          std::make_pair("L1_DoubleMu_10_0_dEta_Max1p8", 35),          std::make_pair("L1_DoubleMu_10_3p5", 26),          std::make_pair("L1_DoubleMu_10_Open", 25),          std::make_pair("L1_DoubleMu_11_4", 27),          std::make_pair("L1_DoubleMu_12_5", 28),          std::make_pair("L1_DoubleMu_12_8", 31),          std::make_pair("L1_DoubleMu_13_6", 29),          std::make_pair("L1_DoubleMu_15_5", 30),          std::make_pair("L1_DoubleTau50er", 114),          std::make_pair("L1_EG25er_HTT125", 175),          std::make_pair("L1_EG27er_HTT200", 174),          std::make_pair("L1_ETM100", 142),          std::make_pair("L1_ETM120", 143),          std::make_pair("L1_ETM30", 137),          std::make_pair("L1_ETM40", 138),          std::make_pair("L1_ETM50", 139),          std::make_pair("L1_ETM60", 140),          std::make_pair("L1_ETM70", 141),          std::make_pair("L1_ETM75", 272),          std::make_pair("L1_ETM75_Jet60_dPhi_Min0p4", 271),          std::make_pair("L1_ETM80", 125),          std::make_pair("L1_ETM85", 273),          std::make_pair("L1_ETM90", 126),          std::make_pair("L1_ETM95", 274),          std::make_pair("L1_ETT25", 135),          std::make_pair("L1_ETT40_BptxAND", 136),          std::make_pair("L1_ETT50_BptxAND", 244),          std::make_pair("L1_ETT60_BptxAND", 245),          std::make_pair("L1_ETT70_BptxAND", 193),          std::make_pair("L1_FirstBunchAfterTrain", 282),          std::make_pair("L1_FirstBunchInTrain", 281),          std::make_pair("L1_HTM100", 130),          std::make_pair("L1_HTM120", 131),          std::make_pair("L1_HTM130", 132),          std::make_pair("L1_HTM140", 133),          std::make_pair("L1_HTM150", 134),          std::make_pair("L1_HTM50", 127),          std::make_pair("L1_HTM60_HTT260", 188),          std::make_pair("L1_HTM70", 128),          std::make_pair("L1_HTM80", 129),          std::make_pair("L1_HTM80_HTT220", 187),          std::make_pair("L1_HTT120", 115),          std::make_pair("L1_HTT160", 116),          std::make_pair("L1_HTT200", 117),          std::make_pair("L1_HTT220", 118),          std::make_pair("L1_HTT240", 119),          std::make_pair("L1_HTT255", 120),          std::make_pair("L1_HTT270", 121),          std::make_pair("L1_HTT280", 122),          std::make_pair("L1_HTT300", 123),          std::make_pair("L1_HTT320", 124),          std::make_pair("L1_IsoEG18er_IsoTau24er_dEta_Min0p2", 269),          std::make_pair("L1_IsoEG20er_IsoTau25er_dEta_Min0p2", 270),          std::make_pair("L1_IsoEG22er_IsoTau26er_dEta_Min0p2", 268),          std::make_pair("L1_IsoEG22er_Tau20er_dEta_Min0p2", 199),          std::make_pair("L1_IsoEG24_DoubleJetC32", 330),          std::make_pair("L1_IsoEG24_DoubleJetC35", 336),          std::make_pair("L1_IsoEG24_DoubleJetC40", 342),          std::make_pair("L1_IsoEG24_DoubleJetC48", 348),          std::make_pair("L1_IsoEG24_DoubleJetC50", 354),          std::make_pair("L1_IsoEG24_JetC32", 300),          std::make_pair("L1_IsoEG24_JetC35", 306),          std::make_pair("L1_IsoEG24_JetC40", 312),          std::make_pair("L1_IsoEG24_JetC48", 318),          std::make_pair("L1_IsoEG24_JetC50", 324),          std::make_pair("L1_IsoEG24_TripleJetC32", 360),          std::make_pair("L1_IsoEG24_TripleJetC35", 366),          std::make_pair("L1_IsoEG24_TripleJetC35_50_60", 390),          std::make_pair("L1_IsoEG24_TripleJetC40", 372),          std::make_pair("L1_IsoEG24_TripleJetC48", 378),          std::make_pair("L1_IsoEG24_TripleJetC50", 384),          std::make_pair("L1_IsoEG26_DoubleJetC32", 331),          std::make_pair("L1_IsoEG26_DoubleJetC35", 337),          std::make_pair("L1_IsoEG26_DoubleJetC40", 343),          std::make_pair("L1_IsoEG26_DoubleJetC48", 349),          std::make_pair("L1_IsoEG26_DoubleJetC50", 355),          std::make_pair("L1_IsoEG26_JetC32", 301),          std::make_pair("L1_IsoEG26_JetC35", 307),          std::make_pair("L1_IsoEG26_JetC40", 313),          std::make_pair("L1_IsoEG26_JetC48", 319),          std::make_pair("L1_IsoEG26_JetC50", 325),          std::make_pair("L1_IsoEG26_TripleJetC32", 361),          std::make_pair("L1_IsoEG26_TripleJetC35", 367),          std::make_pair("L1_IsoEG26_TripleJetC35_50_60", 391),          std::make_pair("L1_IsoEG26_TripleJetC40", 373),          std::make_pair("L1_IsoEG26_TripleJetC48", 379),          std::make_pair("L1_IsoEG26_TripleJetC50", 385),          std::make_pair("L1_IsoEG28_DoubleJetC32", 332),          std::make_pair("L1_IsoEG28_DoubleJetC35", 338),          std::make_pair("L1_IsoEG28_DoubleJetC40", 344),          std::make_pair("L1_IsoEG28_DoubleJetC48", 350),          std::make_pair("L1_IsoEG28_DoubleJetC50", 356),          std::make_pair("L1_IsoEG28_JetC32", 302),          std::make_pair("L1_IsoEG28_JetC35", 308),          std::make_pair("L1_IsoEG28_JetC40", 314),          std::make_pair("L1_IsoEG28_JetC48", 320),          std::make_pair("L1_IsoEG28_JetC50", 326),          std::make_pair("L1_IsoEG28_TripleJetC32", 362),          std::make_pair("L1_IsoEG28_TripleJetC35", 368),          std::make_pair("L1_IsoEG28_TripleJetC35_50_60", 392),          std::make_pair("L1_IsoEG28_TripleJetC40", 374),          std::make_pair("L1_IsoEG28_TripleJetC48", 380),          std::make_pair("L1_IsoEG28_TripleJetC50", 386),          std::make_pair("L1_IsoEG32_DoubleJetC32", 333),          std::make_pair("L1_IsoEG32_DoubleJetC35", 339),          std::make_pair("L1_IsoEG32_DoubleJetC40", 345),          std::make_pair("L1_IsoEG32_DoubleJetC48", 351),          std::make_pair("L1_IsoEG32_DoubleJetC50", 357),          std::make_pair("L1_IsoEG32_JetC32", 303),          std::make_pair("L1_IsoEG32_JetC35", 309),          std::make_pair("L1_IsoEG32_JetC40", 315),          std::make_pair("L1_IsoEG32_JetC48", 321),          std::make_pair("L1_IsoEG32_JetC50", 327),          std::make_pair("L1_IsoEG32_TripleJetC32", 363),          std::make_pair("L1_IsoEG32_TripleJetC35", 369),          std::make_pair("L1_IsoEG32_TripleJetC35_50_60", 393),          std::make_pair("L1_IsoEG32_TripleJetC40", 375),          std::make_pair("L1_IsoEG32_TripleJetC48", 381),          std::make_pair("L1_IsoEG32_TripleJetC50", 387),          std::make_pair("L1_IsoEG34_DoubleJetC32", 334),          std::make_pair("L1_IsoEG34_DoubleJetC35", 340),          std::make_pair("L1_IsoEG34_DoubleJetC40", 346),          std::make_pair("L1_IsoEG34_DoubleJetC48", 352),          std::make_pair("L1_IsoEG34_DoubleJetC50", 358),          std::make_pair("L1_IsoEG34_JetC32", 304),          std::make_pair("L1_IsoEG34_JetC35", 310),          std::make_pair("L1_IsoEG34_JetC40", 316),          std::make_pair("L1_IsoEG34_JetC48", 322),          std::make_pair("L1_IsoEG34_JetC50", 328),          std::make_pair("L1_IsoEG34_TripleJetC32", 364),          std::make_pair("L1_IsoEG34_TripleJetC35", 370),          std::make_pair("L1_IsoEG34_TripleJetC35_50_60", 394),          std::make_pair("L1_IsoEG34_TripleJetC40", 376),          std::make_pair("L1_IsoEG34_TripleJetC48", 382),          std::make_pair("L1_IsoEG34_TripleJetC50", 388),          std::make_pair("L1_IsoEG38_DoubleJetC32", 335),          std::make_pair("L1_IsoEG38_DoubleJetC35", 341),          std::make_pair("L1_IsoEG38_DoubleJetC40", 347),          std::make_pair("L1_IsoEG38_DoubleJetC48", 353),          std::make_pair("L1_IsoEG38_DoubleJetC50", 359),          std::make_pair("L1_IsoEG38_JetC32", 305),          std::make_pair("L1_IsoEG38_JetC35", 311),          std::make_pair("L1_IsoEG38_JetC40", 317),          std::make_pair("L1_IsoEG38_JetC48", 323),          std::make_pair("L1_IsoEG38_JetC50", 329),          std::make_pair("L1_IsoEG38_TripleJetC32", 365),          std::make_pair("L1_IsoEG38_TripleJetC35", 371),          std::make_pair("L1_IsoEG38_TripleJetC35_50_60", 395),          std::make_pair("L1_IsoEG38_TripleJetC40", 377),          std::make_pair("L1_IsoEG38_TripleJetC48", 383),          std::make_pair("L1_IsoEG38_TripleJetC50", 389),          std::make_pair("L1_IsolatedBunch", 219),          std::make_pair("L1_Jet32_DoubleMu_10_0_dPhi_Jet_Mu0_Max0p4_dPhi_Mu_Mu_Min1p0", 179),          std::make_pair("L1_Jet32_Mu0_EG10_dPhi_Jet_Mu_Max0p4_dPhi_Mu_EG_Min1p0", 180),          std::make_pair("L1_MU20_EG15", 198),          std::make_pair("L1_MinimumBiasHF0_AND", 251),          std::make_pair("L1_MinimumBiasHF0_AND_BptxAND", 247),          std::make_pair("L1_MinimumBiasHF0_OR", 250),          std::make_pair("L1_MinimumBiasHF0_OR_BptxAND", 246),          std::make_pair("L1_MinimumBiasHF1_AND", 253),          std::make_pair("L1_MinimumBiasHF1_AND_BptxAND", 249),          std::make_pair("L1_MinimumBiasHF1_OR", 252),          std::make_pair("L1_MinimumBiasHF1_OR_BptxAND", 248),          std::make_pair("L1_Mu0er_ETM40", 182),          std::make_pair("L1_Mu0er_ETM55", 183),          std::make_pair("L1_Mu10er_ETM30", 184),          std::make_pair("L1_Mu10er_ETM50", 185),          std::make_pair("L1_Mu12_EG10", 149),          std::make_pair("L1_Mu14er_ETM30", 186),          std::make_pair("L1_Mu16er_Tau20er", 154),          std::make_pair("L1_Mu16er_Tau24er", 155),          std::make_pair("L1_Mu18er_IsoTau26er", 158),          std::make_pair("L1_Mu18er_Tau20er", 156),          std::make_pair("L1_Mu18er_Tau24er", 157),          std::make_pair("L1_Mu20_EG10", 150),          std::make_pair("L1_Mu20_EG17", 151),          std::make_pair("L1_Mu20_IsoEG6", 275),          std::make_pair("L1_Mu20er_IsoTau26er", 159),          std::make_pair("L1_Mu22er_IsoTau26er", 279),          std::make_pair("L1_Mu23_EG10", 153),          std::make_pair("L1_Mu23_IsoEG10", 152),          std::make_pair("L1_Mu25er_IsoTau26er", 280),          std::make_pair("L1_Mu3_JetC120", 217),          std::make_pair("L1_Mu3_JetC120_dEta_Max0p4_dPhi_Max0p4", 210),          std::make_pair("L1_Mu3_JetC16", 215),          std::make_pair("L1_Mu3_JetC16_dEta_Max0p4_dPhi_Max0p4", 170),          std::make_pair("L1_Mu3_JetC60", 216),          std::make_pair("L1_Mu3_JetC60_dEta_Max0p4_dPhi_Max0p4", 171),          std::make_pair("L1_Mu5_EG15", 144),          std::make_pair("L1_Mu5_EG20", 145),          std::make_pair("L1_Mu5_EG23", 146),          std::make_pair("L1_Mu5_IsoEG18", 147),          std::make_pair("L1_Mu5_IsoEG20", 148),          std::make_pair("L1_Mu6_DoubleEG10", 169),          std::make_pair("L1_Mu6_DoubleEG17", 168),          std::make_pair("L1_Mu6_HTT200", 172),          std::make_pair("L1_Mu8_HTT150", 173),          std::make_pair("L1_NotBptxOR", 254),          std::make_pair("L1_QuadJetC36_Tau52", 177),          std::make_pair("L1_QuadJetC40", 102),          std::make_pair("L1_QuadJetC50", 103),          std::make_pair("L1_QuadJetC60", 104),          std::make_pair("L1_QuadMu0", 38),          std::make_pair("L1_SingleEG10", 40),          std::make_pair("L1_SingleEG15", 41),          std::make_pair("L1_SingleEG18", 42),          std::make_pair("L1_SingleEG24", 43),          std::make_pair("L1_SingleEG26", 44),          std::make_pair("L1_SingleEG28", 45),          std::make_pair("L1_SingleEG2_BptxAND", 192),          std::make_pair("L1_SingleEG30", 46),          std::make_pair("L1_SingleEG32", 258),          std::make_pair("L1_SingleEG34", 48),          std::make_pair("L1_SingleEG36", 259),          std::make_pair("L1_SingleEG38", 260),          std::make_pair("L1_SingleEG40", 50),          std::make_pair("L1_SingleEG45", 52),          std::make_pair("L1_SingleEG5", 39),          std::make_pair("L1_SingleIsoEG18", 53),          std::make_pair("L1_SingleIsoEG18er", 62),          std::make_pair("L1_SingleIsoEG20", 54),          std::make_pair("L1_SingleIsoEG20er", 63),          std::make_pair("L1_SingleIsoEG22", 55),          std::make_pair("L1_SingleIsoEG22er", 64),          std::make_pair("L1_SingleIsoEG24", 56),          std::make_pair("L1_SingleIsoEG24er", 65),          std::make_pair("L1_SingleIsoEG26", 57),          std::make_pair("L1_SingleIsoEG26er", 66),          std::make_pair("L1_SingleIsoEG28", 59),          std::make_pair("L1_SingleIsoEG28er", 68),          std::make_pair("L1_SingleIsoEG30", 60),          std::make_pair("L1_SingleIsoEG30er", 69),          std::make_pair("L1_SingleIsoEG32", 261),          std::make_pair("L1_SingleIsoEG32er", 263),          std::make_pair("L1_SingleIsoEG34", 61),          std::make_pair("L1_SingleIsoEG34er", 70),          std::make_pair("L1_SingleIsoEG36", 262),          std::make_pair("L1_SingleJet120", 85),          std::make_pair("L1_SingleJet12_BptxAND", 195),          std::make_pair("L1_SingleJet140", 86),          std::make_pair("L1_SingleJet150", 87),          std::make_pair("L1_SingleJet16", 80),          std::make_pair("L1_SingleJet160", 88),          std::make_pair("L1_SingleJet170", 89),          std::make_pair("L1_SingleJet180", 90),          std::make_pair("L1_SingleJet20", 81),          std::make_pair("L1_SingleJet200", 91),          std::make_pair("L1_SingleJet35", 82),          std::make_pair("L1_SingleJet60", 83),          std::make_pair("L1_SingleJet8_BptxAND", 194),          std::make_pair("L1_SingleJet90", 84),          std::make_pair("L1_SingleJetC20_NotBptxOR", 191),          std::make_pair("L1_SingleJetC20_NotBptxOR_3BX", 201),          std::make_pair("L1_SingleJetC32_NotBptxOR", 190),          std::make_pair("L1_SingleJetC32_NotBptxOR_3BX", 202),          std::make_pair("L1_SingleJetC40_NotBptxOR_3BX", 203),          std::make_pair("L1_SingleMu10_LowQ", 14),          std::make_pair("L1_SingleMu12", 6),          std::make_pair("L1_SingleMu14", 7),          std::make_pair("L1_SingleMu14er", 16),          std::make_pair("L1_SingleMu16", 8),          std::make_pair("L1_SingleMu16er", 17),          std::make_pair("L1_SingleMu18", 9),          std::make_pair("L1_SingleMu18er", 18),          std::make_pair("L1_SingleMu20", 10),          std::make_pair("L1_SingleMu20er", 19),          std::make_pair("L1_SingleMu22", 11),          std::make_pair("L1_SingleMu22er", 20),          std::make_pair("L1_SingleMu24_DoubleJetC32", 416),          std::make_pair("L1_SingleMu24_DoubleJetC35", 420),          std::make_pair("L1_SingleMu24_DoubleJetC40", 424),          std::make_pair("L1_SingleMu24_DoubleJetC48", 428),          std::make_pair("L1_SingleMu24_DoubleJetC50", 432),          std::make_pair("L1_SingleMu24_JetC32", 396),          std::make_pair("L1_SingleMu24_JetC35", 400),          std::make_pair("L1_SingleMu24_JetC40", 404),          std::make_pair("L1_SingleMu24_JetC48", 408),          std::make_pair("L1_SingleMu24_JetC50", 412),          std::make_pair("L1_SingleMu24_TripleJetC32", 436),          std::make_pair("L1_SingleMu24_TripleJetC35", 440),          std::make_pair("L1_SingleMu24_TripleJetC35_50_60", 456),          std::make_pair("L1_SingleMu24_TripleJetC40", 444),          std::make_pair("L1_SingleMu24_TripleJetC48", 448),          std::make_pair("L1_SingleMu24_TripleJetC50", 452),          std::make_pair("L1_SingleMu25", 12),          std::make_pair("L1_SingleMu25_DoubleJetC30_", 15),          std::make_pair("L1_SingleMu25er", 21),          std::make_pair("L1_SingleMu28_DoubleJetC32", 417),          std::make_pair("L1_SingleMu28_DoubleJetC35", 421),          std::make_pair("L1_SingleMu28_DoubleJetC40", 425),          std::make_pair("L1_SingleMu28_DoubleJetC48", 429),          std::make_pair("L1_SingleMu28_DoubleJetC50", 433),          std::make_pair("L1_SingleMu28_JetC32", 397),          std::make_pair("L1_SingleMu28_JetC35", 401),          std::make_pair("L1_SingleMu28_JetC40", 405),          std::make_pair("L1_SingleMu28_JetC48", 409),          std::make_pair("L1_SingleMu28_JetC50", 413),          std::make_pair("L1_SingleMu28_TripleJetC32", 437),          std::make_pair("L1_SingleMu28_TripleJetC35", 441),          std::make_pair("L1_SingleMu28_TripleJetC35_50_60", 457),          std::make_pair("L1_SingleMu28_TripleJetC40", 445),          std::make_pair("L1_SingleMu28_TripleJetC48", 449),          std::make_pair("L1_SingleMu28_TripleJetC50", 453),          std::make_pair("L1_SingleMu3", 3),          std::make_pair("L1_SingleMu30", 13),          std::make_pair("L1_SingleMu30er", 22),          std::make_pair("L1_SingleMu32_DoubleJetC32", 418),          std::make_pair("L1_SingleMu32_DoubleJetC35", 422),          std::make_pair("L1_SingleMu32_DoubleJetC40", 426),          std::make_pair("L1_SingleMu32_DoubleJetC48", 430),          std::make_pair("L1_SingleMu32_DoubleJetC50", 434),          std::make_pair("L1_SingleMu32_JetC32", 398),          std::make_pair("L1_SingleMu32_JetC35", 402),          std::make_pair("L1_SingleMu32_JetC40", 406),          std::make_pair("L1_SingleMu32_JetC48", 410),          std::make_pair("L1_SingleMu32_JetC50", 414),          std::make_pair("L1_SingleMu32_TripleJetC32", 438),          std::make_pair("L1_SingleMu32_TripleJetC35", 442),          std::make_pair("L1_SingleMu32_TripleJetC35_50_60", 458),          std::make_pair("L1_SingleMu32_TripleJetC40", 446),          std::make_pair("L1_SingleMu32_TripleJetC48", 450),          std::make_pair("L1_SingleMu32_TripleJetC50", 454),          std::make_pair("L1_SingleMu38_DoubleJetC32", 419),          std::make_pair("L1_SingleMu38_DoubleJetC35", 423),          std::make_pair("L1_SingleMu38_DoubleJetC40", 427),          std::make_pair("L1_SingleMu38_DoubleJetC48", 431),          std::make_pair("L1_SingleMu38_DoubleJetC50", 435),          std::make_pair("L1_SingleMu38_JetC32", 399),          std::make_pair("L1_SingleMu38_JetC35", 403),          std::make_pair("L1_SingleMu38_JetC40", 407),          std::make_pair("L1_SingleMu38_JetC48", 411),          std::make_pair("L1_SingleMu38_JetC50", 415),          std::make_pair("L1_SingleMu38_TripleJetC32", 439),          std::make_pair("L1_SingleMu38_TripleJetC35", 443),          std::make_pair("L1_SingleMu38_TripleJetC35_50_60", 459),          std::make_pair("L1_SingleMu38_TripleJetC40", 447),          std::make_pair("L1_SingleMu38_TripleJetC48", 451),          std::make_pair("L1_SingleMu38_TripleJetC50", 455),          std::make_pair("L1_SingleMu5", 4),          std::make_pair("L1_SingleMu7", 5),          std::make_pair("L1_SingleMuCosmics", 218),          std::make_pair("L1_SingleMuOpen", 2),          std::make_pair("L1_SingleMuOpen_NotBptxOR", 189),          std::make_pair("L1_SingleMuOpen_NotBptxOR_3BX", 200),          std::make_pair("L1_SingleTau100er", 106),          std::make_pair("L1_SingleTau120er", 107),          std::make_pair("L1_SingleTau80er", 105),          std::make_pair("L1_TripleEG_14_10_8", 78),          std::make_pair("L1_TripleEG_18_17_8", 79),          std::make_pair("L1_TripleJet_84_68_48_VBF", 99),          std::make_pair("L1_TripleJet_88_72_56_VBF", 100),          std::make_pair("L1_TripleJet_92_76_64_VBF", 101),          std::make_pair("L1_TripleMu0", 36),          std::make_pair("L1_TripleMu_5_0_0", 276),          std::make_pair("L1_TripleMu_5_5_3", 37),          std::make_pair("L1_ZeroBias", 0),          std::make_pair("L1_ZeroBias_FirstCollidingBunch", 211)      };

  static const std::map<std::string, int> Name2Id(name2id, name2id + sizeof(name2id) / sizeof(name2id[0]));
  const std::map<std::string, int>::const_iterator rc = Name2Id.find(name);
  int id = -1;
  if (rc != Name2Id.end()) id = rc->second;
  return id;
}


AlgorithmFunction getFuncFromId(const int index)
{
  static const std::pair<int, AlgorithmFunction> id2func[] = {
          std::make_pair(206, &L1_AlwaysTrue),          std::make_pair(240, &L1_BRIL_TRIG0_AND),          std::make_pair(243, &L1_BRIL_TRIG0_FstBunchInTrain),          std::make_pair(242, &L1_BRIL_TRIG0_OR),          std::make_pair(241, &L1_BRIL_TRIG0_delayedAND),          std::make_pair(223, &L1_BeamGasB1),          std::make_pair(224, &L1_BeamGasB2),          std::make_pair(222, &L1_BeamGasMinus),          std::make_pair(221, &L1_BeamGasPlus),          std::make_pair(208, &L1_BptxMinus),          std::make_pair(209, &L1_BptxOR),          std::make_pair(207, &L1_BptxPlus),          std::make_pair(220, &L1_BptxXOR),          std::make_pair(176, &L1_DoubleEG6_HTT255),          std::make_pair(71, &L1_DoubleEG_15_10),          std::make_pair(72, &L1_DoubleEG_18_17),          std::make_pair(73, &L1_DoubleEG_20_18),          std::make_pair(75, &L1_DoubleEG_22_10),          std::make_pair(196, &L1_DoubleEG_22_12),          std::make_pair(197, &L1_DoubleEG_22_15),          std::make_pair(76, &L1_DoubleEG_23_10),          std::make_pair(77, &L1_DoubleEG_24_17),          std::make_pair(277, &L1_DoubleEG_25_12),          std::make_pair(109, &L1_DoubleIsoTau28er),          std::make_pair(110, &L1_DoubleIsoTau30er),          std::make_pair(111, &L1_DoubleIsoTau32er),          std::make_pair(264, &L1_DoubleIsoTau33er),          std::make_pair(265, &L1_DoubleIsoTau34er),          std::make_pair(266, &L1_DoubleIsoTau35er),          std::make_pair(278, &L1_DoubleIsoTau36er),          std::make_pair(213, &L1_DoubleJet12_ForwardBackward),          std::make_pair(214, &L1_DoubleJet16_ForwardBackward),          std::make_pair(212, &L1_DoubleJet8_ForwardBackward),          std::make_pair(96, &L1_DoubleJetC100),          std::make_pair(97, &L1_DoubleJetC112),          std::make_pair(98, &L1_DoubleJetC120),          std::make_pair(92, &L1_DoubleJetC40),          std::make_pair(93, &L1_DoubleJetC50),          std::make_pair(94, &L1_DoubleJetC60),          std::make_pair(181, &L1_DoubleJetC60_ETM60),          std::make_pair(95, &L1_DoubleJetC80),          std::make_pair(24, &L1_DoubleMu0),          std::make_pair(256, &L1_DoubleMu0_ETM40),          std::make_pair(257, &L1_DoubleMu0_ETM55),          std::make_pair(267, &L1_DoubleMu0er1p4_dEta_Max1p8_OS),          std::make_pair(32, &L1_DoubleMu0er1p6_dEta_Max1p8),          std::make_pair(33, &L1_DoubleMu0er1p6_dEta_Max1p8_OS),          std::make_pair(166, &L1_DoubleMu7_EG14),          std::make_pair(167, &L1_DoubleMu7_EG7),          std::make_pair(23, &L1_DoubleMuOpen),          std::make_pair(35, &L1_DoubleMu_10_0_dEta_Max1p8),          std::make_pair(26, &L1_DoubleMu_10_3p5),          std::make_pair(25, &L1_DoubleMu_10_Open),          std::make_pair(27, &L1_DoubleMu_11_4),          std::make_pair(28, &L1_DoubleMu_12_5),          std::make_pair(31, &L1_DoubleMu_12_8),          std::make_pair(29, &L1_DoubleMu_13_6),          std::make_pair(30, &L1_DoubleMu_15_5),          std::make_pair(114, &L1_DoubleTau50er),          std::make_pair(175, &L1_EG25er_HTT125),          std::make_pair(174, &L1_EG27er_HTT200),          std::make_pair(142, &L1_ETM100),          std::make_pair(143, &L1_ETM120),          std::make_pair(137, &L1_ETM30),          std::make_pair(138, &L1_ETM40),          std::make_pair(139, &L1_ETM50),          std::make_pair(140, &L1_ETM60),          std::make_pair(141, &L1_ETM70),          std::make_pair(272, &L1_ETM75),          std::make_pair(271, &L1_ETM75_Jet60_dPhi_Min0p4),          std::make_pair(125, &L1_ETM80),          std::make_pair(273, &L1_ETM85),          std::make_pair(126, &L1_ETM90),          std::make_pair(274, &L1_ETM95),          std::make_pair(135, &L1_ETT25),          std::make_pair(136, &L1_ETT40_BptxAND),          std::make_pair(244, &L1_ETT50_BptxAND),          std::make_pair(245, &L1_ETT60_BptxAND),          std::make_pair(193, &L1_ETT70_BptxAND),          std::make_pair(282, &L1_FirstBunchAfterTrain),          std::make_pair(281, &L1_FirstBunchInTrain),          std::make_pair(130, &L1_HTM100),          std::make_pair(131, &L1_HTM120),          std::make_pair(132, &L1_HTM130),          std::make_pair(133, &L1_HTM140),          std::make_pair(134, &L1_HTM150),          std::make_pair(127, &L1_HTM50),          std::make_pair(188, &L1_HTM60_HTT260),          std::make_pair(128, &L1_HTM70),          std::make_pair(129, &L1_HTM80),          std::make_pair(187, &L1_HTM80_HTT220),          std::make_pair(115, &L1_HTT120),          std::make_pair(116, &L1_HTT160),          std::make_pair(117, &L1_HTT200),          std::make_pair(118, &L1_HTT220),          std::make_pair(119, &L1_HTT240),          std::make_pair(120, &L1_HTT255),          std::make_pair(121, &L1_HTT270),          std::make_pair(122, &L1_HTT280),          std::make_pair(123, &L1_HTT300),          std::make_pair(124, &L1_HTT320),          std::make_pair(269, &L1_IsoEG18er_IsoTau24er_dEta_Min0p2),          std::make_pair(270, &L1_IsoEG20er_IsoTau25er_dEta_Min0p2),          std::make_pair(268, &L1_IsoEG22er_IsoTau26er_dEta_Min0p2),          std::make_pair(199, &L1_IsoEG22er_Tau20er_dEta_Min0p2),          std::make_pair(330, &L1_IsoEG24_DoubleJetC32),          std::make_pair(336, &L1_IsoEG24_DoubleJetC35),          std::make_pair(342, &L1_IsoEG24_DoubleJetC40),          std::make_pair(348, &L1_IsoEG24_DoubleJetC48),          std::make_pair(354, &L1_IsoEG24_DoubleJetC50),          std::make_pair(300, &L1_IsoEG24_JetC32),          std::make_pair(306, &L1_IsoEG24_JetC35),          std::make_pair(312, &L1_IsoEG24_JetC40),          std::make_pair(318, &L1_IsoEG24_JetC48),          std::make_pair(324, &L1_IsoEG24_JetC50),          std::make_pair(360, &L1_IsoEG24_TripleJetC32),          std::make_pair(366, &L1_IsoEG24_TripleJetC35),          std::make_pair(390, &L1_IsoEG24_TripleJetC35_50_60),          std::make_pair(372, &L1_IsoEG24_TripleJetC40),          std::make_pair(378, &L1_IsoEG24_TripleJetC48),          std::make_pair(384, &L1_IsoEG24_TripleJetC50),          std::make_pair(331, &L1_IsoEG26_DoubleJetC32),          std::make_pair(337, &L1_IsoEG26_DoubleJetC35),          std::make_pair(343, &L1_IsoEG26_DoubleJetC40),          std::make_pair(349, &L1_IsoEG26_DoubleJetC48),          std::make_pair(355, &L1_IsoEG26_DoubleJetC50),          std::make_pair(301, &L1_IsoEG26_JetC32),          std::make_pair(307, &L1_IsoEG26_JetC35),          std::make_pair(313, &L1_IsoEG26_JetC40),          std::make_pair(319, &L1_IsoEG26_JetC48),          std::make_pair(325, &L1_IsoEG26_JetC50),          std::make_pair(361, &L1_IsoEG26_TripleJetC32),          std::make_pair(367, &L1_IsoEG26_TripleJetC35),          std::make_pair(391, &L1_IsoEG26_TripleJetC35_50_60),          std::make_pair(373, &L1_IsoEG26_TripleJetC40),          std::make_pair(379, &L1_IsoEG26_TripleJetC48),          std::make_pair(385, &L1_IsoEG26_TripleJetC50),          std::make_pair(332, &L1_IsoEG28_DoubleJetC32),          std::make_pair(338, &L1_IsoEG28_DoubleJetC35),          std::make_pair(344, &L1_IsoEG28_DoubleJetC40),          std::make_pair(350, &L1_IsoEG28_DoubleJetC48),          std::make_pair(356, &L1_IsoEG28_DoubleJetC50),          std::make_pair(302, &L1_IsoEG28_JetC32),          std::make_pair(308, &L1_IsoEG28_JetC35),          std::make_pair(314, &L1_IsoEG28_JetC40),          std::make_pair(320, &L1_IsoEG28_JetC48),          std::make_pair(326, &L1_IsoEG28_JetC50),          std::make_pair(362, &L1_IsoEG28_TripleJetC32),          std::make_pair(368, &L1_IsoEG28_TripleJetC35),          std::make_pair(392, &L1_IsoEG28_TripleJetC35_50_60),          std::make_pair(374, &L1_IsoEG28_TripleJetC40),          std::make_pair(380, &L1_IsoEG28_TripleJetC48),          std::make_pair(386, &L1_IsoEG28_TripleJetC50),          std::make_pair(333, &L1_IsoEG32_DoubleJetC32),          std::make_pair(339, &L1_IsoEG32_DoubleJetC35),          std::make_pair(345, &L1_IsoEG32_DoubleJetC40),          std::make_pair(351, &L1_IsoEG32_DoubleJetC48),          std::make_pair(357, &L1_IsoEG32_DoubleJetC50),          std::make_pair(303, &L1_IsoEG32_JetC32),          std::make_pair(309, &L1_IsoEG32_JetC35),          std::make_pair(315, &L1_IsoEG32_JetC40),          std::make_pair(321, &L1_IsoEG32_JetC48),          std::make_pair(327, &L1_IsoEG32_JetC50),          std::make_pair(363, &L1_IsoEG32_TripleJetC32),          std::make_pair(369, &L1_IsoEG32_TripleJetC35),          std::make_pair(393, &L1_IsoEG32_TripleJetC35_50_60),          std::make_pair(375, &L1_IsoEG32_TripleJetC40),          std::make_pair(381, &L1_IsoEG32_TripleJetC48),          std::make_pair(387, &L1_IsoEG32_TripleJetC50),          std::make_pair(334, &L1_IsoEG34_DoubleJetC32),          std::make_pair(340, &L1_IsoEG34_DoubleJetC35),          std::make_pair(346, &L1_IsoEG34_DoubleJetC40),          std::make_pair(352, &L1_IsoEG34_DoubleJetC48),          std::make_pair(358, &L1_IsoEG34_DoubleJetC50),          std::make_pair(304, &L1_IsoEG34_JetC32),          std::make_pair(310, &L1_IsoEG34_JetC35),          std::make_pair(316, &L1_IsoEG34_JetC40),          std::make_pair(322, &L1_IsoEG34_JetC48),          std::make_pair(328, &L1_IsoEG34_JetC50),          std::make_pair(364, &L1_IsoEG34_TripleJetC32),          std::make_pair(370, &L1_IsoEG34_TripleJetC35),          std::make_pair(394, &L1_IsoEG34_TripleJetC35_50_60),          std::make_pair(376, &L1_IsoEG34_TripleJetC40),          std::make_pair(382, &L1_IsoEG34_TripleJetC48),          std::make_pair(388, &L1_IsoEG34_TripleJetC50),          std::make_pair(335, &L1_IsoEG38_DoubleJetC32),          std::make_pair(341, &L1_IsoEG38_DoubleJetC35),          std::make_pair(347, &L1_IsoEG38_DoubleJetC40),          std::make_pair(353, &L1_IsoEG38_DoubleJetC48),          std::make_pair(359, &L1_IsoEG38_DoubleJetC50),          std::make_pair(305, &L1_IsoEG38_JetC32),          std::make_pair(311, &L1_IsoEG38_JetC35),          std::make_pair(317, &L1_IsoEG38_JetC40),          std::make_pair(323, &L1_IsoEG38_JetC48),          std::make_pair(329, &L1_IsoEG38_JetC50),          std::make_pair(365, &L1_IsoEG38_TripleJetC32),          std::make_pair(371, &L1_IsoEG38_TripleJetC35),          std::make_pair(395, &L1_IsoEG38_TripleJetC35_50_60),          std::make_pair(377, &L1_IsoEG38_TripleJetC40),          std::make_pair(383, &L1_IsoEG38_TripleJetC48),          std::make_pair(389, &L1_IsoEG38_TripleJetC50),          std::make_pair(219, &L1_IsolatedBunch),          std::make_pair(179, &L1_Jet32_DoubleMu_10_0_dPhi_Jet_Mu0_Max0p4_dPhi_Mu_Mu_Min1p0),          std::make_pair(180, &L1_Jet32_Mu0_EG10_dPhi_Jet_Mu_Max0p4_dPhi_Mu_EG_Min1p0),          std::make_pair(198, &L1_MU20_EG15),          std::make_pair(251, &L1_MinimumBiasHF0_AND),          std::make_pair(247, &L1_MinimumBiasHF0_AND_BptxAND),          std::make_pair(250, &L1_MinimumBiasHF0_OR),          std::make_pair(246, &L1_MinimumBiasHF0_OR_BptxAND),          std::make_pair(253, &L1_MinimumBiasHF1_AND),          std::make_pair(249, &L1_MinimumBiasHF1_AND_BptxAND),          std::make_pair(252, &L1_MinimumBiasHF1_OR),          std::make_pair(248, &L1_MinimumBiasHF1_OR_BptxAND),          std::make_pair(182, &L1_Mu0er_ETM40),          std::make_pair(183, &L1_Mu0er_ETM55),          std::make_pair(184, &L1_Mu10er_ETM30),          std::make_pair(185, &L1_Mu10er_ETM50),          std::make_pair(149, &L1_Mu12_EG10),          std::make_pair(186, &L1_Mu14er_ETM30),          std::make_pair(154, &L1_Mu16er_Tau20er),          std::make_pair(155, &L1_Mu16er_Tau24er),          std::make_pair(158, &L1_Mu18er_IsoTau26er),          std::make_pair(156, &L1_Mu18er_Tau20er),          std::make_pair(157, &L1_Mu18er_Tau24er),          std::make_pair(150, &L1_Mu20_EG10),          std::make_pair(151, &L1_Mu20_EG17),          std::make_pair(275, &L1_Mu20_IsoEG6),          std::make_pair(159, &L1_Mu20er_IsoTau26er),          std::make_pair(279, &L1_Mu22er_IsoTau26er),          std::make_pair(153, &L1_Mu23_EG10),          std::make_pair(152, &L1_Mu23_IsoEG10),          std::make_pair(280, &L1_Mu25er_IsoTau26er),          std::make_pair(217, &L1_Mu3_JetC120),          std::make_pair(210, &L1_Mu3_JetC120_dEta_Max0p4_dPhi_Max0p4),          std::make_pair(215, &L1_Mu3_JetC16),          std::make_pair(170, &L1_Mu3_JetC16_dEta_Max0p4_dPhi_Max0p4),          std::make_pair(216, &L1_Mu3_JetC60),          std::make_pair(171, &L1_Mu3_JetC60_dEta_Max0p4_dPhi_Max0p4),          std::make_pair(144, &L1_Mu5_EG15),          std::make_pair(145, &L1_Mu5_EG20),          std::make_pair(146, &L1_Mu5_EG23),          std::make_pair(147, &L1_Mu5_IsoEG18),          std::make_pair(148, &L1_Mu5_IsoEG20),          std::make_pair(169, &L1_Mu6_DoubleEG10),          std::make_pair(168, &L1_Mu6_DoubleEG17),          std::make_pair(172, &L1_Mu6_HTT200),          std::make_pair(173, &L1_Mu8_HTT150),          std::make_pair(254, &L1_NotBptxOR),          std::make_pair(177, &L1_QuadJetC36_Tau52),          std::make_pair(102, &L1_QuadJetC40),          std::make_pair(103, &L1_QuadJetC50),          std::make_pair(104, &L1_QuadJetC60),          std::make_pair(38, &L1_QuadMu0),          std::make_pair(40, &L1_SingleEG10),          std::make_pair(41, &L1_SingleEG15),          std::make_pair(42, &L1_SingleEG18),          std::make_pair(43, &L1_SingleEG24),          std::make_pair(44, &L1_SingleEG26),          std::make_pair(45, &L1_SingleEG28),          std::make_pair(192, &L1_SingleEG2_BptxAND),          std::make_pair(46, &L1_SingleEG30),          std::make_pair(258, &L1_SingleEG32),          std::make_pair(48, &L1_SingleEG34),          std::make_pair(259, &L1_SingleEG36),          std::make_pair(260, &L1_SingleEG38),          std::make_pair(50, &L1_SingleEG40),          std::make_pair(52, &L1_SingleEG45),          std::make_pair(39, &L1_SingleEG5),          std::make_pair(53, &L1_SingleIsoEG18),          std::make_pair(62, &L1_SingleIsoEG18er),          std::make_pair(54, &L1_SingleIsoEG20),          std::make_pair(63, &L1_SingleIsoEG20er),          std::make_pair(55, &L1_SingleIsoEG22),          std::make_pair(64, &L1_SingleIsoEG22er),          std::make_pair(56, &L1_SingleIsoEG24),          std::make_pair(65, &L1_SingleIsoEG24er),          std::make_pair(57, &L1_SingleIsoEG26),          std::make_pair(66, &L1_SingleIsoEG26er),          std::make_pair(59, &L1_SingleIsoEG28),          std::make_pair(68, &L1_SingleIsoEG28er),          std::make_pair(60, &L1_SingleIsoEG30),          std::make_pair(69, &L1_SingleIsoEG30er),          std::make_pair(261, &L1_SingleIsoEG32),          std::make_pair(263, &L1_SingleIsoEG32er),          std::make_pair(61, &L1_SingleIsoEG34),          std::make_pair(70, &L1_SingleIsoEG34er),          std::make_pair(262, &L1_SingleIsoEG36),          std::make_pair(85, &L1_SingleJet120),          std::make_pair(195, &L1_SingleJet12_BptxAND),          std::make_pair(86, &L1_SingleJet140),          std::make_pair(87, &L1_SingleJet150),          std::make_pair(80, &L1_SingleJet16),          std::make_pair(88, &L1_SingleJet160),          std::make_pair(89, &L1_SingleJet170),          std::make_pair(90, &L1_SingleJet180),          std::make_pair(81, &L1_SingleJet20),          std::make_pair(91, &L1_SingleJet200),          std::make_pair(82, &L1_SingleJet35),          std::make_pair(83, &L1_SingleJet60),          std::make_pair(194, &L1_SingleJet8_BptxAND),          std::make_pair(84, &L1_SingleJet90),          std::make_pair(191, &L1_SingleJetC20_NotBptxOR),          std::make_pair(201, &L1_SingleJetC20_NotBptxOR_3BX),          std::make_pair(190, &L1_SingleJetC32_NotBptxOR),          std::make_pair(202, &L1_SingleJetC32_NotBptxOR_3BX),          std::make_pair(203, &L1_SingleJetC40_NotBptxOR_3BX),          std::make_pair(14, &L1_SingleMu10_LowQ),          std::make_pair(6, &L1_SingleMu12),          std::make_pair(7, &L1_SingleMu14),          std::make_pair(16, &L1_SingleMu14er),          std::make_pair(8, &L1_SingleMu16),          std::make_pair(17, &L1_SingleMu16er),          std::make_pair(9, &L1_SingleMu18),          std::make_pair(18, &L1_SingleMu18er),          std::make_pair(10, &L1_SingleMu20),          std::make_pair(19, &L1_SingleMu20er),          std::make_pair(11, &L1_SingleMu22),          std::make_pair(20, &L1_SingleMu22er),          std::make_pair(416, &L1_SingleMu24_DoubleJetC32),          std::make_pair(420, &L1_SingleMu24_DoubleJetC35),          std::make_pair(424, &L1_SingleMu24_DoubleJetC40),          std::make_pair(428, &L1_SingleMu24_DoubleJetC48),          std::make_pair(432, &L1_SingleMu24_DoubleJetC50),          std::make_pair(396, &L1_SingleMu24_JetC32),          std::make_pair(400, &L1_SingleMu24_JetC35),          std::make_pair(404, &L1_SingleMu24_JetC40),          std::make_pair(408, &L1_SingleMu24_JetC48),          std::make_pair(412, &L1_SingleMu24_JetC50),          std::make_pair(436, &L1_SingleMu24_TripleJetC32),          std::make_pair(440, &L1_SingleMu24_TripleJetC35),          std::make_pair(456, &L1_SingleMu24_TripleJetC35_50_60),          std::make_pair(444, &L1_SingleMu24_TripleJetC40),          std::make_pair(448, &L1_SingleMu24_TripleJetC48),          std::make_pair(452, &L1_SingleMu24_TripleJetC50),          std::make_pair(12, &L1_SingleMu25),          std::make_pair(15, &L1_SingleMu25_DoubleJetC30_),          std::make_pair(21, &L1_SingleMu25er),          std::make_pair(417, &L1_SingleMu28_DoubleJetC32),          std::make_pair(421, &L1_SingleMu28_DoubleJetC35),          std::make_pair(425, &L1_SingleMu28_DoubleJetC40),          std::make_pair(429, &L1_SingleMu28_DoubleJetC48),          std::make_pair(433, &L1_SingleMu28_DoubleJetC50),          std::make_pair(397, &L1_SingleMu28_JetC32),          std::make_pair(401, &L1_SingleMu28_JetC35),          std::make_pair(405, &L1_SingleMu28_JetC40),          std::make_pair(409, &L1_SingleMu28_JetC48),          std::make_pair(413, &L1_SingleMu28_JetC50),          std::make_pair(437, &L1_SingleMu28_TripleJetC32),          std::make_pair(441, &L1_SingleMu28_TripleJetC35),          std::make_pair(457, &L1_SingleMu28_TripleJetC35_50_60),          std::make_pair(445, &L1_SingleMu28_TripleJetC40),          std::make_pair(449, &L1_SingleMu28_TripleJetC48),          std::make_pair(453, &L1_SingleMu28_TripleJetC50),          std::make_pair(3, &L1_SingleMu3),          std::make_pair(13, &L1_SingleMu30),          std::make_pair(22, &L1_SingleMu30er),          std::make_pair(418, &L1_SingleMu32_DoubleJetC32),          std::make_pair(422, &L1_SingleMu32_DoubleJetC35),          std::make_pair(426, &L1_SingleMu32_DoubleJetC40),          std::make_pair(430, &L1_SingleMu32_DoubleJetC48),          std::make_pair(434, &L1_SingleMu32_DoubleJetC50),          std::make_pair(398, &L1_SingleMu32_JetC32),          std::make_pair(402, &L1_SingleMu32_JetC35),          std::make_pair(406, &L1_SingleMu32_JetC40),          std::make_pair(410, &L1_SingleMu32_JetC48),          std::make_pair(414, &L1_SingleMu32_JetC50),          std::make_pair(438, &L1_SingleMu32_TripleJetC32),          std::make_pair(442, &L1_SingleMu32_TripleJetC35),          std::make_pair(458, &L1_SingleMu32_TripleJetC35_50_60),          std::make_pair(446, &L1_SingleMu32_TripleJetC40),          std::make_pair(450, &L1_SingleMu32_TripleJetC48),          std::make_pair(454, &L1_SingleMu32_TripleJetC50),          std::make_pair(419, &L1_SingleMu38_DoubleJetC32),          std::make_pair(423, &L1_SingleMu38_DoubleJetC35),          std::make_pair(427, &L1_SingleMu38_DoubleJetC40),          std::make_pair(431, &L1_SingleMu38_DoubleJetC48),          std::make_pair(435, &L1_SingleMu38_DoubleJetC50),          std::make_pair(399, &L1_SingleMu38_JetC32),          std::make_pair(403, &L1_SingleMu38_JetC35),          std::make_pair(407, &L1_SingleMu38_JetC40),          std::make_pair(411, &L1_SingleMu38_JetC48),          std::make_pair(415, &L1_SingleMu38_JetC50),          std::make_pair(439, &L1_SingleMu38_TripleJetC32),          std::make_pair(443, &L1_SingleMu38_TripleJetC35),          std::make_pair(459, &L1_SingleMu38_TripleJetC35_50_60),          std::make_pair(447, &L1_SingleMu38_TripleJetC40),          std::make_pair(451, &L1_SingleMu38_TripleJetC48),          std::make_pair(455, &L1_SingleMu38_TripleJetC50),          std::make_pair(4, &L1_SingleMu5),          std::make_pair(5, &L1_SingleMu7),          std::make_pair(218, &L1_SingleMuCosmics),          std::make_pair(2, &L1_SingleMuOpen),          std::make_pair(189, &L1_SingleMuOpen_NotBptxOR),          std::make_pair(200, &L1_SingleMuOpen_NotBptxOR_3BX),          std::make_pair(106, &L1_SingleTau100er),          std::make_pair(107, &L1_SingleTau120er),          std::make_pair(105, &L1_SingleTau80er),          std::make_pair(78, &L1_TripleEG_14_10_8),          std::make_pair(79, &L1_TripleEG_18_17_8),          std::make_pair(99, &L1_TripleJet_84_68_48_VBF),          std::make_pair(100, &L1_TripleJet_88_72_56_VBF),          std::make_pair(101, &L1_TripleJet_92_76_64_VBF),          std::make_pair(36, &L1_TripleMu0),          std::make_pair(276, &L1_TripleMu_5_0_0),          std::make_pair(37, &L1_TripleMu_5_5_3),          std::make_pair(0, &L1_ZeroBias),          std::make_pair(211, &L1_ZeroBias_FirstCollidingBunch)      };

  static const std::map<int, AlgorithmFunction> Id2Func(id2func, id2func + sizeof(id2func) / sizeof(id2func[0]));
  const std::map<int, AlgorithmFunction>::const_iterator rc = Id2Func.find(index);
  AlgorithmFunction fp = 0;
  if (rc != Id2Func.end()) fp = rc->second;
  return fp;
}

// eof