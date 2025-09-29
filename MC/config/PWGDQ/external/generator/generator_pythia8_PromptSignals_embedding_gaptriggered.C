#include "FairGenerator.h"
#include "Generators/GeneratorPythia8.h"
#include "Pythia8/Pythia.h"
#include "TRandom.h"
#include "TParticle.h" 

R__ADD_INCLUDE_PATH($O2DPG_MC_CONFIG_ROOT/MC/config/PWGDQ/EvtGen)
#include "GeneratorEvtGen.C"

#include <string>

using namespace o2::eventgen;

namespace o2
{
namespace eventgen
{

class GeneratorPythia8PromptSignalsEmbeddingGapTriggered : public o2::eventgen::GeneratorPythia8 {
public:
  
  /// constructor
  GeneratorPythia8PromptSignalsEmbeddingGapTriggered(int inputTriggerRatio = 5)  {

    mGeneratedEvents = 0;
    mInverseTriggerRatio = inputTriggerRatio;
    // define minimum bias event generator
    auto seed = (gRandom->TRandom::GetSeed() % 900000000);
    TString pathconfigMB = gSystem->ExpandPathName("${O2DPG_MC_CONFIG_ROOT}/MC/config/PWGDQ/pythia8/generator/pythia8_inel_triggerGap.cfg");
    pythiaMBgen.readFile(pathconfigMB.Data());
    pythiaMBgen.readString("Random:setSeed on");
    pythiaMBgen.readString("Random:seed " + std::to_string(seed));
    mConfigMBdecays = "";
    mRapidityMin = -1.;
    mRapidityMax = 1.;
    mVerbose = false; 
  }

  ///  Destructor
  ~GeneratorPythia8PromptSignalsEmbeddingGapTriggered() = default;

  void addHadronPDGs(int pdg) { mHadronsPDGs.push_back(pdg); };

  void setRapidityRange(double valMin, double valMax)
  {
    mRapidityMin = valMin;
    mRapidityMax = valMax;
  };

  void setTriggerGap(int triggerGap) {mInverseTriggerRatio = triggerGap;}

  void setConfigMBdecays(TString val){mConfigMBdecays = val;}

  void setVerbose(bool val) { mVerbose = val; };

protected:

bool generateEvent() override {
  return true;
}

bool Init() override {
        
  if(mConfigMBdecays.Contains("cfg")) {
    pythiaMBgen.readFile(mConfigMBdecays.Data());	
  }
  addSubGenerator(0, "Minimum bias");
  addSubGenerator(1, "Onia injected");
	GeneratorPythia8::Init();
  pythiaMBgen.init();

  return true;
} 

// search for the presence of at least one of the required hadrons in a selected rapidity window
bool findHadrons(Pythia8::Event& event) { 
   
  for (int ipa = 0; ipa < event.size(); ++ipa) {
    
    auto daughterList = event[ipa].daughterList();
  
    for (auto ida : daughterList) {
      for (int pdg : mHadronsPDGs) {   // check that at least one of the pdg code is found in the event
        if (event[ida].id() == pdg) {
          if ((event[ida].y() > mRapidityMin) && (event[ida].y() < mRapidityMax)) {
            cout << "============= Found jpsi y,pt " <<  event[ida].y() << ", " << event[ida].pT() << endl;
            std::vector<int> daughters = event[ida].daughterList();
            for (int d : daughters) {
              cout << "###### daughter " << d << ": code " << event[d].id() << ", pt " << event[d].pT() << endl;
            }
            return true;
          }
        }
      }
    }
  }

  return false;
};

std::vector<int> findAllCharmonia(const Pythia8::Event& event) {
  std::vector<int> out; out.reserve(4);

  for (int ipa = 0; ipa < event.size(); ++ipa) {
    
    auto daughterList = event[ipa].daughterList();
  
    for (auto ida : daughterList) {
      for (int pdg : mHadronsPDGs) {   // check that at least one of the pdg code is found in the event
        if (event[ida].id() == pdg) {
          if ((event[ida].y() > mRapidityMin) && (event[ida].y() < mRapidityMax)) {
            cout << "============= Found jpsi y,pt " <<  event[ida].y() << ", " << event[ida].pT() << endl;
            out.push_back(ida);
          }
        }
      }
    }
  }
  
  return out;
}

void collectDaughters(const Pythia8::Event& event, int idx, std::vector<int>& keep, std::vector<char>& visited) {
  if (idx < 0 || idx >= event.size()) return;
  if (visited[idx]) return;
  visited[idx] = 1;
  
  keep.push_back(idx);

  int daughter1 =  event[idx].daughter1();
  int daughter2 =  event[idx].daughter2();
  if (daughter1<0 || daughter2 < daughter1) return;
  for (int d = daughter1; d <= daughter2; ++d) {
    collectDaughters(event, d, keep, visited);
  }
}

TParticle makeTParticleTemp(const Pythia8::Event& event, int idx) {
  const auto& q = event[idx];
  int status = q.status();
  if (status < 0) {
    return TParticle(0, 0, -1, -1, -1, -1,
                     0.,0.,0.,0., 0.,0.,0.,0.);
    }

  int m1 = q.mother1();
  int m2 = q.mother2();
  int d1 = q.daughter1();
  int d2 = q.daughter2();
  return TParticle(q.id(), status, m1, m2, d1, d2,
                   q.px(), q.py(), q.pz(), q.e(),
                   q.xProd(), q.yProd(), q.zProd(), q.tProd());
}

Bool_t importParticles() override
{
  LOG(info) << "";
  LOG(info) << "*************************************************************";
  LOG(info) << "************** New signal event considered **************";
  LOG(info) << "*************************************************************";
  LOG(info) << "";

  const int nSig = std::max(1, (int)std::lround(mNumSigEvs));
  for (int isig=0; isig<nSig; ++isig) {

    bool genOk=false; 
    std::vector<int> charmonia;
    while (! (genOk && !charmonia.empty())) {
      /// reset event
      mPythia.event.reset();
      genOk = GeneratorPythia8::generateEvent(); 
      if (!genOk) continue;
      charmonia = findAllCharmonia(mPythia.event);
    }

    std::vector<int> decayChains;
    std::vector<char> visited(mPythia.event.size(), 0);
    decayChains.reserve(128);
    for (int c : charmonia) collectDaughters(mPythia.event, c, decayChains, visited);
    std::sort(decayChains.begin(), decayChains.end());
    decayChains.erase(std::unique(decayChains.begin(), decayChains.end()), decayChains.end());

    const int offset = (int)mParticles.size();
    mParticles.reserve(offset + (int)decayChains.size());
    
    std::vector<int> idxMap(mPythia.event.size(), -1);
    for (int i = 0; i < (int)decayChains.size(); ++i) {
      if (decayChains[i] < 0 || decayChains[i] >= mPythia.event.size()) continue; 
      idxMap[decayChains[i]] = offset + i;
    }

    for (int idx : decayChains){
      TParticle part = makeTParticleTemp(mPythia.event, idx);
      if (part.GetPdgCode() != 0) {
        mParticles.push_back(part);
      }
    }

    for (int iLoc = 0; iLoc < (int)decayChains.size(); ++iLoc) {
      const int srcIdx = decayChains[iLoc];
      const auto& src = mPythia.event[srcIdx];

      const int mother1 = (src.mother1()  >= 0 ? idxMap[src.mother1()]  : -1);
      const int mother2 = (src.mother2()  >= 0 ? idxMap[src.mother2()]  : -1);
      const int daughter1 = (src.daughter1()>= 0 ? idxMap[src.daughter1()] : -1);
      const int daughter2 = (src.daughter2()>= 0 ? idxMap[src.daughter2()] : -1);

      // update TParticle
      TParticle& particle = mParticles[offset + iLoc];
      particle.SetFirstMother(mother1);
      particle.SetLastMother (mother2);
      particle.SetFirstDaughter(daughter1);
      particle.SetLastDaughter (daughter2);
    }
    LOG(info) << "-----------------------------------------------";
    LOG(info) << "============ After event " << isig << " (size " << decayChains.size() << ")";
    LOG(info) << "Full stack (size " << mParticles.size() << "):";
    LOG(info) << "-----------------------------------------------";
    // printParticleVector(mParticles);

    notifySubGenerator(1);
  }

  if (mVerbose) mOutputEvent.list();
  
  return kTRUE;
}

void notifyEmbedding(const o2::dataformats::MCEventHeader* bkgHeader) override {
  LOG(info) << "[notifyEmbedding] ----- Function called";
  
  /// Impact parameter between the two nuclei
  const float x = bkgHeader->GetB();
  LOG(info) << "[notifyEmbedding] ----- Collision impact parameter: " << x;

  /// number of events to be embedded in a background event
  mNumSigEvs = std::max(1.,120.*(x<5.)+80.*(1.-x/20.)*(x>5.)*(x<11.)+240.*(1.-x/13.)*(x>11.));
  LOG(info) << "[notifyEmbedding] ----- generating " << mNumSigEvs << " signal events " << std::endl;
};

private:
  // Interface to override import particles
  Pythia8::Event mOutputEvent;

  // Control gap-triggering
  unsigned long long mGeneratedEvents;
  int mInverseTriggerRatio;
  Pythia8::Pythia pythiaMBgen; // minimum bias event  
  TString mConfigMBdecays;		
  std::vector<int> mHadronsPDGs;
  double mRapidityMin; 
  double mRapidityMax;
  bool mVerbose;

  // number of signal events to be embedded in a background event
  int mNumSigEvs{1};
};    

}

}

// Predefined generators:
FairGenerator*
  GeneratorPromptJpsi_EvtGenMidY(int triggerGap, double rapidityMin = -1.5, double rapidityMax = 1.5, bool verbose = false, bool embedding = false)
{
  auto gen = new o2::eventgen::GeneratorEvtGen<o2::eventgen::GeneratorPythia8PromptSignalsEmbeddingGapTriggered>();
  gen->setTriggerGap(triggerGap);
  gen->setRapidityRange(rapidityMin, rapidityMax);
  gen->addHadronPDGs(443);
  gen->setVerbose(verbose);

  TString pathO2table = gSystem->ExpandPathName("${O2DPG_MC_CONFIG_ROOT}/MC/config/PWGDQ/pythia8/decayer/switchOffJpsi.cfg");
  gen->readFile(pathO2table.Data());
  gen->setConfigMBdecays(pathO2table);
  gen->PrintDebug(true);

  gen->SetSizePdg(1);
  gen->AddPdg(443, 0);

  gen->SetForceDecay(kEvtDiElectron);

  // set random seed
  gen->readString("Random:setSeed on");
  uint random_seed;
  unsigned long long int random_value = 0;
  ifstream urandom("/dev/urandom", ios::in|ios::binary);
  urandom.read(reinterpret_cast<char*>(&random_value), sizeof(random_seed));
  gen->readString(Form("Random:seed = %llu", random_value % 900000001));

  // print debug
  // gen->PrintDebug();

  return gen;
}

FairGenerator*
  GeneratorPromptJpsiPsi2S_EvtGenMidY(int triggerGap, double rapidityMin = -1.5, double rapidityMax = 1.5, bool verbose = false, bool embedding = false)
{
  auto gen = new o2::eventgen::GeneratorEvtGen<o2::eventgen::GeneratorPythia8PromptSignalsEmbeddingGapTriggered>();
  gen->setTriggerGap(triggerGap);
  gen->setRapidityRange(rapidityMin, rapidityMax);
  gen->addHadronPDGs(443);
  gen->addHadronPDGs(100443);
  gen->setVerbose(verbose);
  TString pathO2table = gSystem->ExpandPathName("${O2DPG_MC_CONFIG_ROOT}/MC/config/PWGDQ/pythia8/decayer/switchOffJpsi.cfg");
  gen->readFile(pathO2table.Data());
  gen->setConfigMBdecays(pathO2table);
  gen->PrintDebug(true);
  gen->SetSizePdg(2);
  gen->AddPdg(443, 0);
  gen->AddPdg(100443, 1);
  gen->SetForceDecay(kEvtDiElectron);
  // set random seed
  gen->readString("Random:setSeed on");
  uint random_seed;
  unsigned long long int random_value = 0;
  ifstream urandom("/dev/urandom", ios::in|ios::binary);
  urandom.read(reinterpret_cast<char*>(&random_value), sizeof(random_seed));
  gen->readString(Form("Random:seed = %llu", random_value % 900000001));
  // print debug
  // gen->PrintDebug();
  return gen;
}

FairGenerator*
  GeneratorPromptJpsi_EvtGenFwdy(int triggerGap, double rapidityMin = -4.3, double rapidityMax = -2.3, bool verbose = false, bool embedding = false)
{
  auto gen = new o2::eventgen::GeneratorEvtGen<o2::eventgen::GeneratorPythia8PromptSignalsEmbeddingGapTriggered>();
  gen->setTriggerGap(triggerGap);
  gen->setRapidityRange(rapidityMin, rapidityMax);
  gen->addHadronPDGs(443);
  gen->setVerbose(verbose);

  TString pathO2table = gSystem->ExpandPathName("${O2DPG_MC_CONFIG_ROOT}/MC/config/PWGDQ/pythia8/decayer/switchOffJpsi.cfg");
  gen->readFile(pathO2table.Data());
  gen->setConfigMBdecays(pathO2table);
  gen->PrintDebug(true);

  gen->SetSizePdg(1);
  gen->AddPdg(443, 0);

  gen->SetForceDecay(kEvtDiMuon);

  // set random seed
  gen->readString("Random:setSeed on");
  uint random_seed;
  unsigned long long int random_value = 0;
  ifstream urandom("/dev/urandom", ios::in|ios::binary);
  urandom.read(reinterpret_cast<char*>(&random_value), sizeof(random_seed));
  gen->readString(Form("Random:seed = %llu", random_value % 900000001));

  // print debug
  // gen->PrintDebug();

  return gen;
}

FairGenerator*
  GeneratorPromptJpsiPsi2S_EvtGenFwdY(int triggerGap, double rapidityMin = -4.3, double rapidityMax = -2.3, bool verbose = false, bool embedding = false)
{
  auto gen = new o2::eventgen::GeneratorEvtGen<o2::eventgen::GeneratorPythia8PromptSignalsEmbeddingGapTriggered>();
  gen->setTriggerGap(triggerGap);
  gen->setRapidityRange(rapidityMin, rapidityMax);
  gen->addHadronPDGs(443);
  gen->addHadronPDGs(100443);
  gen->setVerbose(verbose);
  TString pathO2table = gSystem->ExpandPathName("${O2DPG_MC_CONFIG_ROOT}/MC/config/PWGDQ/pythia8/decayer/switchOffJpsi.cfg");
  gen->readFile(pathO2table.Data());
  gen->setConfigMBdecays(pathO2table);
  gen->PrintDebug(true);
  gen->SetSizePdg(2);
  gen->AddPdg(443, 0);
  gen->AddPdg(100443, 1);
  gen->SetForceDecay(kEvtDiMuon);
  // set random seed
  gen->readString("Random:setSeed on");
  uint random_seed;
  unsigned long long int random_value = 0;
  ifstream urandom("/dev/urandom", ios::in|ios::binary);
  urandom.read(reinterpret_cast<char*>(&random_value), sizeof(random_seed));
  gen->readString(Form("Random:seed = %llu", random_value % 900000001));
  // print debug
  // gen->PrintDebug();
  return gen;
}
