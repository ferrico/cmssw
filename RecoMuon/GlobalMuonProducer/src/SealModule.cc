#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoMuon/GlobalMuonProducer/src/GlobalMuonProducer.h"
#include "RecoMuon/GlobalMuonProducer/src/GlobalMuonProducerFLORIDA.h"
#include "RecoMuon/GlobalMuonProducer/src/TevMuonProducer.h"


DEFINE_FWK_MODULE(GlobalMuonProducer);
DEFINE_FWK_MODULE(GlobalMuonProducerFLORIDA);
DEFINE_FWK_MODULE(TevMuonProducer);

