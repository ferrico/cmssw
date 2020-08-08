import FWCore.ParameterSet.Config as cms

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from RecoMuon.TrackingTools.MuonTrackLoader_cff import *
from RecoMuon.GlobalTrackingTools.GlobalTrajectoryBuilderCommon_cff import *
globalMuonsFLORIDA = cms.EDProducer("GlobalMuonProducerFLORIDA",
    MuonServiceProxy,
    MuonTrackLoaderForGLBFLORIDA,
    GLBTrajBuilderParameters = cms.PSet(
        GlobalTrajectoryBuilderCommon
    ),
    TrackerCollectionLabel = cms.InputTag("generalTracks"),
    MuonCollectionLabel = cms.InputTag("standAloneMuons","UpdatedAtVtx")
)

globalMuonsFLORIDA.GLBTrajBuilderParameters.GlobalMuonTrackMatcher.Propagator = cms.string('SmartPropagatorRK')
globalMuonsFLORIDA.GLBTrajBuilderParameters.TrackTransformer.Propagator = cms.string('SmartPropagatorAnyRK') 

# FastSim has no template fit on tracker hits
# FastSim doesn't use Runge Kute for propagation
from Configuration.Eras.Modifier_fastSim_cff import fastSim
fastSim.toModify(globalMuonsFLORIDA,
                 GLBTrajBuilderParameters =dict(GlbRefitterParameters = dict(TrackerRecHitBuilder = 'WithoutRefit',
                                                                             Propagator = 'SmartPropagatorAny'),
                                                TrackerRecHitBuilder = 'WithoutRefit',
                                                TrackTransformer = dict(TrackerRecHitBuilder = 'WithoutRefit',
                                                                        Propagator = 'SmartPropagatorAny'),
                                                GlobalMuonTrackMatcher = dict(Propagator = 'SmartPropagator') 
                                                )
)

