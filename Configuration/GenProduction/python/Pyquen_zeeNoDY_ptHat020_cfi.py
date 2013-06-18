import FWCore.ParameterSet.Config as cms

from Configuration.Generator.PyquenDefaultSettings_cff import *

hiSignal = cms.EDFilter("PyquenGeneratorFilter",
                        collisionParameters2760GeV,
                        qgpParameters,
                        pyquenParameters,
                        doQuench = cms.bool(False),
                        bFixed = cms.double(0.0), ## fixed impact param (fm); valid only if cflag_=0
                        cFlag = cms.int32(0), ## centrality flag
                        bMin = cms.double(0.0), ## min impact param (fm); valid only if cflag_!=0
                        bMax = cms.double(0.0), ## max impact param (fm); valid only if cflag_!=0
                        PythiaParameters = cms.PSet(pyquenPythiaDefaultBlock,
                                                    parameterSets = cms.vstring('pythiaUESettings',
                                                                                'pythiaZmumu',
                                                                                'kinematics'
                                                                                ),
                                                    pythiaZmumu = cms.vstring('MSEL = 0 ! users defined processes only',
                                                                              'MSUB(15)=1          !qq->Z0/gamma*+g',
	                                                                      'MSUB(30)=1          !qg->Z0/gamma*+q',                                                                                             'MSTP(43)=2          !Only Z0',
                                                                              'MDME( 174,1) = 0    !Z decay into d dbar', 
                                                                              'MDME( 175,1) = 0    !Z decay into u ubar', 
                                                                              'MDME( 176,1) = 0    !Z decay into s sbar', 
                                                                              'MDME( 177,1) = 0    !Z decay into c cbar', 
                                                                              'MDME( 178,1) = 0    !Z decay into b bbar', 
                                                                              'MDME( 179,1) = 0    !Z decay into t tbar', 
                                                                              'MDME( 182,1) = 1    !Z decay into e- e+', 
                                                                              'MDME( 183,1) = 0    !Z decay into nu_e nu_ebar', 
                                                                              'MDME( 184,1) = 0    !Z decay into mu- mu+', 
                                                                              'MDME( 185,1) = 0    !Z decay into nu_mu nu_mubar', 
                                                                              'MDME( 186,1) = 0    !Z decay into tau- tau+', 
                                                                              'MDME( 187,1) = 0    !Z decay into nu_tau nu_taubar', 
                                                                              ),
                                                    kinematics = cms.vstring ('CKIN(3) = 0.       !(D=0 GeV) lower lim pT_hat',
                                                                              'CKIN(4) = 20.       !(D=-1 GeV) upper lim pT_hat, if < 0 innactive',
                                                                              "CKIN(7)=-2.5",  #min rapidity
                                                                              "CKIN(8)=2.5"    #max rapidity
                                                                              )
                                                    )
                       
                        )

hiSignal.embeddingMode = True

ProductionFilterSequence = cms.Sequence(hiSignal)
