
**UL 2016, 2017 & 2018 High Pt Photon ID for EB**

1. First, XGBoost has to be built either from source (https://xgboost.readthedocs.io/en/latest/build.html - version release_1.2.0 or higher) or using pip

		 pip3/pip --user install xgboost

2.  Indentify the installation directory of XGBoost. 

3. Include the xgboost library and helper cc in your ROOT macro:

		 R__ADD_INCLUDE_PATH(path_to_xgboost/include/xgboost/) 
		 R__LOAD_LIBRARY(path_to_xgboost/lib/libxgboost.so) 
		 #include <path_to_xgboost/include/xgboost/c_api.h>
		 #include "helpers.cc"

4. Load the BDT model and isolation corrections:
		
		BDT needs to be loaded only once (no need to repeat for each event)

			2017/2018 BDT model file: aNTGC_photon_BDT_EB_2021_08_26_09_39_52.model

			2016 BDT model file: aNTGC_photon_BDT_EB_2022_01_10_23_54_43.model
			
				    DMatrixHandle 		dTest;
				    BoosterHandle 		phoBDT_h;
				    XGBoosterCreate(NULL, 0, &phoBDT_h); 
				    XGBoosterSetParam(phoBDT_h, "seed", "0"); 
				    Int_t mLdSuccess = XGBoosterLoadModel(phoBDT_h, "BDT_model_path");
				    if(mLdSuccess !=0) std::cout<<"Failed to load model"<<std::endl;


	    Isolation corrections

		    For all years:
				isoCorrMap ecalIsoRhoCorrMap("PATH/phoPFClusEcalIso_PtCorrections.txt", 2);
				isoCorrMap ecalIsoPtCorrMap("PATH/phoPFClusEcalIso_PtCorrections.txt", 2);

			For 2017/2018:
				isoCorrMap tkrIsoRhoCorrMap("PATH/phoTrkSumPtHollowConeDR03_RhoCorrections.txt", 2);

			For 2016:
				isoCorrMap worstChHadIsoRhoCorrMap("PATH/phoPFChWorstIso_RhoCorrections.txt", 2);

5. Predict the BDT score per photon.

		std::vector<Float_t> feats{_phoE2x2Full5x5[iPho] / (_phoR9Full5x5[iPho] * _ecalSC_RawEn[phoSCindex]),
		                               std::abs(_phoSCeta),
		                               _phoE1x3Full5x5[iPho] / _ecalSC_RawEn[phoSCindex],
		                               _phoE2ndFull5x5[iPho] / _ecalSC_RawEn[phoSCindex],
		                               _phoE2x5Full5x5[iPho] / _ecalSC_RawEn[phoSCindex],
		                               _phoMaxEnergyXtal[iPho] / _ecalSC_RawEn[phoSCindex],
		                               _ecalSC_etaWidth[phoSCindex] / _ecalSC_phiWidth[phoSCindex],
		                               _ecalSC_etaWidth[phoSCindex],
		                               _ecalSC_phiWidth[phoSCindex],
		                               _phoCalibEt[phoSCindex],
		                               _phoR9Full5x5[iPho],
		                               _phoE2x2Full5x5[iPho] / _ecalSC_RawEn[phoSCindex],
		                               _phoSigmaIEtaIEtaFull5x5[iPho] / _phoSigmaIPhiIPhiFull5x5[iPho],
		                               _phoSigmaIEtaIEtaFull5x5[iPho],
		                               _phoSigmaIEtaIPhiFull5x5[iPho],
		                               _phoSigmaIPhiIPhiFull5x5[iPho]};


		XGDMatrixCreateFromMat((float*)feats.data(), 1, feats.size(), -9999999999, &dTest);
		bst_ulong out_len;
		const float *prediction;
		XGBoosterPredict(phoBDT_h, dTest, 0, 0, 0, &out_len, &prediction);
		assert(out_len == 1);
		XGDMatrixFree(dTest);
		Float_t phoBDTpred_ = prediction[0]; /// This is the prediciton
	
	The shower shape variables can be extracted from these two methods of PAT::photon

				PAT::Photon::full5x5_showerShapeVariables()
				PAT::Photon::superCluster()
	
	Refer to CMSSW documentation:				[full5x5_showerShapeVariables](https://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_10_6_24/doc/html/d0/d08/structreco_1_1Photon_1_1ShowerShape.html) and	[superCluster](https://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_10_6_24/doc/html/d2/de8/classreco_1_1SuperCluster.html)

6. Apply pileup and pT corrections to these isolations:
				
		ECAL isolation:
			phoPFClusEcalIso_ = PAT::Photon::ecalPFClusterIso()

		Two different tracker isolations:
			
			For 2017/2018
				phoTrkSumPtHollowConeDR03_ = PAT::Photon::trkSumPtHollowConeDR03()

			For 2016 
				worstChHadIso_ 			   = PAT::Photon::chargedHadronWorstVtxIso() 


		Additional inputs for correction:
				
				rho_ 						: "fixedGridRhoFastjetAll"

				phoPt_     					= PAT::Photon::et() * PAT::Photon::userFloat("ecalEnergyPostCorr") / PAT::Photon::energy());

				phoAbsSCeta_				= PAT::Photon::superCluster()->eta();


		Apply corrections:

				For all years:
					phoPFECALClusIsoCorr_      	= phoPFClusEcalIso_ - ecalIsoRhoCorrMap.getIsoCorr(phoAbsSCeta_, rho_, 0) - ecalIsoPtCorrMap.getIsoCorr(phoAbsSCeta_, phoPt_, 1);
				
				For 2017/2018:
					phoTkrIsoCorr_             	= phoTrkSumPtHollowConeDR03_ - tkrIsoRhoCorrMap.getIsoCorr(phoAbsSCeta_, rho_, 1);
				
				For 2016 only:
					worstChHadIsoCorr_          = worstChHadIso_ - worstChHadIsoRhoCorrMap.getIsoCorr(phoAbsSCeta_, rho_, 1);
	

7. Get ID decision:

				**2016**: passBDTid_ = (phoBDTpred_ > 0.8493 && phoHoverE_ <  0.0222 &&  phoPFECALClusIsoCorr92_ <  2.16 && worstChHadIsoCorr_ < 2.19);

				**2017**: passBDTid_ = (phoBDTpred_ > 0.8361 && phoHoverE_ <  0.04012 &&  phoPFECALClusIsoCorr_ <  1.84 && phoTkrIsoCorr_ < 1.63);

				**2018**: passBDTid_ = (phoBDTpred_ > 0.8466 && phoHoverE_ <  0.0402 &&  phoPFECALClusIsoCorr_ <  1.84 && phoTkrIsoCorr_ < 1.58);
