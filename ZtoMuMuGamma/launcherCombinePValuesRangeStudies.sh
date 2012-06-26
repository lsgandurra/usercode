#./bin/bash
date="Fall11"

for isMC in '0' '1' ##`seq 0 1`
do
     if [ "$isMC" = "0" ]
     then
     	date="16Jan"
     fi
     if [ "$isMC" = "1" ]
     then
	date="Fall11"
     fi


	for Endcaps in '0' '1' ##`seq 0 1`
	do
	        for r9sup in '0' '1' '2' ##`seq 0 2`
	        do
	                for function in 'RooVoigtian2' ##'RooLandau2' 'RooLandauConvGaussian' 'RooCrystalBall' 'RooBifurcatedGauss'
	                do
	                        for variableX in 'Photon_Et' ##'Photon_SC_rawEt' 'Photon_Et' 'Photon_E' 'Photon_SC_Eta' 'Photon_SC_brem'
	                        do
	                                for MuonCorrection in 'Rochester' ##'NoMuonCorrection' 'Rochester'
	                                do
	                                	for Category in 'OneBin' ##'Vgamma24' 'Vgamma8' 'OneBin'
	                                        do                
							for SurfaceMethod in 'Profile_Surface'
	                                                do
	                                                                for MZbinning in '5GeV' ##'05GeV' '1GeV' '2GeV' '5GeV'
	                                                                do
	                                                                         ##root -l -b -q "CombinePValuesRangeStudiesF_Juin_2012.C($Endcaps,$r9sup,\"Photon_E_o_MC_E\",\"Regression\",\"${function}\",\"2011\",\"40\",\"80\",1,1,\"$variableX\",0,50,60,\"${MZbinning}\",\"${MuonCorrection}\",\"${SurfaceMethod}\",\"${Category}\",\"${date}\")"   
	                                                                         ##root -l -b -q "CombinePValuesRangeStudiesF_Juin_2012.C($Endcaps,$r9sup,\"mmg_s\",\"Regression\",\"${function}\",\"2011\",\"40\",\"80\",1,1,\"$variableX\",${isMC},60,70,\"${MZbinning}\",\"${MuonCorrection}\",\"${SurfaceMethod}\",\"${Category}\",\"${date}\")"
	                                                                         ##root -l -b -q "CombinePValuesRangeStudiesF_Juin_2012.C($Endcaps,$r9sup,\"mmg_s\",\"Regression\",\"${function}\",\"2011\",\"40\",\"80\",1,1,\"$variableX\",${isMC},70,80,\"${MZbinning}\",\"${MuonCorrection}\",\"${SurfaceMethod}\",\"${Category}\",\"${date}\")"
	                                                                         root -l -b -q "CombinePValuesRangeStudiesF_Juin_2012.C($Endcaps,$r9sup,\"mmg_s\",\"Regression\",\"${function}\",\"2011\",\"40\",\"80\",1,1,\"$variableX\",${isMC},80,90,\"${MZbinning}\",\"${MuonCorrection}\",\"${SurfaceMethod}\",\"${Category}\",\"${date}\")"
	                                                                         root -l -b -q "CombinePValuesRangeStudiesF_Juin_2012.C($Endcaps,$r9sup,\"mmg_s\",\"Regression\",\"${function}\",\"2011\",\"40\",\"80\",1,1,\"$variableX\",${isMC},90,101,\"${MZbinning}\",\"${MuonCorrection}\",\"${SurfaceMethod}\",\"${Category}\",\"${date}\")"
									done
							done	
						done
					done
				done
			done
		done
	done
done

###for isMC in '1' ##`seq 0 1`
###do
###     if [ "$isMC" = "0" ]
###     then
###     	date="16Jan"
###     fi
###     if [ "$isMC" = "1" ]
###     then
###	date="Fall11"
###     fi
###
###
###	for Endcaps in '0' '1' ##`seq 0 1`
###	do
###	        for r9sup in '0' '1' '2' ##`seq 0 2`
###	        do
###	                for function in 'RooVoigtian2' ##'RooLandau2' 'RooLandauConvGaussian' 'RooCrystalBall' 'RooBifurcatedGauss'
###	                do
###	                        for variableX in 'Photon_Et' ##'Photon_SC_rawEt' 'Photon_Et' 'Photon_E' 'Photon_SC_Eta' 'Photon_SC_brem'
###	                        do
###	                                for MuonCorrection in 'Rochester' ##'NoMuonCorrection' 'Rochester'
###	                                do
###	                                	for Category in 'OneBin' ##'Vgamma24' 'Vgamma8' 'OneBin'
###	                                        do                
###							for SurfaceMethod in 'Profile_Surface'
###	                                                do
###	                                                                for MZbinning in '5GeV' ##'05GeV' '1GeV' '2GeV' '5GeV'
###	                                                                do
###	                                                                         ##root -l -b -q "CombinePValuesRangeStudiesF_Juin_2012.C($Endcaps,$r9sup,\"Photon_E_o_MC_E\",\"Regression\",\"${function}\",\"2011\",\"40\",\"80\",1,1,\"$variableX\",0,50,60,\"${MZbinning}\",\"${MuonCorrection}\",\"${SurfaceMethod}\",\"${Category}\",\"${date}\")"   
###	                                                                         root -l -b -q "CombinePValuesRangeStudiesF_Juin_2012.C($Endcaps,$r9sup,\"Photon_E_o_MC_E\",\"Regression\",\"${function}\",\"2011\",\"40\",\"80\",1,1,\"$variableX\",${isMC},60,70,\"${MZbinning}\",\"${MuonCorrection}\",\"${SurfaceMethod}\",\"${Category}\",\"${date}\")"
###	                                                                         root -l -b -q "CombinePValuesRangeStudiesF_Juin_2012.C($Endcaps,$r9sup,\"Photon_E_o_MC_E\",\"Regression\",\"${function}\",\"2011\",\"40\",\"80\",1,1,\"$variableX\",${isMC},70,80,\"${MZbinning}\",\"${MuonCorrection}\",\"${SurfaceMethod}\",\"${Category}\",\"${date}\")"
###	                                                                         root -l -b -q "CombinePValuesRangeStudiesF_Juin_2012.C($Endcaps,$r9sup,\"Photon_E_o_MC_E\",\"Regression\",\"${function}\",\"2011\",\"40\",\"80\",1,1,\"$variableX\",${isMC},80,90,\"${MZbinning}\",\"${MuonCorrection}\",\"${SurfaceMethod}\",\"${Category}\",\"${date}\")"
###	                                                                         root -l -b -q "CombinePValuesRangeStudiesF_Juin_2012.C($Endcaps,$r9sup,\"Photon_E_o_MC_E\",\"Regression\",\"${function}\",\"2011\",\"40\",\"80\",1,1,\"$variableX\",${isMC},90,101,\"${MZbinning}\",\"${MuonCorrection}\",\"${SurfaceMethod}\",\"${Category}\",\"${date}\")"
###									done
###							done	
###						done
###					done
###				done
###			done
###		done
###	done
###done
###




###int CombinePValuesRangeStudiesF_Avril_2012(int EndCaps = 0, int r9sup = 1, string scale = "mmg_s_MZ_Surface", string SetOfCorrections = "Regression", string nomFitMethode = "RooCrystalBall", string PileUpVersion = "2011", string LowMmumuLim = "40", string HightMmumuLim = "80", bool phiCracks = true, bool etaCracks = true, string variableX = "Photon_Et", int isMC = 1, int LowPourcent = 80, int HighPourcent = 90, string MZbinning = "05GeV", string MuonCorrection = "NoMuonCorrection", string SurfaceMethod = "Fitted_Surface", string Category = "Vgamma8", string date = "16Jan")
