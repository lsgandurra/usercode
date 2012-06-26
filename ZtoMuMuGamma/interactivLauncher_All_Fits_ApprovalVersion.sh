#./bin/bash
date="Fall11"

for FitPourcentage in `seq 80 100`
do
	for isMC in '1' '0'
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
			for r9sup in '2' ##`seq 0 2`
			do
		
				for SetOfCorrections in 'Regression' ##'ElectronTunedCorrections'
		                do
					for phiCracks in '1' ##'0' '1'
		                        do
		                                for etaCracks in '1' ##'0' '1'
		                                do
		                                        for variableX in 'Photon_Et' ##'Photon_SC_rawEt' 'Photon_Et' 'Photon_E' 'Photon_SC_Eta' 'Photon_SC_brem' 
		                                        do
		                                                for MuonCorrection in 'Rochester' ##'NoMuonCorrection' 'Rochester'
		                                                do
		                                                        for Category in 'OneBin' ##'Vgamma8' 'OneBin' 'Vgamma24'
		                                                        do
										for HighMmumuLim in '80'
		                                                                do
		
											./mmg_s_Fits_ApprovalVersion.exe $Endcaps $r9sup $SetOfCorrections RooVoigtian2 2011 40 $HighMmumuLim ${FitPourcentage} $phiCracks $etaCracks $variableX ${isMC} 5GeV ${MuonCorrection} Profile_Surface ${Category} ${date}
										done
									done
								done
							done
						done
					done
				done
			done
		done
	done
done


for FitPourcentage in `seq 80 100`
do
	for isMC in '1' '0'
	do
	
		if [ "$isMC" = "0" ]
	        then
	                date="16Jan"
	        fi
	        if [ "$isMC" = "1" ]
	        then
	                date="Fall11"
	        fi
	
		for Endcaps in '0' ##`seq 0 1`
		do
			for r9sup in '1' ##`seq 0 2`
			do
		
				for SetOfCorrections in 'Regression' ##'ElectronTunedCorrections'
		                do
					for phiCracks in '1' ##'0' '1'
		                        do
		                                for etaCracks in '1' ##'0' '1'
		                                do
		                                        for variableX in 'Photon_Et' ##'Photon_SC_rawEt' 'Photon_Et' 'Photon_E' 'Photon_SC_Eta' 'Photon_SC_brem' 
		                                        do
		                                                for MuonCorrection in 'Rochester' ##'NoMuonCorrection' 'Rochester'
		                                                do
		                                                        for Category in 'OneBin' ##'Vgamma8' 'OneBin' 'Vgamma24'
		                                                        do
										for HighMmumuLim in '80'
		                                                                do
		
											./mmg_s_Fits_ApprovalVersion.exe $Endcaps $r9sup $SetOfCorrections RooVoigtian2 2011 40 $HighMmumuLim ${FitPourcentage} $phiCracks $etaCracks $variableX ${isMC} 5GeV ${MuonCorrection} Profile_Surface ${Category} ${date}
										done
									done
								done
							done
						done
					done
				done
			done
		done
	done
done


###for FitPourcentage in `seq 60 79`
###do
###	for isMC in '1' 
###	do
###	
###		if [ "$isMC" = "0" ]
###	        then
###	                date="16Jan"
###	        fi
###	        if [ "$isMC" = "1" ]
###	        then
###	                date="Fall11"
###	        fi
###	
###		for Endcaps in '0' '1' ##`seq 0 1`
###		do
###			for r9sup in '2' ##`seq 0 2`
###			do
###		
###				for SetOfCorrections in 'Regression' ##'ElectronTunedCorrections'
###		                do
###					for phiCracks in '1' ##'0' '1'
###		                        do
###		                                for etaCracks in '1' ##'0' '1'
###		                                do
###		                                        for variableX in 'Photon_Et' ##'Photon_SC_rawEt' 'Photon_Et' 'Photon_E' 'Photon_SC_Eta' 'Photon_SC_brem' 
###		                                        do
###		                                                for MuonCorrection in 'Rochester' ##'NoMuonCorrection' 'Rochester'
###		                                                do
###		                                                        for Category in 'OneBin' ##'Vgamma8' 'OneBin' 'Vgamma24'
###		                                                        do
###										for HighMmumuLim in '80'
###		                                                                do
###		
###											./Photon_E_o_MC_Fits_ApprovalVersion.exe $Endcaps $r9sup $SetOfCorrections RooVoigtian2 2011 40 $HighMmumuLim ${FitPourcentage} $phiCracks $etaCracks $variableX ${isMC} 5GeV ${MuonCorrection} Profile_Surface ${Category} ${date}
###											./Photon_E_o_MC_Fits_ApprovalVersion.exe $Endcaps $r9sup $SetOfCorrections RooBreitWigner2 2011 40 $HighMmumuLim ${FitPourcentage} $phiCracks $etaCracks $variableX ${isMC} 5GeV ${MuonCorrection} Profile_Surface ${Category} ${date}
###                                                                                        ./Photon_E_o_MC_Fits_ApprovalVersion.exe $Endcaps $r9sup $SetOfCorrections RooBifurcatedGauss 2011 40 $HighMmumuLim ${FitPourcentage} $phiCracks $etaCracks $variableX ${isMC} 5GeV ${MuonCorrection} Profile_Surface ${Category} ${date}
###   ###                                                                                     ./Photon_E_o_MC_Fits_ApprovalVersion.exe $Endcaps $r9sup $SetOfCorrections RooCrystalBall 2011 40 $HighMmumuLim ${FitPourcentage} $phiCracks $etaCracks $variableX ${isMC} 5GeV ${MuonCorrection} Profile_Surface ${Category} ${date}
###										done
###									done
###								done
###							done
###						done
###					done
###				done
###			done
###		done
###	done
###done
###
###
###for FitPourcentage in `seq 60 79`
###do
###	for isMC in '1' 
###	do
###	
###		if [ "$isMC" = "0" ]
###	        then
###	                date="16Jan"
###	        fi
###	        if [ "$isMC" = "1" ]
###	        then
###	                date="Fall11"
###	        fi
###	
###		for Endcaps in '0' ##`seq 0 1`
###		do
###			for r9sup in '1' ##`seq 0 2`
###			do
###		
###				for SetOfCorrections in 'Regression' ##'ElectronTunedCorrections'
###		                do
###					for phiCracks in '1' ##'0' '1'
###		                        do
###		                                for etaCracks in '1' ##'0' '1'
###		                                do
###		                                        for variableX in 'Photon_Et' ##'Photon_SC_rawEt' 'Photon_Et' 'Photon_E' 'Photon_SC_Eta' 'Photon_SC_brem' 
###		                                        do
###		                                                for MuonCorrection in 'Rochester' ##'NoMuonCorrection' 'Rochester'
###		                                                do
###		                                                        for Category in 'OneBin' ##'Vgamma8' 'OneBin' 'Vgamma24'
###		                                                        do
###										for HighMmumuLim in '80'
###		                                                                do
###		
###											./Photon_E_o_MC_Fits_ApprovalVersion.exe $Endcaps $r9sup $SetOfCorrections RooVoigtian2 2011 40 $HighMmumuLim ${FitPourcentage} $phiCracks $etaCracks $variableX ${isMC} 5GeV ${MuonCorrection} Profile_Surface ${Category} ${date}
###											./Photon_E_o_MC_Fits_ApprovalVersion.exe $Endcaps $r9sup $SetOfCorrections RooBreitWigner2 2011 40 $HighMmumuLim ${FitPourcentage} $phiCracks $etaCracks $variableX ${isMC} 5GeV ${MuonCorrection} Profile_Surface ${Category} ${date}
###											./Photon_E_o_MC_Fits_ApprovalVersion.exe $Endcaps $r9sup $SetOfCorrections RooBifurcatedGauss 2011 40 $HighMmumuLim ${FitPourcentage} $phiCracks $etaCracks $variableX ${isMC} 5GeV ${MuonCorrection} Profile_Surface ${Category} ${date}
######											./Photon_E_o_MC_Fits_ApprovalVersion.exe $Endcaps $r9sup $SetOfCorrections RooCrystalBall 2011 40 $HighMmumuLim ${FitPourcentage} $phiCracks $etaCracks $variableX ${isMC} 5GeV ${MuonCorrection} Profile_Surface ${Category} ${date}	
###										done
###									done
###								done
###							done
###						done
###					done
###				done
###			done
###		done
###	done
###done
###
##RooBreitWigner2 RooBifurcatedGauss RooCrystalBall
