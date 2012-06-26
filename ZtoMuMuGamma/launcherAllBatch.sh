#./bin/bash
for isMC in '1' ##'0' '1'
do
        for FitPourcentage in `seq 60 100`
        ##for FitPourcentage in `seq 50 60`	
	do
                for Function in 'RooBifurcatedGauss' ##'RooLandau2' 'RooLandauConvGaussian' 'RooCrystalBall' 'RooBifurcatedGauss'
                do
				##for MuonCorrection in 'NoMuonCorrection' 'Rochester'
	
			for SurfaceMethod in 'Profile_Surface' ##'Fitted_Surface' 'Profile_Surface'
			do
				if [ "$SurfaceMethod" = "Fitted_Surface" ]
				then
					for MZbinning in '2GeV' '5GeV'
                        		do
						if [ "$isMC" = "0" ]
                                		then	
							qsub batchJob.sh $isMC $FitPourcentage $Function $MZbinning $SurfaceMethod 16Jan
                					echo "isMC = $isMC, FitPourcentage = $FitPourcentage, Function= $Function, MZbinning= $MZbinning, SurfaceMethod= $SurfaceMethod date= 16Jan"
						fi
						if [ "$isMC" = "1" ]
                                                then
                                                        qsub batchJob.sh $isMC $FitPourcentage $Function $MZbinning $SurfaceMethod Fall11
                                                        echo "isMC = $isMC, FitPourcentage = $FitPourcentage, Function= $Function, MZbinning= $MZbinning, SurfaceMethod= $SurfaceMethod date= Fall11"
                                                fi
					done	
				fi
				if [ "$SurfaceMethod" = "Profile_Surface" ]
                                then
                                        for MZbinning in '05GeV' '1GeV' '2GeV' '5GeV'
                                        do
						if [ "$isMC" = "0" ]
                                                then
                                                        qsub batchJob.sh $isMC $FitPourcentage $Function $MZbinning $SurfaceMethod 16Jan
                                                        echo "isMC = $isMC, FitPourcentage = $FitPourcentage, Function= $Function, MZbinning= $MZbinning, SurfaceMethod= $SurfaceMethod date= 16Jan"
                                                fi
                                                if [ "$isMC" = "1" ]
                                                then
                                                        qsub batchJob.sh $isMC $FitPourcentage $Function $MZbinning $SurfaceMethod Fall11
                                                        echo "isMC = $isMC, FitPourcentage = $FitPourcentage, Function= $Function, MZbinning= $MZbinning, SurfaceMethod= $SurfaceMethod date= Fall11"
                                                fi
                                        done
                                fi
			done
		done
        done
done


###qsub batchJob.sh 0 88 RooBifurcatedGauss Profile_Surface 05GeV 16Jan
