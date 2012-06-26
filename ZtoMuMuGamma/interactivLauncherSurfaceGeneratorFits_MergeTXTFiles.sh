#! /usr/local/bin/bash -l


for typeOf in 'Data'
do
	for MuonCorrections in 'NoMuonCor' 'RochCor'
	do
		for date in '16Jan2012' '30Nov2011' ##'16Jan2012' '30Nov2011' 'Vgamma'
		do 
			for choix in 'Mmumu.txt' ##'plot3D' 'plot2D' 'Mmumu.txt'
			do
				for SurfaceBinning in '2GeV' '5GeV' ##'05GeV' '1GeV' '2GeV' '5GeV'
				do 
					./SurfaceGeneratorFit_v2_MergeTXTFiles.exe ${typeOf} ${MuonCorrections} ${date} ${choix} ${SurfaceBinning} 
				done
			done
		done
	done
done


###for typeOf in 'MC'
###do
###        for MuonCorrections in 'NoMuonCor' 'RochCor'
###        do  
###                for date in 'Fall11'
###                do
###                	for choix in 'Mmumu.txt' ##'plot3D' 'plot2D' 'Mmumu.txt'
###                        do  
###                                for SurfaceBinning in '2GeV' '5GeV' ##'05GeV' '1GeV' '2GeV' '5GeV'
###                                do  
###                                        ./SurfaceGeneratorFit_v2_MergeTXTFiles.exe ${typeOf} ${MuonCorrections} ${date} ${choix} ${SurfaceBinning} 
###                                done
###                        done	
###
###		done
###        done
###done
###


###int SurfaceGenerator(string type = "Data", string RochesterCorrections = "NoMuonCor", string date = "16Jan2012", string choix = "Mmumu.txt", string SurfaceBinning = "5GeV") 

