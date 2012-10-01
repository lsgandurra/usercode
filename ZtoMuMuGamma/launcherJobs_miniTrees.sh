#./bin/bash

for sample2 in 'TTJets_TuneZ2_7TeV-madgraph-tauola' 'WJetsToLNu_TuneZ2_7TeV-madgraph-tauola' ###'Run2011A-ZMu-May10ReReco-v1' 'Run2011A-ZMu-PromptSkim-v4' 'Run2011A-ZMu-05Aug2011-v1_V04' 'Run2011A-ZMu-03Oct2011-v1' 'Run2011B-ZMu-PromptSkim-v1_finalJson' 'Run2011A-ZMu-05Jul2011ReReco-ECAL-v1' 'DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Fall11-PU_S6_START42_V14B-v1' 'DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_v2' 'Run2011A-ZMu-05Jul2011ReReco-ECAL-v1'
do
	for Zgamma2 in '3' ##'0' '1' '2' '3'
	do
		for SetOfCorrections2 in 'MITregression' ##'1.0' 'ETHZ' 'START42_V11'
		do
			for muonsCor in '3' ##'0' 
			do
			##for HighMmumuLim in '76' '84' ##36
			##do 
			##done		
			##for HighMmumuLim in '78' '82' '84' '86' '88' ##38
			##do
			##done
			##for HighMmumuLim in '80' ##40
			##do
			##done
			##for HighMmumuLim in '82' ##42
			##do
			##done
			##for HighMmumuLim in '84' ##44
			##do
			##done
				qsub batchJob_miniTrees.sh ${sample2} ${sample2}_40_80_${Zgamma2}_v11 $Zgamma2 40 80 $SetOfCorrections2 1.0 $muonsCor
				##hadd miniTree_${sample2}_40_80_${Zgamma2}_v11_partALL.root miniTree_${sample2}_40_80_${Zgamma2}_v11_part[0-9]*root
                        	##echo "miniTree_${sample2}_${Zgamma2}_${SetOfCorrections2}_40_80_v11_part[0-9].root"
				##mv miniTree_${sample2}_40_80_${Zgamma2}_v11_part[0-9]*root stored_miniTree/
                        	##mv ${sample2}_40_80_${Zgamma2}_v11_part[0-9]*err stored_errput/
                        	##mv ${sample2}_40_80_${Zgamma2}_v11_part[0-9]*out stored_output/
			done
		done


	done


done


###./Selection_miniTree.exe ${sample2} ${sample2}_${SetOfCorrections2}_Rochester_40_80_5GeV_partALL 9999 -1 0 2011 PU_S6 40 80 ${SetOfCorrections2} 1.0 0

###./Selection_miniTree.exe Run2011A-ZMu-May10ReReco-v1 Run2011A-ZMu-May10ReReco-v1 9999 0 0 2011 PU_S6 40 80 ETHZ 1.0 3
