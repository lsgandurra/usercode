#./bin/bash

###for MuonCorrection in 'NoMuonCor' 'RochCor'
###do
###	for ijob in `seq 0 44`
###	do
###		qsub batchJobSurfaceGeneratorFit_v2b.sh $MuonCorrection $ijob MC Fall11 2GeV 45
###        	echo "MuonCorrection= $MuonCorrection, ijob= $ijob, MC, 2GeV"
###	done
###done
###
for date in '16Jan2012' '30Nov2011'
do

	for MuonCorrection in 'NoMuonCor' 'RochCor'
	do
	        for ijob in `seq 0 44`
	        do
	                qsub batchJobSurfaceGeneratorFits.sh $MuonCorrection $ijob Data $date 2GeV 45
	                echo "MuonCorrection= $MuonCorrection, ijob= $ijob, Data, 2GeV"
	        done
	done
	
	for MuonCorrection in 'NoMuonCor' 'RochCor'
	do
	        for ijob in `seq 0 17`
	        do
	                qsub batchJobSurfaceGeneratorFits.sh $MuonCorrection $ijob Data $date 5GeV 17
	                echo "MuonCorrection= $MuonCorrection, ijob= $ijob, Data, 5GeV"
	        done
	done
done
