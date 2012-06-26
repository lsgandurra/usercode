#! /usr/local/bin/bash -l
#$ -l ct=40000        ###time in seconds
#$ -P P_cmsf         
#$ -l vmem=4G
#$ -l fsize=30G
#$ -q long
#$ -l sps=1
###$ -l hpss=1
#$ -N FinalResults_Mai_2012_v3_SurfaceGen
### Merge the stdout et stderr in a single file
#$ -j y
### fichiers .e et .o copied to current working directory
#$ -cwd
###$ -m be
### set array job indices 'min-max:interval'
####$ -t 1-16


syntax="${0} {parameter}"
#if [[ -z ${6} ]]
if [[ -z ${2} ]]
then
        echo ${syntax}
        exit 1
fi

###ijob=`echo "${SGE_TASK_ID} - 1" | bc -ql`

echo "USER=${USER}"

# LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS
echo "LOAD CORRECT ENVIRONMENT VARIABLES FROM SPS"
###export HOMEDIR=/afs/in2p3.fr/home/o/obondu
###source ${HOMEDIR}/428v2.sh
export HOMEDIR=/afs/in2p3.fr/home/s/sgandurr
source ${HOMEDIR}/423.sh
SPSDIR=`pwd`
WORKDIR=${TMPDIR}

echo "USER=${USER}"

# CHECK THE ENVIRONMENT VARIABLES
echo "CHECK THE ENVIRONMENT VARIABLES"
echo "ROOTSYS :"
echo ${ROOTSYS}
echo ""

echo "USER=${USER}"

cd ${TMPDIR}/
echo "pwd; ls -als"
pwd; ls -als
echo ""


# COPY HEADER FILES TO WORKER
echo "COPY HEADER FILES TO WORKER"
#mkdir ${TMPDIR}/interface
#cp ${SPSDIR}/UserCode/IpnTreeProducer/interface/*h ${TMPDIR}/interface/
if [[ ! -e ${TMPDIR}/interface ]]
then
  mkdir ${TMPDIR}/interface
        cp ${SPSDIR}/UserCode/IpnTreeProducer/interface/*h ${TMPDIR}/interface/
fi

echo "USER=${USER}"

# COPY IpnTree LIB FILE TO WORKER
#mkdir ${TMPDIR}/lib
#cp ${SPSDIR}/UserCode/IpnTreeProducer/src/libToto.so ${TMPDIR}/lib/
echo "COPY IpnTree LIB FILE TO WORKER"
if [[ ! -e ${TMPDIR}/lib ]]
then
        mkdir ${TMPDIR}/lib
        cp ${SPSDIR}/UserCode/IpnTreeProducer/src/libToto.so ${TMPDIR}/lib/
fi

echo "USER=${USER}"

# ADD CURRENT DIRECTORY AND LIB TO LIRARY PATH
LD_LIBRARY_PATH=`echo "${LD_LIBRARY_PATH}:${WORKDIR}/lib:${WORKDIR}"`
echo "LD_LIBRARY_PATH"
echo ${LD_LIBRARY_PATH}
echo ""

echo "USER=${USER}"


# COPY EXECUTABLE TO WORKER
echo "COPY EXECUTABLE TO WORKER"
cp ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/mmg_ik_MZ_Surface_Fits.exe ${TMPDIR}/
cp ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/CrystalBall.C ${TMPDIR}/
cp ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/setTDRStyle.C ${TMPDIR}/
cp ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/CMSStyle.C ${TMPDIR}/
cp ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/*.txt ${TMPDIR}/
##cp /sps/cms/sgandurr/CMSSW_4_2_8/src/UserCode/IpnTreeProducer/OlivierMiniTrees/*.root ${TMPDIR}/
##cp /sps/cms/sgandurr/CMSSW_4_2_8/src/UserCode/IpnTreeProducer/ZMuMuGammaMinitrees_Mars_2012/*.root ${TMPDIR}/
cp /sps/cms/sgandurr/CMSSW_4_2_8/src/UserCode/IpnTreeProducer/*.root ${TMPDIR}/

###cp ${SPSDIR}/UserCode/IpnTreeProducer/ComparaisonFits1overKrecoDiffMmumuJan/mmg_ik_MZ_Surface_Fits.exe ${TMPDIR}/
###cp ${SPSDIR}/UserCode/IpnTreeProducer/ComparaisonFits1overKrecoDiffMmumuJan/CrystalBall.C ${TMPDIR}/
###cp ${SPSDIR}/UserCode/IpnTreeProducer/ComparaisonFits1overKrecoDiffMmumuJan/setTDRStyle.C ${TMPDIR}/
###cp ${SPSDIR}/UserCode/IpnTreeProducer/ComparaisonFits1overKrecoDiffMmumuJan/CMSStyle.C ${TMPDIR}/


echo "pwd; ls -als"
pwd; ls -als
echo ""

echo "USER=${USER}"

# EXECUTE JOB
echo "EXECUTE JOB"
cd ${TMPDIR}/


for Endcaps in `seq 0 1`
do
        for r9sup in `seq 0 2`
        do
                for SetOfCorrections in 'ETHZCorrections' ##'ElectronTunedCorrections'
                do
                        for phiCracks in '1' ##'0' '1'
                        do
                                for etaCracks in '1' ##'0' '1'
                                do     
                                        for variableX in 'Photon_Et' ##'Photon_SC_rawEt' 'Photon_Et' 'Photon_E' 'Photon_SC_Eta' 'Photon_SC_brem' 
                                        do
						for MuonCorrection in 'NoMuonCorrection' 'Rochester'
						do 
							for Category in 'Vgamma24' ##'Vgamma8' 'OneBin' 'Vgamma24'
							do					
		###						for HighMmumuLim in '76' '84'
		###			                        do
		###							echo "$Endcaps $r9sup $SetOfCorrections ${3} 2011 36 $HighMmumuLim ${2} $phiCracks $etaCracks $variableX ${1}"
		###							./mmg_ik_MZ_Surface_Fits.exe $Endcaps $r9sup $SetOfCorrections ${3} 2011 36 $HighMmumuLim ${2} $phiCracks $etaCracks $variableX ${1} >> sortie_${Endcaps}_${r9sup}_${SetOfCorrections}_${3}_2011_36_${HighMmumuLim}_${2}_${phiCracks}_${etaCracks}_${variableX}_${1}.out 2> sortie_${Endcaps}_${r9sup}_${SetOfCorrections}_${3}_2011_36_${HighMmumuLim}_${2}_${phiCracks}_${etaCracks}_${variableX}_${1}.err
		###							##cp ${TMPDIR}/*.err ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/
		###							##cp ${TMPDIR}/*.out ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/
		###
		###			                        done
		###			                        for HighMmumuLim in '78' '82' '84' '86' '88'
		###			                        do
		###							echo "$Endcaps $r9sup $SetOfCorrections ${3} 2011 38 $HighMmumuLim ${2} $phiCracks $etaCracks $variableX ${1}"
		###							./mmg_ik_MZ_Surface_Fits.exe $Endcaps $r9sup $SetOfCorrections ${3} 2011 38 $HighMmumuLim ${2} $phiCracks $etaCracks $variableX ${1} >> sortie_${Endcaps}_${r9sup}_${SetOfCorrections}_${3}_2011_38_${HighMmumuLim}_${2}_${phiCracks}_${etaCracks}_${variableX}_${1}.out 2> sortie_${Endcaps}_${r9sup}_${SetOfCorrections}_${3}_2011_38_${HighMmumuLim}_${2}_${phiCracks}_${etaCracks}_${variableX}_${1}.err
		###			     
		###			                        done
					                        for HighMmumuLim in '80'
					                        do
									echo "$Endcaps $r9sup $SetOfCorrections ${3} 2011 40 $HighMmumuLim ${2} $phiCracks $etaCracks $variableX ${1} ${4} ${MuonCorrection} ${5} ${Category} ${6}"
									./mmg_ik_MZ_Surface_Fits.exe $Endcaps $r9sup $SetOfCorrections ${3} 2011 40 $HighMmumuLim ${2} $phiCracks $etaCracks $variableX ${1} ${4} ${MuonCorrection} ${5} ${Category} ${6} >> sortie_${Endcaps}_${r9sup}_${SetOfCorrections}_${3}_2011_40_${HighMmumuLim}_${2}_${phiCracks}_${etaCracks}_${variableX}_${1}_${4}_${MuonCorrection}_${5}_${Category}_${6}.out 2> sortie_${Endcaps}_${r9sup}_${SetOfCorrections}_${3}_2011_40_${HighMmumuLim}_${2}_${phiCracks}_${etaCracks}_${variableX}_${1}_${4}_${MuonCorrection}_${5}_${Category}_${6}.err
					    			##./mmg_ik_MZ_Surface_Fits.exe 0 1 ETHZCorrections RooBifurcatedGauss 2011 40 80 89 1 1 Photon_Et 1 05GeV ${i} 
					                        done
		###			                        for HighMmumuLim in '82'
		###			                        do
		###							echo "$Endcaps $r9sup $SetOfCorrections ${3} 2011 42 $HighMmumuLim ${2} $phiCracks $etaCracks $variableX ${1}"
		###							./mmg_ik_MZ_Surface_Fits.exe $Endcaps $r9sup $SetOfCorrections ${3} 2011 42 $HighMmumuLim ${2} $phiCracks $etaCracks $variableX ${1} >> sortie_${Endcaps}_${r9sup}_${SetOfCorrections}_${3}_2011_42_${HighMmumuLim}_${2}_${phiCracks}_${etaCracks}_${variableX}_${1}.out 2> sortie_${Endcaps}_${r9sup}_${SetOfCorrections}_${3}_2011_42_${HighMmumuLim}_${2}_${phiCracks}_${etaCracks}_${variableX}_${1}.err
		###			
		###			                        done
		###			                        for HighMmumuLim in '84'
		###			                        do
		###							echo "$Endcaps $r9sup $SetOfCorrections ${3} 2011 44 $HighMmumuLim ${2} $phiCracks $etaCracks $variableX ${1}"
		###							./mmg_ik_MZ_Surface_Fits.exe $Endcaps $r9sup $SetOfCorrections ${3} 2011 44 $HighMmumuLim ${2} $phiCracks $etaCracks $variableX ${1} >> sortie_${Endcaps}_${r9sup}_${SetOfCorrections}_${3}_2011_44_${HighMmumuLim}_${2}_${phiCracks}_${etaCracks}_${variableX}_${1}.out 2> sortie_${Endcaps}_${r9sup}_${SetOfCorrections}_${3}_2011_44_${HighMmumuLim}_${2}_${phiCracks}_${etaCracks}_${variableX}_${1}.err
		###			                        done
							done
						done	
					done
				done
			done
		done
	done
done






echo "pwd; ls -als"
pwd; ls -als
echo ""

# GET BACK OUTPUT FILES TO SPS
echo "GET BACK OUTPUT FILES TO SPS AND REMOVE THEM FROM DISTANT DIR"
mkdir ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/Results_Mai_2012_v3_SurfaceGen_err/
mkdir ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/Results_Mai_2012_v3_SurfaceGen_out/
mv ${TMPDIR}/sortie_*_*_*_${3}_2011_40_*_${2}_*_*_*_${1}_${4}_*_${5}_*_${6}.err ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/Results_Mai_2012_v3_SurfaceGen_err/
mv ${TMPDIR}/sortie_*_*_*_${3}_2011_40_*_${2}_*_*_*_${1}_${4}_*_${5}_*_${6}.out ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/Results_Mai_2012_v3_SurfaceGen_out/
cp -r ${TMPDIR}/Results_Mai_2012_v3_SurfaceGen/ ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/
rm -rf Results_Mai_2012_v3_SurfaceGen/
rm ${TMPDIR}/mmg_ik_MZ_Surface_Fits.exe
rm ${TMPDIR}/CrystalBall.C
rm ${TMPDIR}/setTDRStyle.C
rm ${TMPDIR}/CMSStyle.C
rm ${TMPDIR}/*.txt
rm ${TMPDIR}/*.root

#"cd ${SPSDIR}/UserCode/IpnTreeProducer/OlivierMiniTrees/"
#cd ${SPSDIR}/UserCode/IpnTreeProducer/OlivierMiniTrees/
#echo "pwd; ls -als"
#pwd; ls -als
#echo ""

exit 0


