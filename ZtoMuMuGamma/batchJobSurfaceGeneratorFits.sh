#! /usr/local/bin/bash -l
#$ -l ct=30000        ###time in seconds
#$ -P P_cmsf         
#$ -l vmem=4G
#$ -l fsize=30G
#$ -q long
#$ -l sps=1
###$ -l hpss=1
#$ -N FinalSurfacePlots_v2
### Merge the stdout et stderr in a single file
#$ -j y
### fichiers .e et .o copied to current working directory
#$ -cwd
###$ -m be
### set array job indices 'min-max:interval'
###$ -t 1-18


syntax="${0} {parameter}"
#if [[ -z ${6} ]]
if [[ -z ${1} ]]
then
        echo ${syntax}
        exit 1
fi

##ijob=`echo "${SGE_TASK_ID} - 1" | bc -ql`

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
mkdir Dossier_${3}_${1}_${4}_${5}_${2}
cd Dossier_${3}_${1}_${4}_${5}_${2}
echo "pwd; ls -als"
pwd; ls -als
echo ""


# COPY HEADER FILES TO WORKER
echo "COPY HEADER FILES TO WORKER"
#mkdir ${TMPDIR}/interface
#cp ${SPSDIR}/UserCode/IpnTreeProducer/interface/*h ${TMPDIR}/interface/
if [[ ! -e ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/interface ]]
then
  mkdir ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/interface
        cp ${SPSDIR}/UserCode/IpnTreeProducer/interface/*h ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/interface/
fi

echo "USER=${USER}"

# COPY IpnTree LIB FILE TO WORKER
#mkdir ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/lib
#cp ${SPSDIR}/UserCode/IpnTreeProducer/src/libToto.so ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/lib/
echo "COPY IpnTree LIB FILE TO WORKER"
if [[ ! -e ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/lib ]]
then
        mkdir ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/lib
        cp ${SPSDIR}/UserCode/IpnTreeProducer/src/libToto.so ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/lib/
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
cp ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/SurfaceGeneratorFit_v2.exe ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/
cp ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/CrystalBall.C ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/
cp ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/setTDRStyle.C ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/
cp ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/CMSStyle.C ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/
cp ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/*.txt ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/
##cp /sps/cms/sgandurr/CMSSW_4_2_8/src/UserCode/IpnTreeProducer/OlivierMiniTrees/*.root ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/
##cp /sps/cms/sgandurr/CMSSW_4_2_8/src/UserCode/IpnTreeProducer/ZMuMuGammaMinitrees_Mars_2012/*.root ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/
cp /sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v08*.root ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/
cp /sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v09*.root ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/
cp /sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/Selection/miniTreeMuons_v10*.root ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/


###cp ${SPSDIR}/UserCode/IpnTreeProducer/ComparaisonFits1overKrecoDiffMmumuJan/SurfaceGeneratorFit_v2.exe ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/
###cp ${SPSDIR}/UserCode/IpnTreeProducer/ComparaisonFits1overKrecoDiffMmumuJan/CrystalBall.C ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/
###cp ${SPSDIR}/UserCode/IpnTreeProducer/ComparaisonFits1overKrecoDiffMmumuJan/setTDRStyle.C ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/
###cp ${SPSDIR}/UserCode/IpnTreeProducer/ComparaisonFits1overKrecoDiffMmumuJan/CMSStyle.C ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/


echo "pwd; ls -als"
pwd; ls -als
echo ""

echo "USER=${USER}"

# EXECUTE JOB
echo "EXECUTE JOB"
cd ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/

	##./SurfaceGeneratorFit_v2.exe MC ${1} Fall11 Mmumu.txt 5GeV RooCrystalBall 18 ${ijob} >> SurfaceGeneratorFit_v2_MC_${1}_Fall11_Mmumu_5GeV_CB_${ijob}.out 2> SurfaceGeneratorFit_v2_MC_${1}_Fall11_Mmumu_5GeV_CB_${ijob}.err
./SurfaceGeneratorFit_v2.exe ${3} ${1} ${4} Mmumu.txt ${5} RooCrystalBall ${6} ${2} >> SurfaceGeneratorFit_v2_${3}_${1}_${4}_Mmumu_${5}_CB_${2}.out 2> SurfaceGeneratorFit_v2_${3}_${1}_${4}_Mmumu_${5}_CB_${2}.err


echo "pwd; ls -als"
pwd; ls -als
echo ""

# GET BACK OUTPUT FILES TO SPS
echo "GET BACK OUTPUT FILES TO SPS AND REMOVE THEM FROM DISTANT DIR"
mkdir ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/SurfacePlots_v2_err/
mkdir ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/SurfacePlots_v2_out/
##mv SurfaceGeneratorFit_v2_MC_${1}_Fall11_Mmumu_5GeV_CB_${ijob}.err ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/SurfacePlots_v2_err/
##mv SurfaceGeneratorFit_v2_MC_${1}_Fall11_Mmumu_5GeV_CB_${ijob}.out ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/SurfacePlots_v2_out/
mv SurfaceGeneratorFit_v2_${3}_${1}_${4}_Mmumu_${5}_CB_${2}.err ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/SurfacePlots_v2_err/
mv SurfaceGeneratorFit_v2_${3}_${1}_${4}_Mmumu_${5}_CB_${2}.out ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/SurfacePlots_v2_out/
cp -r ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/SurfacePlots_v2/ ${SPSDIR}/UserCode/IpnTreeProducer/mmg_ik_MZ_Surface/
rm -rf SurfacePlots_v2/
rm ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/SurfaceGeneratorFit_v2.exe
rm ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/CrystalBall.C
rm ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/setTDRStyle.C
rm ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/CMSStyle.C
rm ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/*.txt
rm ${TMPDIR}/Dossier_${3}_${1}_${4}_${5}_${2}/*.root

#"cd ${SPSDIR}/UserCode/IpnTreeProducer/OlivierMiniTrees/"
#cd ${SPSDIR}/UserCode/IpnTreeProducer/OlivierMiniTrees/
#echo "pwd; ls -als"
#pwd; ls -als
#echo ""

exit 0


