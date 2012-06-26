#!/bin/bash
syntax="Syntax ${0} {executableName}"
if [[ -z ${1} ]]
then
  echo ${syntax}
#  exit 1
fi
exeFile=${1}


##g++ mmg_ik_MZ_Surface_Fits.C -L${ROOTSYS}lib -lTMVA -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMLP -lTreePlayer -lMinuit -pthread -lm -ldl -rdynamic -pthread -m64 -I${ROOTSYS}include -L`pwd` `root-config --libs --cflags` -o ${1}

g++ ${1} -L`pwd` -lRooFitCore -lRooFit `root-config --libs --cflags` -o ${2}


