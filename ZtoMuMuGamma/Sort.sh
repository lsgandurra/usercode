#./bin/bash

echo ${1}

for i in `seq 0 ${2}` 
do 
 	echo MmumuFit_${i}.txt >> ${1}liste.txt 
done


