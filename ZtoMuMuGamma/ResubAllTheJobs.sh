#!/usr/local/bin/bash
# Delete jobs on ccage
# Written by Louis Sgandurra (January 2012)

for file in `qstat | awk '{print $1}' | tail -n +3`
do
        qresub "$file"
	qdel "$file"
	##echo "$file"
done

