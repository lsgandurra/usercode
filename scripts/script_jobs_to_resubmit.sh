#!/bin/bash

## -------------------------------------------------------
## program implemented by Louis Sgandurra (October 2012)
## -------------------------------------------------------

cat $1 | awk '{print $2}' | sed 's/$/,/g' | sed ':a;N;$!ba;s/\n//g'
