#!/bin/bash

if [ "$1" == "R" ]
then
 echo "copying to bfl_RIKEN"
 scp -r .git bfl_RIKEN:~/LQCD/X_VOJTA_code/
 ssh bfl_RIKEN 'cd LQCD/X_VOJTA_code; git checkout .'
else
 echo "copying to bfl"
 scp -r .git bfl:~/LQCD/X_VOJTA_code/
 ssh bfl 'cd LQCD/X_VOJTA_code; git checkout .'
fi
