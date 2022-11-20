#!/bin/bash

# Put here the program (maybe with path)
umbralPot=$1
umbralZCR=$2
umbralR1=$3
umbralRMax=$4




if [ -z "$umbralPot" ]
then
    umbralPot=-45
fi
if [ -z "$umbralZCR" ]
then
    umbralZCR=164
fi

if [ -z "$umbralR1" ]
then
    umbralR1=0.52
fi

if [ -z "$umbralRMax" ]
then 
    umbralRMax=0.396
fi


GETF0="/home/alejandro/PAV/bin/get_pitch -p $umbralPot -z $umbralZCR -n $umbralR1 -u $umbralRMax"
for fwav in pitch_db/train/*.wav; do
    ff0=${fwav/.wav/.f0}
    echo "$GETF0 $fwav $ff0 ----"
    $GETF0 $fwav $ff0 > /dev/null || (echo "Error in $GETF0 $fwav $ff0"; exit 1)
done




exit 0
