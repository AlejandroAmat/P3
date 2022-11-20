#!/bin/bash

for umbralRMax in $(seq 0 1 10);do
    for umbralZCR in $(seq 154 1 154);do
        for umbralR1 in $(seq 0.82 1 0.82);do
            echo -n "umbralPot=-45 umbralZCR=$umbralZCR umbralR1=$umbralR1 umbralRMax=$umbralRMax "
            ./scripts/run_get_pitch.sh -45 162 0.52 0.396  > /dev/null
             ~/PAV/bin/pitch_evaluate pitch_db/train/*f0ref | fgrep TOTAL
             
        done
    done 
done | sort -t: -k 2n;
exit 0