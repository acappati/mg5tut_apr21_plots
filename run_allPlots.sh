#! /usr/bin/env sh

OPERATOR=("FM0" "FM1" "FM2" "FM3" "FM4" "FM5" "FM7" "FS0" "FS1" "FS2")
VALUES=("0" "minus20" "minus10" "minus5" "minus2" "2" "5" "10" "20")

for o in ${OPERATOR[@]}
do
    for v in ${VALUES[@]}
    do
        echo "running operator $o value $v ..."
        python3 plot_and_compute_fractions_checkCuts.py out_process1_ggTozzh $o $v 
    done
done
