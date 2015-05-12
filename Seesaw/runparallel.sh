#!/usr/bin/env bash

k=1
for i in 50 1050 2050 3050 4050 5050 6050
do
    if [ ! -d tmp$k ];then 
        mkdir tmp$k
    fi
    cd tmp$k
    ln -s ../main . 2> /dev/null
    ln -s ../schmidt_Neu_Res_Spread.py . 2> /dev/null
    ln -s ../schmidt_spread.py . 2> /dev/null
    cat schmidt_Neu_Res_Spread.py | sed -r 's/(mn1min=)[0-9]+/\1'$i'/' > kk.py 
    cat kk.py | sed -r 's/(mn1max=)[0-9]+/\1'$(($i+1000))'/' > kkk.py 
    cat kkk.py | sed -r 's/(mn1pts=)[0-9]+/\1 1000/' > kk.py 
    echo run $i $j... 
    python kk.py # > log$i &
    k=$(($k+1))
    cd ../
done


#echo -n "open pdfs? (y/n): "
#read open
#if [ "$open" = "y" ]; then
#    acroread $(ls R[0-9][0-9][0-0]t*_c.pdf)
#fi