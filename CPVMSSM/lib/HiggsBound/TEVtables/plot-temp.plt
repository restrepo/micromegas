#command: plot "CDF_SM_2_exp.txt CDF_SM_2_obs.txt" "" "" "" "" 
set timestamp
set term postscript portrait color solid
set output "plot-temp.ps" 



plot \
"CDF_SM_2_exp.txt"  wi li lw 4 ,\
"CDF_SM_2_obs.txt"  wi li lw 4 ,\
0 t "" wi li 7
