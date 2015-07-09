#!/bin/zsh
for ((th = 0.01 ; th < 3.14159 ; th += 0.01))
do
	typeset -F 3 th
	echo $th
	./ibothros2d $th 230.e9 /home/gammie/iHARM2d_021613/highres/dumps/dump030 6.e18
	mv ibothros2d_fnu.ppm fnu.$th.ppm
	mv ibothros2d_lfnu.ppm lfnu.$th.ppm
done

