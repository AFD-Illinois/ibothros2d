#!/bin/sh
#
rm -rf images
mkdir images
#
imax=100
for ((i = 0 ; i < $imax ; i++))
do
	printf -v n "%03d" $i
	echo $n
	# convert n to the x2 (colatitude-like) variable
	x2=$(bc -l <<< "($n+0.5)*1.0/$imax")
	echo $x2
	#
	./ibothros2d $x2 230.e9 ../iharm2d_v3/prob:kerr_torus/dumps/dump040 6.e18
	mv ibothros2d_fnu.ppm images/fnu.$n.ppm
	mv ibothros2d_lfnu.ppm images/lfnu.$n.ppm
	mv ibothros2d.dat images/fnu.$n.dat
done

