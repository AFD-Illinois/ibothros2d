#!/bin/zsh
n=0
for ((a = 0. ; a <= 1. ; a+=0.01))
do
	((n++))
	typeset -Z3 n
	echo $n
	./mibothros $a 1.57
	mv shadow.ppm shadow_$n.ppm
done
for ((i = 1.57 ; i > 0. ; i-=0.01))
do
	((n++))
	typeset -Z3 n
	echo $n
	./mibothros $a $i
	mv shadow.ppm shadow_$n.ppm
done
