#!/bin/sh
#
# .tga files are needed by mencoder; mencoder is part
# of mplayer package
#
# possible to resize, label figures, etc.  remove resize
# if you want the full resolution of original figures
#
# here images files are assumed to have the form ang.000.ppm,
# ang.001.ppm, # etc.
#
for fil in `/bin/ls shadow*.ppm`
do
	echo $fil 
	convert $fil s_$fil.tga
done
#
# ffv1 is ffmpeg internal lossles codec
# this makes movie run 3 times, once fast and then 2 times slowly
#
mencoder "mf://*.tga" -mf fps=30 "mf://*.tga" -mf fps=10 "mf://*.tga" \
	-mf fps=10 -o tmp.avi -ovc lavc -lavcopts vcodec=ffv1
#
# 
# sameq gives same quality (lossless) as original.
# remove to reduce size of movie -- which will look less good.
ffmpeg -sameq -i tmp.avi 2d_shbox_ang_small.mpg
#
