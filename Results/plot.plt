set term pngcairo
set output 'test10.png'
set autoscale cbfix
set pm3d map
splot 'map.val' u 1:2:3
