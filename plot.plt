set term pngcairo
set output 'test4.png'
set autoscale cbfix
set pm3d map
splot 'map.val' matrix notitle
