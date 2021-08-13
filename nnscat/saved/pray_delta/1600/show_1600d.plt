
#1d2

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '11_1d2.eps'
set xlabel 'setn=11 k_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^1D_2' at 100,6 font "Time-Roman,30"

chn = "1d2"
lmbd = "1600"
plot '../empirical_value/1d2.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:4 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:5 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3d2

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '11_3d2.eps'
set xlabel 'setn=11 k_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3D_2' at 100,20 font "Time-Roman,30"

chn = "3d2"
lmbd = "1600"
plot '../empirical_value/3d2.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:4 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:5 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3d3

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '11_3d3.eps'
set xlabel 'setn=11 k_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3D_3' at 100,4 font "Time-Roman,30"

chn = "3d3"
lmbd = "1600"
plot '../empirical_value/3d3.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:9 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:12 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:15 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3g3

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '11_3g3.eps'
set xlabel 'setn=11 k_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3G_3' at 100,-2 font "Time-Roman,30"

chn = "3d3"
lmbd = "1600"
plot '../empirical_value/3g3.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:10 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:13 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:16 title 'UpToQn=4' w l lt 7 lw 3 
reset

#e3

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '11_e3.eps'
set xlabel 'setn=11 k_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel 'Mixing angle [deg]'
set mytics 5
set nokey
set label '{/Symbol e}_3' at 100,6 font "Time-Roman,30"

chn = "3d3"
lmbd = "1600"
plot '../empirical_value/e3.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:8 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:11 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:14 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_11/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:17 title 'UpToQn=4' w l lt 7 lw 3 
reset


#1d2

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '5_1d2.eps'
set xlabel 'setn=5 k_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^1D_2' at 100,6 font "Time-Roman,30"

chn = "1d2"
lmbd = "1600"
plot '../empirical_value/1d2.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:4 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:5 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3d2

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '5_3d2.eps'
set xlabel 'setn=5 k_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3D_2' at 100,20 font "Time-Roman,30"

chn = "3d2"
lmbd = "1600"
plot '../empirical_value/3d2.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:4 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:5 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3d3

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '5_3d3.eps'
set xlabel 'setn=5 k_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3D_3' at 100,4 font "Time-Roman,30"

chn = "3d3"
lmbd = "1600"
plot '../empirical_value/3d3.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:9 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:12 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:15 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3g3

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '5_3g3.eps'
set xlabel 'setn=5 k_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3G_3' at 100,-2 font "Time-Roman,30"

chn = "3d3"
lmbd = "1600"
plot '../empirical_value/3g3.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:10 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:13 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:16 title 'UpToQn=4' w l lt 7 lw 3 
reset

#e3

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '5_e3.eps'
set xlabel 'setn=5 k_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel 'Mixing angle [deg]'
set mytics 5
set nokey
set label '{/Symbol e}_3' at 100,6 font "Time-Roman,30"

chn = "3d3"
lmbd = "1600"
plot '../empirical_value/e3.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:8 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:11 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:14 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_5/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:17 title 'UpToQn=4' w l lt 7 lw 3 
reset


#1d2

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '8_1d2.eps'
set xlabel 'setn=8 k_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^1D_2' at 100,6 font "Time-Roman,30"

chn = "1d2"
lmbd = "1600"
plot '../empirical_value/1d2.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:4 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:5 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3d2

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '8_3d2.eps'
set xlabel 'setn=8 k_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3D_2' at 100,20 font "Time-Roman,30"

chn = "3d2"
lmbd = "1600"
plot '../empirical_value/3d2.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:4 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:5 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3d3

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '8_3d3.eps'
set xlabel 'setn=8 k_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3D_3' at 100,4 font "Time-Roman,30"

chn = "3d3"
lmbd = "1600"
plot '../empirical_value/3d3.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:9 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:12 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:15 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3g3

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '8_3g3.eps'
set xlabel 'setn=8 k_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3G_3' at 100,-2 font "Time-Roman,30"

chn = "3d3"
lmbd = "1600"
plot '../empirical_value/3g3.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:10 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:13 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:16 title 'UpToQn=4' w l lt 7 lw 3 
reset

#e3

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '8_e3.eps'
set xlabel 'setn=8 k_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel 'Mixing angle [deg]'
set mytics 5
set nokey
set label '{/Symbol e}_3' at 100,6 font "Time-Roman,30"

chn = "3d3"
lmbd = "1600"
plot '../empirical_value/e3.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:8 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:11 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:14 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:17 title 'UpToQn=4' w l lt 7 lw 3 
reset

