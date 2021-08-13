#3p2

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_3p2.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3P_2' at 100,20 font "Time-Roman,30"

chn = "3p2"
lmbd = "0800"
plot '../empirical_value/3p2.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:9 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:12 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:15 title 'UpToQn=4' w l lt 7 lw 3 
reset


#3f2

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_3f2.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3F_2' at 100,1.5 font "Time-Roman,30"

chn = "3p2"
lmbd = "0800"
plot '../empirical_value/3f2.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:10 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:13 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:16 title 'UpToQn=4' w l lt 7 lw 3 
reset

#e2

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_e2.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel 'Mixing angle [deg]'
set mytics 5
set nokey
set label '{/Symbol e}_2' at 100,-3 font "Time-Roman,30"

chn = "3p2"
lmbd = "0800"
plot '../empirical_value/e2.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:8 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:11 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:14 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:17 title 'UpToQn=4' w l lt 7 lw 3 
reset


#3p1

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_3p1.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3P_1' at 100,-25 font "Time-Roman,30"

chn = "3p1"
lmbd = "0800"
plot '../empirical_value/3p1.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:4 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:5 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=4' w l lt 7 lw 3 
reset

#1p1

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_1p1.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^1P_1' at 100,-20 font "Time-Roman,30"

chn = "1p1"
lmbd = "0800"
plot '../empirical_value/1p1.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:4 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:5 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=4' w l lt 7 lw 3 
reset

#1d2

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_1d2.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^1D_2' at 100,6 font "Time-Roman,30"

chn = "1d2"
lmbd = "0800"
plot '../empirical_value/1d2.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:4 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:5 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3d2

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_3d2.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3D_2' at 100,20 font "Time-Roman,30"

chn = "3d2"
lmbd = "0800"
plot '../empirical_value/3d2.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:4 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:5 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3d3

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_3d3.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3D_3' at 100,4 font "Time-Roman,30"

chn = "3d3"
lmbd = "0800"
plot '../empirical_value/3d3.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:9 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:12 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:15 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3g3

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_3g3.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3G_3' at 100,-2 font "Time-Roman,30"

chn = "3d3"
lmbd = "0800"
plot '../empirical_value/3g3.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:10 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:13 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:16 title 'UpToQn=4' w l lt 7 lw 3 
reset

#e3

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_e3.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel 'Mixing angle [deg]'
set mytics 5
set nokey
set label '{/Symbol e}_3' at 100,6 font "Time-Roman,30"

chn = "3d3"
lmbd = "0800"
plot '../empirical_value/e3.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:8 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:11 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:14 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:17 title 'UpToQn=4' w l lt 7 lw 3 
reset

#1f3

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_1f3.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^1F_3' at 100,-2 font "Time-Roman,30"

chn = "1f3"
lmbd = "0800"
plot '../empirical_value/1f3.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:4 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:5 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3f3

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_3f3.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3F_3' at 100,-1 font "Time-Roman,30"

chn = "3f3"
lmbd = "0800"
plot '../empirical_value/3f3.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:4 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:5 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3f4

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_3f4.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3F_4' at 100,2 font "Time-Roman,30"

chn = "3f4"
lmbd = "0800"
plot '../empirical_value/3f4.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:9 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:12 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:15 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3h4

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_3h4.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3H_4' at 100,0.35 font "Time-Roman,30"

chn = "3f4"
lmbd = "0800"
plot '../empirical_value/3h4.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:10 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:13 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:16 title 'UpToQn=4' w l lt 7 lw 3 
reset

#e4

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_e4.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel 'Mixing angle [deg]'
set mytics 5
set nokey
set label '{/Symbol e}_4' at 100,-1 font "Time-Roman,30"

chn = "3f4"
lmbd = "0800"
plot '../empirical_value/e4.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:8 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:11 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:14 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:17 title 'UpToQn=4' w l lt 7 lw 3 
reset

#1g4

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_1g4.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^1G_4' at 100,1.4 font "Time-Roman,30"

chn = "1g4"
lmbd = "0800"
plot '../empirical_value/1g4.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:4 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:5 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3g4

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_3g4.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3G_4' at 100,5 font "Time-Roman,30"

chn = "3g4"
lmbd = "0800"
plot '../empirical_value/3g4.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:4 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:5 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:7 title 'UpToQn=4' w l lt 7 lw 3 
reset

#3g5

set terminal postscript eps color enhanced solid "Time-Roman,30"
set output '../picture_saved/8_3g5.eps'
set xlabel 'K_c_._m_.[MeV]'
set xtics 100
set mxtics 4
set ylabel '{/Symbol d} [deg]'
set mytics 5
set nokey
set label '^3G_5' at 100,-0.8 font "Time-Roman,30"

chn = "3g5"
lmbd = "0800"
plot '../empirical_value/3g5.nij' using 1:3 title 'empirical value' pt 7 ps 1.5, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:6 title 'UpToQn=1' w l linestyle 0 lc 3 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:9 title 'UpToQn=2' with line dashtype 4 lc 2 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:12 title 'UpToQn=3' with line dashtype 5 lc 4 lw 3, 'set_8/run_'.chn.'_lmbd_'.lmbd.'.out' using 1:15 title 'UpToQn=4' w l lt 7 lw 3 
reset