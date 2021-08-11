source /opt/intel/compilers_and_libraries_2020.0.166/linux/mkl/bin/mklvars.sh intel64
source /opt/intel/compilers_and_libraries_2020.0.166/linux/bin/compilervars.sh intel64

cd "/home/miemie/My_Project/lapack_exercise/" && g++ -llapacke -fopenmp det.cpp -o det && "/home/miemie/My_Project/lapack_exercise/"det > data1.txt
sort -g data1.txt > data.txt
#rm -f data1.txt
python -u "/home/miemie/My_Project/lapack_exercise/plot.py"