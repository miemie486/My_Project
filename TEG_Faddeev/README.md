# Installation

### Compilers

g++ and gfortran (7 or higher)

### LAPACK package

Depending on your OS, installation of LAPACK may vary a lot. If you are on a
Linux box, main-stream distros should include a package manager that installs
LAPACK straightforwardly. E.g., on Ubuntu, the following line can do the job:

$ sudo apt-get install libblas-dev liblapack-dev

Whatever OS you are running, you can always compile the library from source.
Download the source of LAPACK from

http://performance.netlib.org/lapack/

and compile it by following the instructions. In the end, you must have the
static lib files, typically named as liblapack.a, be placed in the directories
where the linker can find them. More on this later.

### Intel MKL package

At the following link,

https://software.intel.com/en-us/mkl/choose-download

choose the OS and download the MKL package. You must follow the package's
instructions to install it. Towards the end of installation, you are likely
to be notified of two scripts, looking like what follows,

/opt/intel/compilers_and_libraries_2019.1.144/linux/mkl/bin/mklvars.sh intel64
/opt/intel/compilers_and_libraries_2019.1.144/linux/bin/compilervars.sh intel64

Careful!!!  2019.1.144 is the version number of MKL, if you have install other MKL versions, you should change it.

You need to run them before compiling. Alternatively, you can put the following lines in in ~/.bashrc (Linux) or ~/.bash_profile (Mac)

source /opt/intel/compilers_and_libraries_2019.1.144/linux/mkl/bin/mklvars.sh intel64

source /opt/intel/compilers_and_libraries_2019.1.144/linux/bin/compilervars.sh intel64

so that they are automatically executed whenever a terminal is opened.

### Clone both repositories TEG_Faddeev and nnscat

Put both repos at the same level of directories

nnscat/

TEG_Faddeev/

### Compile the code

Navigate to src/ and run

$ make all

By looking into src/makefile, you can find a specific target to make, e.g.,

$ make objs

The source of executables is located in src/exec. To understand their usage,
please refer to the comments in the source.


#### Env vars TEGFC, TEGCXXC, and MACCXXC

It is difficult to have the source packaged in such a way that it just compiles
right out of the box. So one may need to tweak a lot to compile successfully.
make will use "g++" and "gfortran" to invoke C++/Fortran compilers. One needs
to make sure they point to the compilers with wanted versions. If you have
trouble directing make to find g++ and/or gfortran, you can use environment
variables TEGFC (gfortran) and TEGCXXC (g++) to specify the path to
the wanted GCC compilers. E.g., one can run, or put in .bashrc, the following
lines before compiling.

$ export TEGFC=/opt/gcc/bin/gfortran

$ export TEGCXXC=/opt/gcc/bin/gfortran

#### MAC

On MAC, make needs to be aware of mac osx so that appropriate MKL compiler flags can
be passed on to g++. This is done by setting environment variable

$ export MACCXXC=g++

or whatever path your system needs to access g++.

#### Env var TEGLIB

It may happen that with everything having been tried out, the linker still can
not find the required libraries, such as liblapack or libgfortran (used for C/
Fortran interface), even though the packages have been properly installed. If
the absolute path to those libs is known, you can specify it explicitly by

$ export TEGLIB=-L/path_to_lib/

-L is a linker option gcc used to direct the linker to a specific directory to
look for libraries.


*The following sections are no longer up to date as of Jun 6, 2020.*

# File system

Faddeev and related code created by the Effective Group at Sichuan University.

bin/  : Executable binaries
obj/  : Object files (\*.o)
src/  : Source code
inc/  : head files

src/kMat  : Caculate K_matrix
src/tMat  : Two-body t_matrix
src/exec  : executable programs

makefile.sh: make and run findEigen.exe. Use 'sh makefile.sh' commmand in terminal.

Note: Many codes are commended in programs. They may apear to be ugly but pretty useful for checking, so please don't touch them. I will delete them in appropriate time. --zeyuan

# Understanding kMat

The core function of kMat is finding the binding energy, which only takes two steps.

1. Initialization

Default initializing of kMatrix require one to input all necessary inputs, with the following consequences and data type:

kMatrix(int NP, int nChannels, REAL mambda, void (\*_tMatGenFunc)(int, int, int, int, REAL\*, REAL, int, REAL\*, REAL\*, REAL\*, int), REAL \*\*_channels)
where NP is the mesh points of the p momentum, which typically ranges from 12 to 50. Larger cutoff would lead to larger NP.
nChannels is the number of the channels.
mambda is the cutoff. The defination of REAL in included in TEG_Faddev/inc/def.h, which is double by default.
tMatGenFunc is a function which generate the t-matrix.
channels is a two dimentional array, with row the index of the channel and column the specific quantum number.
There are other parameters related to the numerical accuracy, such as NQ, NX ... They can be modified manually (see other inputs section), but we encourage you to use the default values, which is given by running the function
kMat.smartUpdate();


2. Finding the binding energy

After the initializing, one can find the binding energy by calling function

BE = kMatR.findBE(a, b, accuracy, "Secant");

where [a, b] is the interval where the binding energy lives in. If there is no binding energy in this interval, function findBE would return a or b. Note that a, b and BE are all negtive. Accuracy is the absolute error. "Secant" is the root finding method, alternative option is "Newton", which is slower.

3. Example

Here's a simple example

'''
#include "kMatGen.h"

using namespace std;

void tMatGenFunc(int intL, int intS, int intJ, int NQint, REAL* mEArray, REAL mambda, int NPint, REAL* pn, REAL* wpn, REAL* tMat, int NCH){
  get_FOST_tmtrx_DLSPR(&intL, &intS, &intJ, &NQint, mEArray, &mambda, &NPint, pn, wpn, tMat, &NCH);
} // Calling tMatrix from nnscat project

int main(){

  // Initialize 5 channels by function iniChannel
  int nChannels = 5;
  REAL **channels = new REAL*[nChannels];
  for(int i = 0; i < nChannels; i++){
    channels[i] = new REAL[COL_CHANNELs];
  }
  iniChannel(channels, nChannels);

  // Initialize k matrix, note all configurations can be changed later.
  int NP = 15;
  REAL mambda = 600;
  kMatrix kMatR {NP, nChannels, mambda, tMatGenFunc, channels};
  kMatR.smartUpdate();
  BE = kMatR.findBE(-5, -10, 0.01, "Secant");
'''

# Input: Channels

1. Routine
Channels is a two dimentional array. Every row contains all required quantum numbers.

'''
channel[row][] = {l, s, j, lambda, J, T, Coupled Channel}
'''

the meaning of those quantum numbers are listed below

'''
here    Physics meaning
l       Orbit momentem of the pairs
T       Isospin of the pairs
lambda  Orbit momentem of the third paticle
s       Spin of the pair
j       Total momentem of the pair
J       Total momentem of the third paticle
CC      Coupled channels' pair's orbit angular momentem. Negtive number means no coupling in this channel
opT     Total isospin = 1 / 2
'''

For example, the first two channels in Elster's notes is given by array

'''
{
  {0,   0,  0,  0,      .5,   1., -1},
  {0,   1,  1,  0,      .5,   0., 2}
}
'''

Details see channels.h and channels.cpp

2. Input method

The first one is use the default initializer mentioned in Understanding kMatrix section. Typical 6 channels can be initialized by function iniChannel. One can also write a two dimentional array manyally.

Alternatively, one can call function

  kMat.setChannel(REAL** _channels, int _nChannels)

Here's an example

'''
#include "kMatGen.h"

using namespace std;

void tMatGenFunc(int intL, int intS, int intJ, int NQint, REAL* mEArray, REAL mambda, int NPint, REAL* pn, REAL* wpn, REAL* tMat, int NCH){
  get_FOST_tmtrx_DLSPR(&intL, &intS, &intJ, &NQint, mEArray, &mambda, &NPint, pn, wpn, tMat, &NCH);
} // Calling tMatrix from nnscat project

int main(){

  // Initialize 5 channels by function iniChannel
  int nChannels = 5;
  REAL **channels = new REAL*[nChannels];
  for(int i = 0; i < nChannels; i++){
    channels[i] = new REAL[COL_CHANNELs];
  }

  // Initialize k matrix, note all configurations can be changed later.
  int NP = 15;
  REAL mambda = 600;
  kMatrix kMatR {NP, nChannels, mambda, tMatGenFunc, channels};

  iniChannel(channels, nChannels);
  kMatR.setChannel(channels, nChannels);
  kMatR.smartUpdate();
  BE = kMatR.findBE(-5, -10, 0.01, "Secant");
'''

With the help of class Chs (see channels.h), one can easily find all possible channels. The example is showed in TEG_Faddev/test/testTemp/testChs.cpp. Only three steps are important,

1. write the upper limit and parity

'''
  channels.writeMaxL(2);
  channels.writeMaxOpT(0.5);
  channels.writeOpJ(0.5);
  channels.writeParity(0);
'''

2. Generate all possible channels

  channels.genChs();

3. Output

Output the channels into a two dimentional array, which can be used by kMatrix. In this method, there's an additonal column which indicates the total isospin opT. Therefore the total size of the array is nChannels * 8, where nChannels is the required number of channels

channels.outChannels(double **outChs, int nChannels);

Print out the possible channels on the screen

  channels.showChannels();
