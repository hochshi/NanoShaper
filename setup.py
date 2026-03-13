import subprocess,os
import platform
import os
import shutil
import sys
my_env = os.environ.copy()
# Universal NanoShaper installer for so and stand-alone 
# version 1.5

print ("--------------------------------------") 
print ("Welcome to NanoShaper Setup")
print ("--------------------------------------")
print ("This script will setup NanoShaper executable/so module for Linux.")
print ("At the end you can check the build results in make_ns.txt.")
print ("")
print ("For any help contact sergio.decherchi@iit.it")
print ("")
print ("Pre-requisites are: gmp, mpfr, boost Threads/FileSystem/Chrono")
print ("together with a working Internet connection.")
print("If gmp/mpfr/boost are not in standard locations you can set")
print("the following env variables (before running this script) to provide custom locations:")
print("MPFR_LIBRARIES path to mpfr .so (excluding file name)")
print("MPFR_INCLUDE_DIR for mpfr include dir")
print("BOOST_LIBRARYDIR path to boost .so files")
print("BOOST_INCLUDEDIR for boost include dir")
print("GMP_LIBRARIES path for gmp .so file (excluding file name)")
print("GMP_INCLUDE_DIR for gmp include dir")

print ("")

################### cmake, tbb, cgal versions ##################
cgalVersion = '5.6.2'
cgalFile = 'cgal-'+cgalVersion+'.tar.gz'
# relative path to NS root
cgalDir = 'cgal-'+cgalVersion
cmakeVersion = "3.22.0-rc2" 
cmake_name = "cmake-%s-linux-x86_64"%cmakeVersion
tbb_ver_short = '2021.4.0'
tbb_ver = 'oneapi-tbb-%s'%(tbb_ver_short)
#################################################################

compMode = 0
str = 'exe'
if sys.version_info[0] < 3:
  str = raw_input('Do you want to compile NanoShaper as a shared object or as a stand-alone executable? [so/exe]').lower()
else:
  str = input('Do you want to compile NanoShaper as a shared object or as a stand-alone executable? [so/exe]').lower()

if (str=='exe'):
  print ("Stand alone mode selected")
  compMode = 0
elif (str=='so'):
  print ("Shared object mode selected")
  compMode = 2
else:
  print ("Option not recognised, assuming stand-alone executable")

pwd = os.getcwd()
cmake_path = pwd + '/' + cmake_name
tbb_path = pwd + '/%s'%(tbb_ver)

#CMAKE (at least 3.1 for tbb and at least 3.14 for cgal)
if(os.path.isdir(cmake_path)==False):
  subprocess.call('wget https://github.com/Kitware/CMake/releases/download/v%s/%s.tar.gz'%(cmakeVersion,cmake_name),shell=True) 
  subprocess.call('tar -vxzf %s.tar.gz'%cmake_name,shell=True)
  print ("")

#TBB 
if(os.path.isdir(tbb_ver)==False):
  subprocess.call('wget https://github.com/oneapi-src/oneTBB/releases/download/v%s/%s-lin.tgz'%(tbb_ver_short,tbb_ver),shell=True)    
  subprocess.call('tar -xvf %s-lin.tgz'%(tbb_ver),shell=True)    

#CGAL
if(os.path.isdir(cgalDir)==False):
  subprocess.call('wget https://github.com/CGAL/cgal/archive/v'+cgalVersion+'.tar.gz',shell=True) 
  subprocess.call('sync',shell=True)
  subprocess.call('tar -xvf v'+cgalVersion+'.tar.gz',shell=True)
  subprocess.call('sync',shell=True)
  subprocess.call('rm v'+cgalVersion+'.tar.gz',shell=True)

  print ("")
  print ("Building CGAL...")
  print ("")
  
build_cgal = "cd ./%s \n mkdir build \n cd build \n rm -f CMakeCache.txt \n sync \n"%(cgalDir) 

temp = cmake_path+'/bin/cmake .. -DBUILD_SHARED_LIBS=ON  -DWITH_examples=false  -DWITH_CGAL_Qt5=false  -DWITH_CGAL_ImageIO=false -DCMAKE_BUILD_TYPE="Release" '

build_cgal2 = ""

# add custom env variables if present in the environment
if ("MPFR_LIBRARIES" in my_env):
	build_cgal2 = build_cgal2+ " -DMPFR_LIBRARIES=$MPFR_LIBRARIES/libmpfr.so"
if ("MPFR_INCLUDE_DIR" in my_env):
	build_cgal2 = build_cgal2+ " -DMPFR_INCLUDE_DIR=$MPFR_INCLUDE_DIR"
if ("BOOST_INCLUDEDIR" in my_env):
	build_cgal2 = build_cgal2+ " -DBoost_INCLUDE_DIR=$BOOST_INCLUDEDIR"
if ("BOOST_LIBRARYDIR" in my_env):
	build_cgal2 = build_cgal2+ " -DBOOST_LIBRARYDIR=$BOOST_LIBRARYDIR"
if ("GMP_LIBRARIES" in my_env):
	build_cgal2 = build_cgal2+ " -DGMP_LIBRARIES=$GMP_LIBRARIES/libgmp.so"
if ("GMP_INCLUDE_DIR" in my_env):
	build_cgal2 = build_cgal2+ " -DGMP_INCLUDE_DIR=$GMP_INCLUDE_DIR"
  
build_cgal = build_cgal+temp+build_cgal2+"\n sync \n cd .. \n cd .."
 
print("Cmake command cgal")
print(build_cgal)
subprocess.call(build_cgal,env=my_env,shell=True)

print ("")

where = ''

if (compMode==0):
  print ("Building NanoShaper Stand-Alone please wait...")
  shutil.copy2('CMakeLists_standalone.txt','CMakeLists.txt')
  where = 'build'
elif (compMode==2):
  print ("Building NanoShaper shared object please wait...")
  shutil.copy2('CMakeLists_so.txt','CMakeLists.txt')
  where = 'build_so'
        
subprocess.call('mkdir {}\n'.format(where),shell=True) 
print ("")

build_ns = """export NS_ROOT=$PWD
cd %s
rm -fr * """%(where)
build_ns = build_ns+ '\nexport TBB_DIR='+tbb_path+'/lib/cmake/tbb \n %s/bin/cmake .. -G Ninja -DCGAL_DIR=$NS_ROOT/%s -DCMAKE_BUILD_TYPE="Release" '%(cmake_path,cgalDir)

# add custom env variables if present in the environment
if ("MPFR_LIBRARIES" in my_env):
	build_ns = build_ns+ " -DMPFR_LIBRARIES=$MPFR_LIBRARIES/libmpfr.so"
if ("MPFR_INCLUDE_DIR" in my_env):
	build_ns = build_ns+ " -DMPFR_INCLUDE_DIR=$MPFR_INCLUDE_DIR"
if ("BOOST_INCLUDEDIR" in my_env):
	build_ns = build_ns+ " -DBOOST_INCLUDEDIR=$BOOST_INCLUDEDIR"
if ("BOOST_LIBRARYDIR" in my_env):
	build_ns = build_ns+ " -DBOOST_LIBRARYDIR=$BOOST_LIBRARYDIR"
if ("GMP_LIBRARIES" in my_env):
	build_ns = build_ns+ " -DGMP_LIBRARIES=$GMP_LIBRARIES/libgmp.so"
if ("GMP_INCLUDE_DIR" in my_env):
	build_ns = build_ns+ " -DGMP_INCLUDE_DIR=$GMP_INCLUDE_DIR"

build_ns = build_ns + "\n make clean \n make -j 8"

print("Cmake command")
print(build_ns)

subprocess.call(build_ns,shell=True,env=my_env)
  
if (compMode==1):

  print ("For using the lib with DelPhi please assure that the lib is reachable by the OS")
  print ("Once you have root privileges you can copy it on /usr/lib or /usr/lib64 folder.")
  print ("If you have not root privileges you can update your LD_LIBRARY_PATH to the path where the lib is located")
  
print ("")
print ("")
