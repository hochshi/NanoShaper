import subprocess
import platform
import os
import shutil

# Universal NanoShaper installer for lib, stand-alone and python module
# version 0.7

print "--------------------------------------"
print "Welcome to NanoShaper Setup"
print "--------------------------------------"
print "This script will setup NanoShaper executable/lib/python module for Linux/Mac."
print "For Windows pre-compiled binaries are available."
print "At the end you can check the build results in make_ns.txt."
print ""
print "For any help contact sergio.decherchi@iit.it"
print ""
print "If you additionally have root privileges this script can install for you also the missing packages"
print "and if NanoShaper is compiled as a library/python module this will be installed. "
print "To get the full automation of the installation process please get root privileges."
print ""
print "wget/curl are needed in Linux/Mac respectively"
print "together with a working Internet connection if packages are installed."
print "--------------------------------------"
print ""

root = False
str = raw_input('Do you have root privileges? [y/n] ').lower()

#cgalFile = 'CGAL-3.9.tar.gz'
#cgalDir = 'CGAL-3.9'
#cgalN = '/29125/'

cgalFile = 'CGAL-4.2-beta1.tar.gz'
cgalDir = 'CGAL-4.2-beta1'
cgalN = '/32183/'

if (str=='y'):
  root = True
else:
  print ""
  print "Without root priviliges you cannot install packets and the compiled libraries"
  print "The needed packets are: boost, gmp, mpfr and cmake."
  print ""
  str = raw_input('Press any key to continue...')
  print ""

compMode = 0
str = raw_input('Do you want to compile as a DelPhi-Module, as Stand-Alone or as Python module? [lib/exe/py] ').lower()

if (str=='lib'):
  print 'Lib mode selected'
  compMode = 1
elif (str=='py'):
  print 'Python mode selected'
  compMode = 2
elif (str=='exe'):
  print 'Stand alone mode selected'
  compMode = 0
else:
  print "Option not recognised, assuming stand-alone executable"

pwd = os.getcwd()
#print "Trying to guess OS ..."

# mac flag
mac = 0
# win flag
win = 0
# snowleopard and leopard force 32 bit compilation (that OS is a bit confused about being 32 or 64 bits...)
mode32 = 0
# enable chrono switch
enable_chrono = False

plat = platform.system()
plat = plat.split('-')
plat = plat[0]
plat = plat.lower()

if (plat == 'darwin' and root==True):
	mac = 1
	print "Detected Mac, installing required packages (boost, gmp, mpfr)"
	print "This script assumes that you have fink software manager installed (it ships with Xcode)\n"
	str = raw_input('Do you want to continue? [y/n] ').lower()
	if (str == 'y'):
		pass
	else:
		quit()
	
	v = platform.mac_ver()[0].split('.')
	
	if (v[0]=='10' and v[1]=='6'):
		print "Detected Snowleopard. For compatibility reasons the build is forced to 32 bits"
		str = raw_input('Do you want to install packets? [y/n] ').lower()
		if (str == 'y'):
			subprocess.call('fink -y install boost1.46.1.cmake',shell=True)
			subprocess.call('fink -y install gmp5',shell=True)
			subprocess.call('fink -y install libmpfr4',shell=True)
			
		mode32 = 1
	
	elif (v[0]=='10' and (v[1]=='7' or v[1]=='8')):
		print "Please perform packets installation by hand before proceeding"
		str = raw_input('Do you want to continue? [y/n] ').lower()
		if (str == 'y'):
			pass
		else:
			quit()
	
	elif (v[0]=='10' and v[1]<='5'):
		print "Detected Leopard or previous OS. For compatibility reasons the build is forced to 32 bits"
		print "Please perform packets installation by hand before proceeding"
		mode32 = 1
		str = raw_input('Do you want to continue? [y/n] ').lower()
		if (str == 'y'):
			pass
		else:
			quit()
	
elif (plat == 'windows'):
	print ''
	print '-----------------------------------------'
	print 'Windows binaries are availale on \\bin'
	print 'To compile on Windows please follow user guide instructions; these are summarized here'
	print '1) Download and compile boost libraries. Set BOOST_DIR to your boost root'
	print '2) Download CGAL. Patch the two include files in CGALPatch dir. Compile CGAL. Set CGAL_DIR to your CGAL root '
	print '3) Cmake of NanoShaper'
	print '4) Open Visual Studio project and compile'
	print '-----------------------------------------'
	quit()
	
elif (root==True):
	distro = platform.linux_distribution()[0]
	distro = distro.lower()
	print "Detected Linux os"

	if (('suse' in distro) or ('opensuse' in distro)):
		print "Detected SUSE Linux distribution, installing required packages: boost, gmp, mpfr, cmake"
		print ""
		str = raw_input('Do you want to install packets? [y/n] ').lower()
		if (str=='y'):
			subprocess.call('yast -i boost boost-devel gmp gmp-devel mpfr mpfr-devel cmake',shell=True)
	
	elif (('ubuntu' in distro) or ('debian' in distro)):
		print "Detected Debian based Linux distribution, installing required packages: boost, gmp, mpfr, cmake, g++" 
		str = raw_input('Do you want to install packets? [y/n] ').lower()
		if (str=='y'):
			subprocess.call('apt-get -y install g++',shell=True)
			subprocess.call('apt-get -y install libboost-thread-dev',shell=True)
			#subprocess.call('apt-get -y install libboost-chrono-dev',shell=True)
			subprocess.call('apt-get -y install libmpfr-dev',shell=True)
			subprocess.call('apt-get -y install cmake',shell=True)
	
	elif (('red' in distro) or ('hat' in distro) or ('redhat' in distro) or ('centos' in distro) or ('fedora' in distro)):
		print "Detected RedHat based Linux distribution, installing required packages: boost, gmp, mpfr, cmake"
		str = raw_input('Do you want to install packets? [y/n] ').lower()
		if (str=='y'):
			subprocess.call('yum install boost',shell=True)
			subprocess.call('yum install boost-devel',shell=True)
			subprocess.call('yum install gmp',shell=True)
			subprocess.call('yum install gmp-devel',shell=True)
			subprocess.call('yum install mpfr',shell=True)
			subprocess.call('yum install mpfr-devel',shell=True)
			subprocess.call('yum install cmake',shell=True)

	else:
		print "I am not able to retrieve a known distro."
		print "You have to install packets manually; these are:"
		print "boost, gmp, mpfr, cmake. If you find boost-devel, gmp-devel etc.."
		print "please also install them to make header files available."
		print ""
		str = raw_input('Do you want to continue? [y/n] ').lower()
		if (str == 'y'):
			pass
		else:
			quit()

# if cgal not found download and compile
if(os.path.isdir('./'+cgalDir)==False):
  print ""
  print "Downloading and extracting CGAL..."
  print ""

  if (mac==1):
    subprocess.call('curl -C - -O https://gforge.inria.fr/frs/download.php'+cgalN+cgalFile,shell=True)
  else:
    subprocess.call('wget https://gforge.inria.fr/frs/download.php'+cgalN+cgalFile,shell=True)

  subprocess.call('sync',shell=True)
  subprocess.call('tar xvf '+cgalFile,shell=True)
  subprocess.call('sync',shell=True)

  # patch the build file: enable chrono and add a index to vertices class
  print ""
  print "Patching CGAL..."
  print ""

  subprocess.call('rm -f ./'+cgalDir+'/include/CGAL/Weighted_point.h',shell=True)
  subprocess.call('sync',shell=True)
  subprocess.call('cp ./CGALPatch/Weighted_point.h ./'+cgalDir+'/include/CGAL/',shell=True)
  subprocess.call('sync',shell=True)

  '''  
  str = raw_input('\nBoost Chrono library will allow to get very accurate execution timings but may be not available in older distros. \nDo you want to compile with Boost Chrono library support? [y/n] ').lower()
  
  flag = False
  if (str == 'y'):
    if (os.path.isfile('./'+cgalDir+'/cmake/modules/CGAL_SetupBoost.cmake.old')):
      print "CGAL already patched, reverting to original file and repatching"
      shutil.copy2('./'+cgalDir+'/cmake/modules/CGAL_SetupBoost.cmake.old','./'+cgalDir+'/cmake/modules/CGAL_SetupBoost.cmake')
      subprocess.call('sync',shell=True)

    shutil.copy2('./'+cgalDir+'/cmake/modules/CGAL_SetupBoost.cmake','./'+cgalDir+'/cmake/modules/CGAL_SetupBoost.cmake.old')
    subprocess.call('sync',shell=True)  
  
    boostDepOld  = open ('./'+cgalDir+'/cmake/modules/CGAL_SetupBoost.cmake.old','r')
    boostDepNew  = open ('./'+cgalDir+'/cmake/modules/CGAL_SetupBoost.cmake','w')
  
    while(1):
    
      line = boostDepOld.readline()

      if (not line):
        boostDepOld.close()
        boostDepNew.close()
        break
    
      v = line.split()
    
      if (len(v)==0):
        continue
    
      if (  line.startswith('  find_package')):
        for j in range(len(v)):
          if (v[j].startswith(')')):
            found = True
            enable_chrono = True
            v[j]=' chrono)\n'
          boostDepNew.write('%s '%v[j])
      else:
        boostDepNew.write(line)
  else:
    found = True
  '''
  found = True
  if (found==True):
    print "\nCGAL correctly patched"
  else:
    print "\n Error in CGAL patching. I need to add 'chrono' dependency into CGAL_SetupBoost.cmake but I cannot find the expected strings"
    quit()

  print ""
  print "Building CGAL..."
  print ""

  if (mode32==0):
   build_cgal = 'cd '+cgalDir+' \n rm -f CMakeCache.txt \n sync \n cmake . -DWITH_examples=false -DWITH_CGAL_Qt4=false -DWITH_CGAL_Qt3=false -DWITH_CGAL_ImageIO=false > cmake_cgal.txt \n sync \n make clean \n make > make_cgal.txt \n sync'
  else:
    # forcing 32 bits mode
    build_cgal = 'cd '+cgalDir+' \n rm -f CMakeCache.txt \n sync \n cmake . -DWITH_examples=false -DWITH_CGAL_Qt4=false -DWITH_CGAL_Qt3=false -DWITH_CGAL_ImageIO=false -DCMAKE_CXX_FLAGS="-arch i386" > cmake_cgal.txt \n sync \n make clean \n make > make_cgal.txt \n sync'

  subprocess.call(build_cgal,shell=True)

print ""

where = ''

if (compMode==0):
  print "Building NanoShaper Stand-Alone please wait..."
  shutil.copy2('CMakeLists_standalone.txt','CMakeLists.txt')
  where = 'build'
elif (compMode==1):
  print "Building NanoShaper DelPhi-lib please wait..."
  shutil.copy2('CMakeLists_lib.txt','CMakeLists.txt')
  where = 'build_lib'
else:
  print "Building NanoShaper Python Module please wait..."
  shutil.copy2('CMakeLists_python.txt','CMakeLists.txt')
  where = 'build_python'

print ""

  
if (mode32==0):    
	build_ns = 'cd '+where+' \n rm -fr * \n cmake .. -DCGAL_DIR=%s/'%(pwd)+cgalDir+' > cmake_ns.txt \n make clean \n make > make_ns.txt'
else:
	# forcing 32 bits mode
	build_ns = 'cd '+where+' \n rm -fr * \n cmake .. -DCGAL_DIR=%s/'%(pwd)+cgalDir+' -DCMAKE_CXX_FLAGS="-arch i386" > cmake_ns.txt \n make clean \n make > make_ns.txt'

subprocess.call(build_ns,shell=True)

if (compMode==0):  
  
  print 'Copying NanoShaper executable in NanoShaper root'
  shutil.copy2('./build/NanoShaper','./')
  
elif (compMode==1 and mac!=1):

  if (root==True):  
    print 'Installing NanoShaper library in /usr/lib or /usr/lib64'
    install_ns = 'cd build_lib \n chmod 777 ./libDelphiSurface.so \n rm -fr /usr/lib64/libDelphiSurface.so \n rm -fr /usr/lib/libDelphiSurface.so \n cp libDelphiSurface.so /usr/lib64 \n cp libDelphiSurface.so /usr/lib \n '
    subprocess.call(install_ns,shell=True)
  else:
    print 'For using the lib with DelPhi please assure that the lib is reachable by the OS'
    print 'Once you have root privileges you can copy it on /usr/lib or /usr/lib64 folder.'
    print 'If you have not root privileges you can update your LD_LIBRARY_PATH to the path where the lib is located'
  
elif (compMode==1 and mac==1):
  
  if (root==True):  
    print 'Installing NanoShaper library in /usr/local/lib'
    install_ns = 'cd build_lib \n chmod 777 ./libDelphiSurface.dylib \n rm -fr /usr/local/lib/libDelphiSurface.dylib \n cp libDelphiSurface.dylib /usr/local/lib \n'
    subprocess.call(install_ns,shell=True)
  else:
    print 'For using the lib with DelPhi please assure that the lib is reachable by the OS'
    print 'Once you have root privileges you can copy it on /usr/lib or /usr/local/lib folder.'
    print 'If you have not root privileges you can update your LD_LIBRARY_PATH to the path where the lib is located'
   

elif (compMode==2):
  
  if (root==True):  
    subprocess.call('cd build_python \n mkdir NanoShaper',shell=True)
    subprocess.call('cd build_python \n cp _NanoShaper.so ./NanoShaper',shell=True)
    subprocess.call('cd build_python \n cp ../setupInstallation.py ./',shell=True)
    subprocess.call('cd build_python \n python setupInstallation.py install > log.txt',shell=True)
    ff = open('./build_python/log.txt','r')
    strs = ff.read().split()
    the_path = ''
    index = 0
    for ss in strs:
      if (ss=='Writing'):
        ss = strs[index+1]
        strs = ss.split('/')
        for ss in strs:
          if (ss.startswith('NanoShaper')==True):
            break 
          else:
            the_path = the_path+ss+'/'
        break
      index = index + 1

    print "Installation was done in ",the_path
    shutil.copy2('./build_python/_NanoShaper.so',the_path)
  else:
    print ""
    print "Without root privileges I cannot perform installation"
    print "Once you get root privileges you can rerun this script or perform installation following these steps:"
    print "1) Go to \\build_python and mkdir NanoShaper"
    print "2) cp _NanoShaper.so ./NanoShaper"
    print "3) cp ../setupInstallation.py ./"
    print "4) python setupInstallation.py install > log.txt"
    print '5) On the log.txt file identify the path that appears after the "Writing" string'
    print '6) Copy _NanoShaper.so into one level before that path. For instance:'
    print '   If in the log you get "/usr/dist-packages/NanoShaper-0.5.egg" then copy into "/usr/dist-packages/"'
    print "   Now you can import NanoShaper in python as per any other package"
    print ""

print ''
print ''