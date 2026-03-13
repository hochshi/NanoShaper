# NanoShaper 1.5

NanoShaper is a tool aimed at computing the molecular surface and pockets of a biomolecular system. It does that through a ray tracing
technique which renders NanoShaper compatible with patch-based, implicit surface or triangulated surfaces.
NanoShaper is GPL licensed software. The 1.5 version parallelizes completely the SES build-up and bears many optimizations which support 
increased accuracy, improved speed, and a significantly reduced memory footprint. 


1. Contents

	A. Files and folders to compile NanoShaper
	- "\src" C++ sources
	- "\src_client" C++ sources and example main for API usage
	- "setup.py" script to support the compilation

	B. Example input files for standalone NanoShaper
	- "\example" 

	C. Docker file and image
	- "Dockerfile" 	

	

2. Installation

	NanoShaper can be compiled on Linux for x86-64 with gcc, tested from 8 to 9.5. 
	Pre-requisites libraries are: boost, gmp, mpfr.
	To install please run the setup.py script and select if a standalone, or .so (for API usage) 
	is required.
	
	```bash
	git clone https://gitlab.iit.it/sdecherchi/nanoshaper
 	git checkout NanoShaper_PatchBased
	python setup.py
	```	


3. Execution	
	
	Running the NanoShaper executable file if you installed as described above
	```console
	cd <Path of the directory where input files are>
	NanoShaper  <Surface Configuration file>
	```
	To use NanoShaper in VMD please add NanoShaper path to your $PATH variable.      
	In order to use nanoshaper docker image inside VMD, please
	rename NanoShaper_Docker_VMD.sh into NanoShaper and add its path to your $PATH environment variable.
	Also please be sure that in the .sh file the correct Docker image in your system is referenced.
	In some OS X systems (e.g., Big Sur) you may need to change the VMD default temporary folder using the VMDTMPDIR env variable because of permission issues.		