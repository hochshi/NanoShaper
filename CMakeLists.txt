
# Cmake config for DelphiSurface module
project(DelphiSurface) 
cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

include_directories (BEFORE include)
include_directories(src)
include_directories(src/ConfigFile)
include_directories(src/jama)
include_directories(src/tnt)
include_directories(src/sturm)
include_directories(src/octree)

find_package(Boost COMPONENTS thread)
if (Boost_FOUND)
  link_directories ( ${Boost_LIBRARY_DIRS} )
  link_libraries(${Boost_LIBRARIES})
  include_directories ( ${Boost_INCLUDE_DIRS} )
  include_directories (BEFORE include)
  add_definitions(-DENABLE_BOOST_THREADS)
else (Boost_FOUND)
  message(STATUS "Boost libs not found, Multithreading will be disabled")  
endif (Boost_FOUND)
	
# if cgal is found, boost is automatically found
 find_package(CGAL)
 
 file(GLOB SRC
    "./src/*.h"
    "./src/*.cpp")

 if (CGAL_FOUND)

	include( ${CGAL_USE_FILE} )		
	add_definitions(-DENABLE_CGAL)
	add_definitions(-DENABLE_BOOST_THREADS)
	
	if (${UNIX}) 
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mtune=native") 
	endif (${UNIX}) 
	
	#######################
	### delphi bind mode ##
	#######################
	add_definitions(-DDELPHI_BIND)
	#delphi module  
	if (MSVC)
	add_library(libDelphiSurface SHARED
								./src/sturm/sturm.cpp
								./src/sturm/util.cpp
								./src/sturm/solve.h
								${SRC}
								./src/ply/ply.cpp
								./src/ply/ply.h
								./src/ConfigFile/ConfigFile.cpp
								./src/ConfigFile/ConfigFile.h 
								./src/octree/octree.h 
								)    
	set_target_properties(libDelphiSurface PROPERTIES COMPILE_FLAGS "/wd4503")
	else (MSVC)
	add_library(DelphiSurface SHARED
								./src/sturm/sturm.cpp
								./src/sturm/util.cpp
								./src/sturm/solve.h
								${SRC}
								./src/ply/ply.cpp
								./src/ply/ply.h								
								./src/ConfigFile/ConfigFile.cpp
								./src/ConfigFile/ConfigFile.h 
								./src/octree/octree.h 
								)    
	endif (MSVC)

 else(CGAL_FOUND)
		
	message(STATUS "CGAL not found, Skin Surface and Connolly Surface will be disabled")
	
	if (${UNIX}) 
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive") 
	endif (${UNIX}) 

	add_definitions(-DDELPHI_BIND)
	
	if (MSVC)
		add_definitions(-D_CRT_SECURE_NO_WARNINGS)	
	endif (MSVC)  
	
	if (CMAKE_COMPILER_IS_GNUCXX)		
		#...
	endif (CMAKE_COMPILER_IS_GNUCXX)
	
	if(CMAKE_C_COMPILER MATCHES "icc")
		#..
	endif ()
	
	message(STATUS "C++ Release flags ${CMAKE_CXX_FLAGS_RELEASE}")
	
	#######################
	### delphi bind mode ##
	#######################
	add_definitions(-DDELPHI_BIND)	
	#delphi module  
	add_library(DelphiSurface SHARED 
								./src/sturm/sturm.cpp
								./src/sturm/util.cpp
								./src/sturm/solve.h
								${SRC}
								./src/ConfigFile/ConfigFile.cpp
								./src/ConfigFile/ConfigFile.h 
								./src/octree/octree.h 
								)    
								
 endif(CGAL_FOUND)