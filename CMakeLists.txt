# Cmake config for NanoShaper
project(NanoShaper) 
cmake_minimum_required(VERSION 3.1)
cmake_policy(VERSION 2.8.4)

if (NOT TBB_DISABLED)
	find_package(TBB REQUIRED)
endif()

include_directories (BEFORE include)
include_directories(src)
include_directories(src/ConfigFile)
include_directories(src/jama)
include_directories(src/tnt)
include_directories(src/sturm)
include_directories(src/octree)

option(ENABLE_PLY "Enable binary PLY mesh output via tinyply" OFF)
set(PLY_SOURCES "")
if (ENABLE_PLY)
	add_definitions(-DPLY_ENABLED)

	set(USE_LOCAL_TINYPLY_FALLBACK OFF)
	if (CMAKE_VERSION VERSION_LESS 3.11)
		message(WARNING "CMake < 3.11: FetchContent unavailable, using local tinyply fallback")
		set(USE_LOCAL_TINYPLY_FALLBACK ON)
	else()
		include(FetchContent)
		FetchContent_Declare(
			tinyply
			GIT_REPOSITORY https://github.com/ddiakopoulos/tinyply.git
			GIT_TAG 40aa4a0ae9e9c203e11893f78b8bcaf8a50e65f0
			GIT_SHALLOW TRUE
		)
		FetchContent_GetProperties(tinyply)
		if (NOT tinyply_POPULATED)
			FetchContent_Populate(tinyply)
		endif()

		if (EXISTS "${tinyply_SOURCE_DIR}/source/tinyply.cpp" AND EXISTS "${tinyply_SOURCE_DIR}/source/tinyply.h")
			include_directories(${tinyply_SOURCE_DIR}/source)
			list(APPEND PLY_SOURCES
				${tinyply_SOURCE_DIR}/source/tinyply.cpp
				${tinyply_SOURCE_DIR}/source/tinyply.h)
			message(STATUS "NanoShaper PLY output enabled (tinyply via FetchContent)")
		else()
			message(WARNING "tinyply fetch incomplete, using local tinyply fallback")
			set(USE_LOCAL_TINYPLY_FALLBACK ON)
		endif()
	endif()

	if (USE_LOCAL_TINYPLY_FALLBACK)
		include_directories(src/tinyply)
		list(APPEND PLY_SOURCES
			./src/tinyply/tinyply.cpp
			./src/tinyply/tinyply.h)
		message(STATUS "NanoShaper PLY output enabled (local fallback)")
	endif()
endif()

# if cgal is found, boost is automatically found
find_package(CGAL)
  
file(GLOB SRC
	"./src/*.h"
	"./src/*.cpp")
	
if (CGAL_FOUND)

	include( ${CGAL_USE_FILE} )		
	add_definitions(-DENABLE_CGAL)

	if(NOT BOOST_DISABLED)
		add_definitions(-DENABLE_BOOST_THREADS)
		message(STATUS "NanoShaper BOOST enabled")
	else()
		message(STATUS "NanoShaper BOOST not enabled (only used in cgal)")
	endif()
	
	find_package( Boost COMPONENTS thread system filesystem REQUIRED ) #whatever libs you need
	INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
		
	add_executable(NanoShaper
								./src/sturm/sturm.cpp
								./src/sturm/util.cpp
								./src/sturm/solve.h
								${SRC}
								./src/ply/ply.cpp
								./src/ply/ply.h								
								${PLY_SOURCES}
	
								./src/ConfigFile/ConfigFile.cpp
								./src/ConfigFile/ConfigFile.h 
								./src/octree/octree.h 
								)

	if(NOT BOOST_DISABLED)
		TARGET_LINK_LIBRARIES(NanoShaper PUBLIC ${Boost_LIBRARIES})
	endif(NOT BOOST_DISABLED)

	if (NOT TBB_DISABLED)
		include(CGAL_TBB_support)
		if(TARGET CGAL::TBB_support)
			target_link_libraries(NanoShaper PUBLIC CGAL::TBB_support)
			message(STATUS "TBB found")
		else()
			message(STATUS "TBB not found")
		endif()
	else()
		message(STATUS "TBB not enabled")
	endif()
			
else(CGAL_FOUND)
		
	message(STATUS "CGAL not found, compilation stopped")  
	
endif(CGAL_FOUND)
