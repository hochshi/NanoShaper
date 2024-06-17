# Cmake config for NanoShaper
cmake_minimum_required(VERSION 3.6)
project(NanoShaper LANGUAGES CXX C)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

cmake_policy(SET CMP0042 NEW) # MacOS RPATH on by default
cmake_policy(SET CMP0048 NEW) # project() command manages the VERSION variables
cmake_policy(SET CMP0054 NEW) # Only interpret if() arguments as variables or keywords when unquoted
cmake_policy(SET CMP0077 NEW) # option() honors normal variables (i.e. does nothing if a normal variable with the same name exists)
cmake_policy(SET CMP0083 NEW) # Pass flags needed for position-independent executables
cmake_policy(SET CMP0091 NEW) # MSVC runtime library flags are selected by an abstraction (i.e. CMAKE_MSVC_RUNTIME_LIBRARY)
set(CMAKE_MACOSX_RPATH 1)

include(FetchContent)
include(ExternalProject)

################################################################################
# Build type
################################################################################
set(default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")

  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)

  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "RelWithDebInfo")
endif ()

################################################################################
# SPDLog
################################################################################
option(ENABLE_SPDLOG "SPDLog" OFF)

if(ENABLE_SPDLOG)
  FetchContent_Declare( spdlog
      GIT_REPOSITORY https://github.com/gabime/spdlog.git
      GIT_TAG 27cb4c76708608465c413f6d0e6b8d99a4d84302 # v1.14.1
  )
  FetchContent_MakeAvailable( spdlog )
  add_definitions(-DSPDLOG_ENABLED)
  link_libraries(
    spdlog::spdlog_header_only
  )
endif()

################################################################################
# NLOHMANN_JSON
################################################################################
option(ENABLE_JSON "JSON" OFF)

if(ENABLE_JSON)
  FetchContent_Declare( json
      GIT_REPOSITORY https://github.com/nlohmann/json.git
      GIT_TAG bc889afb4c5bf1c0d8ee29ef35eaaf4c8bef8a5d # v3.11.2
  )
  FetchContent_MakeAvailable( json )
  set(JSON_BuildTests OFF CACHE INTERNAL "")
  add_definitions(-DJSON_ENABLED)
  link_libraries(
    nlohmann_json::nlohmann_json
  )
endif()

################################################################################
# TinyPLY
################################################################################
option(ENABLE_PLY "PLY" OFF)

if(ENABLE_PLY)

  FetchContent_Declare( tinyply
      GIT_REPOSITORY https://github.com/ddiakopoulos/tinyply.git
      GIT_TAG 40aa4a0ae9e9c203e11893f78b8bcaf8a50e65f0 # 2.3.4
  )

  FetchContent_GetProperties(tinyply)
  if(NOT tinyply_POPULATED)
    FetchContent_Populate(tinyply)
    add_library(tinyply STATIC
      ${tinyply_SOURCE_DIR}/source/tinyply.cpp
      ${tinyply_SOURCE_DIR}/source/tinyply.h
    )
  endif()
  
  add_definitions(-DPLY_ENABLED)
  link_libraries(
    tinyply
  )
endif()

################################################################################
# BOOST
################################################################################
option(ENABLE_BOOST "BOOST" OFF)

if(ENABLE_BOOST)
  find_package(Boost REQUIRED COMPONENTS thread)
endif()

if (Boost_FOUND)
  link_directories ( ${Boost_LIBRARY_DIRS} )
  link_libraries(${Boost_LIBRARIES})
  include_directories ( ${Boost_INCLUDE_DIRS} )
  include_directories (BEFORE include)
  add_definitions(-DENABLE_BOOST_THREADS)
else (Boost_FOUND)
  message(STATUS "Boost libs not found, Multithreading will be disabled")  
endif (Boost_FOUND)


################################################################################
# CGAL
################################################################################
option(ENABLE_CGAL "CGAL" OFF)

if(ENABLE_CGAL)
  find_package(CGAL)
endif()

if (CGAL_FOUND)

  set(CMAKE_CXX_STANDARD 14)
  set(CXX_STANDARD_REQUIRED ON)
	include( ${CGAL_USE_FILE} )		
	add_definitions(-DENABLE_CGAL)
	
	if (MSVC)
		set_target_properties(NanoShaper PROPERTIES COMPILE_FLAGS "/wd4503")	
		#message(STATUS "Adding optimizations for processor ${CMAKE_SYSTEM_PROCESSOR}")
	endif (MSVC)  

else(CGAL_FOUND)
		
	message(STATUS "CGAL not found, Skin Surface and Connolly Surface will be disabled")  
	
	if (MSVC)
		add_definitions(-D_CRT_SECURE_NO_WARNINGS)	
	endif (MSVC)  
	
	if (CMAKE_COMPILER_IS_GNUCXX)		
		#...
	endif (CMAKE_COMPILER_IS_GNUCXX)
	
	if(CMAKE_C_COMPILER MATCHES "icc")
		#..
	endif ()
	
endif(CGAL_FOUND)

################################################################################
# Setting up src builds
################################################################################
add_subdirectory(src)
