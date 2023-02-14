# Cmake config for NanoShaper
project(NanoShaper LANGUAGES CXX C)
cmake_minimum_required(VERSION 3.6)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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
# Git submodules
################################################################################
option(ENABLE_GIT_SUBMODULES "init and update git submodules" OFF)

if(ENABLE_GIT_SUBMODULES)
  include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists_git.txt)
endif()


################################################################################
# SPDLog
################################################################################
option(ENABLE_SPDLOG "SPDLog" OFF)

if(ENABLE_SPDLOG)
  add_definitions(-DSPDLOG_ENABLED)
  link_libraries(
    spdlog::spdlog_header_only
  )
  if(ENABLE_GIT_SUBMODULES)
    add_subdirectory(lib/spdlog)
  endif()
endif()

################################################################################
# NLOHMANN_JSON
################################################################################
option(ENABLE_JSON "JSON" OFF)

if(ENABLE_JSON)
  set(JSON_BuildTests OFF CACHE INTERNAL "")
  add_definitions(-DJSON_ENABLED)
  link_libraries(
    nlohmann_json::nlohmann_json
  )
  if(ENABLE_GIT_SUBMODULES)
    add_subdirectory(lib/json)
  endif()
endif()

# if (${UNIX}) 
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mtune=native") 
# endif (${UNIX}) 

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
