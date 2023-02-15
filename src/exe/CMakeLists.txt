add_executable(NanoShaper
  ./main.cpp
  ${SOURCE_FILES}
  ../sturm/sturm.cpp
  ../sturm/util.cpp
  ../sturm/solve.h
  ../ply/ply.cpp
  ../ply/ply.h								
  ../ConfigFile/ConfigFile.cpp
  ../ConfigFile/ConfigFile.h 
  ../octree/octree.h 
  )

set_target_properties(NanoShaper
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
