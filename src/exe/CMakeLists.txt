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
if(ENABLE_SPDLOG)
  target_link_libraries(NanoShaper
    PRIVATE spdlog::spdlog_header_only
  )
endif()
if(ENABLE_NLOHMANN_JSON)
  target_link_libraries(NanoShaper
    PRIVATE nlohmann_json::nlohmann_json
  )
endif()
