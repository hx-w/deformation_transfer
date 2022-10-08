add_subdirectory(meshlib)
list(APPEND MESH_LIB)
set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS "-pthread" LINK_FLAGS "-pthread")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

add_definitions(-D_CRT_SECURE_NO_WARNINGS)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} MESH_LIB) 
target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_BINARY_DIR}")