# add_subdirectory(meshlib)
set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS "-pthread" LINK_FLAGS "-pthread")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC ${LIBS}) 
target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_BINARY_DIR}")