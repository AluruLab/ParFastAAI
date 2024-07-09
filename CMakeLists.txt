cmake_minimum_required(VERSION 3.5)
# Project Name
set(PROJECT_NAME_STR par_fastaai)
set(PARFASTAAI_VERSION 1.0)
project(${PROJECT_NAME_STR} C CXX)

#
set(CMAKE_VERBOSE_MAKEFILE ON)
set(BUILD_SHARED_LIBS ON)

# option to build test
option(BUILD_TESTING "Build tests" ON) # build tests

set(SQLITE_SRC_DIR "${PROJECT_SOURCE_DIR}/ext/sqlite")
set(FMT_SRC_DIR "${PROJECT_SOURCE_DIR}/ext/fmt")
set(CATCH2_SRC_DIR "${PROJECT_SOURCE_DIR}/ext/Catch2")
set(CEREAL_SRC_DIR "${PROJECT_SOURCE_DIR}/ext/cereal")
# set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_BUILD_TYPE Release)

# include directories
include_directories("${PROJECT_SOURCE_DIR}/include/")

# Additional cmake files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
message(STATUS "MODULE PATH   ${CMAKE_MODULE_PATH}")

# External headers for sqlite
add_library(sqlite INTERFACE)
target_include_directories(sqlite INTERFACE ${SQLITE_SRC_DIR})
add_library(fmt INTERFACE)
target_include_directories(fmt INTERFACE ${FMT_SRC_DIR}/include)
add_library(cereal INTERFACE)
target_include_directories(cereal INTERFACE ${CEREAL_SRC_DIR}/include)

# Check if OpenMP is available
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP version:  ${OpenMP_CXX_VERSION}")
else()
  message(STATUS "OpenMP NOT FOUND!")
endif()

# Compile
set(PROJECT_EXEC ${PROJECT_NAME_STR}.x)
set(PROJECT_TEST_EXEC ${PROJECT_NAME_STR}_tests.x)

# Executable
add_executable(${PROJECT_EXEC} src/main.cpp ${SQLITE_SRC_DIR}/sqlite3.c
    ${FMT_SRC_DIR}/src/format.cc)
#
add_compile_definitions(PARFASTAAI_VERSION=${PARFASTAAI_VERSION})
# OpenMP
include_directories(${OpenMP_CXX_INCLUDE_DIRS})
target_link_libraries(${PROJECT_EXEC} PUBLIC OpenMP::OpenMP_CXX sqlite fmt
                      cereal ${CMAKE_DL_LIBS})


if(${BUILD_TESTING})
    add_subdirectory(ext/Catch2)
    # 
    add_executable(${PROJECT_TEST_EXEC} tests/pfaai_tests.cpp 
        ${SQLITE_SRC_DIR}/sqlite3.c ${FMT_SRC_DIR}/src/format.cc)
    target_include_directories(${PROJECT_TEST_EXEC} PUBLIC 
                              ${OpenMP_CXX_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_TEST_EXEC} PRIVATE OpenMP::OpenMP_CXX sqlite
                                                       fmt cereal Catch2::Catch2WithMain)
    file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif(${BUILD_TESTING})
