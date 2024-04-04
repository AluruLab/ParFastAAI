cmake_minimum_required(VERSION 3.0)
# Project Name
set(PROJECT_NAME_STR par_fastaai)
set(PARFASTAAI_VERSION 1.0)
project(${PROJECT_NAME_STR} C CXX)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(BUILD_SHARED_LIBS ON)
set(SQLITE_SRC_DIR "${PROJECT_SOURCE_DIR}/ext/sqlite-amalgamation-3450100")
# set(CMAKE_BUILD_TYPE Debug) set(CMAKE_BUILD_TYPE Release)

# include directories
include_directories("${PROJECT_SOURCE_DIR}/include/")

# Additional cmake files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
message(STATUS "MODULE PATH   ${CMAKE_MODULE_PATH}")

# External headers for sqlite
add_library(sqlite INTERFACE)
target_include_directories(sqlite INTERFACE ${SQLITE_SRC_DIR})

# Check if OpenMP is available
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP version:  ${OpenMP_CXX_VERSION}")
else()
  message(STATUS "OpenMP NOT FOUND!")
endif()

# Compile
set(PROJECT_EXEC ${PROJECT_NAME_STR}.x)

# Executable
add_executable(${PROJECT_EXEC} src/main.cpp ${SQLITE_SRC_DIR}/sqlite3.c)
#
add_compile_definitions(PARFASTAAI_VERSION=${PARFASTAAI_VERSION})
# OpenMP
include_directories(${OpenMP_CXX_INCLUDE_DIRS})
target_link_libraries(${PROJECT_EXEC} PUBLIC OpenMP::OpenMP_CXX sqlite
                                             ${CMAKE_DL_LIBS})
