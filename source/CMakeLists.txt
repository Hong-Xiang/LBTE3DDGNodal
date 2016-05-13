cmake_minimum_required (VERSION 3.0)
project (LBTE3DDGNodal)

#========Add Libraries:=========================#

get_filename_component(CMAKE_MY_MODULE $ENV{USERPROFILE}/OneDrive/Template/cmake_modules ABSOLUTE)
get_filename_component(XLIB_PATH $ENV{USERPROFILE}/OneDrive/xlib ABSOLUTE)
set(CMAKE_MODULE_PATH "${CMAKE_MY_MODULE}")
set(LIBS)
find_package(MKL)
set(LIBS ${MKL_LIBRARIES} ${LIBS})
include_directories(${MKL_INCLUDE_DIR})

find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY MAT_LIBRARY ENG_LIBRARY)
set(LIBS ${Matlab_LIBRARIES} ${LIBS})
include_directories(${Matlab_INCLUDE_DIRS})

find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})
set(LIBS ${Boost_LIBRARIES} ${LIBS})
link_directories(${Boost_LIBRARY_DIRS})

include_directories(${XLIB_PATH}//include)
find_library(XLib_LIB xlib ${XLIB_PATH}//lib)
set(LIBS ${XLib_LIB})

#================================================#

#========Compile options:========================#
                                                 
#add_compile_options(/DMKL_ILP64)
#add_compile_options(/DMKL_INT=int)
add_compile_options(/Qopenmp)
add_compile_options(/DUSE_DOUBLE_PRECISION)
#================================================#

#========Add source files:=======================#

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(SOURCEFILES)
set(HEADERFILES)
add_subdirectory(include)
add_subdirectory(source)
message(STATUS "Headers: ${HEADERFILES}")
message(STATUS "Sources: ${SOURCEFILES}")


#================================================#

#========Build and link:=========================#
add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp ${SOURCEFILES} ${HEADERFILES})
target_link_libraries(${PROJECT_NAME} ${LIBS})

#================================================#