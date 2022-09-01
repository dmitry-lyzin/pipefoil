cmake_minimum_required(VERSION 2.8.7)

project(pipefoil)
add_executable(pipefoil main.cc)

include(CheckLibraryExists)
CHECK_LIBRARY_EXISTS(m sin "" HAVE_LIB_M)                                                                                                
                                                                                                                                         
if (HAVE_LIB_M)
    target_link_libraries(pipefoil m)
endif (HAVE_LIB_M)
