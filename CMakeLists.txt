cmake_minimum_required (VERSION 3.8)

project(pipefoil)
add_executable(pipefoil main.cc)

include(CheckLibraryExists)
CHECK_LIBRARY_EXISTS(m sin "" HAVE_LIB_M)                                                                                                
                                                                                                                                         
if (HAVE_LIB_M)
    target_link_libraries(pipefoil m)
endif (HAVE_LIB_M)

if (MSVC)
    target_compile_options(pipefoil PUBLIC /GR- "/Zc:__cplusplus")
endif()

if (CMAKE_COMPILER_IS_GNUCC)
    target_compile_options(pipefoil PUBLIC -fno-rtti -finput-charset=UTF-8 -Wno-normalized -Wno-multichar -Wno-literal-suffix)
endif()
#set(CMAKE_CXX_EXTENSIONS OFF)      # no vendor-specific extensions
#set_property(TARGET pipefoil PROPERTY CXX_STANDARD 20)
