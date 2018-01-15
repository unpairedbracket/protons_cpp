if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(/usr/lib/x86_64-linux-gnu/hdf5/serial/include)
link_directories(/usr/lib/x86_64-linux-gnu/hdf5/serial/lib)

file(GLOB GRAPHICS_SOURCES "graphics/*.cpp")
file(GLOB PRIMARY_SOURCES "config/*.cpp" "detectors/*.cpp" "fields/*.cpp" "integrators/*.cpp" "interpolation/*.cpp" "particles/*.cpp" "sources/*.cpp" "util/*.cpp")

add_library(Code ${PRIMARY_SOURCES})
target_link_libraries(Code -ldl -lz -lm -lyaml-cpp -lCGAL -lCGAL_Core -lgmp -lpthread)
target_link_libraries(Code libhdf5_hl_cpp.so libhdf5_hl.so libhdf5_cpp.so libhdf5.so)

add_library(Graphics ${GRAPHICS_SOURCES})
target_link_libraries(Graphics -lGL -lX11 -lGLU -lGLEW -lglfw -lpng)

add_executable(protons main.cpp)
target_compile_definitions(protons PUBLIC USE_GL)
target_link_libraries(protons Code)
target_link_libraries(protons Graphics)

add_executable(protons-headless main.cpp)

target_link_libraries(protons-headless Code)

install(TARGETS protons protons-headless DESTINATION bin)
