file(GLOB GRAPHICS_SOURCES "graphics/*.cpp")
file(GLOB PRIMARY_SOURCES "config/*.cpp" "detectors/*.cpp" "fields/*.cpp" "integrators/*.cpp" "interpolation/*.cpp" "particles/*.cpp" "sources/*.cpp" "util/*.cpp")

include_directories(/usr/lib/x86_64-linux-gnu/hdf5/serial/include)
link_directories(/usr/lib/x86_64-linux-gnu/hdf5/serial/lib)

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
