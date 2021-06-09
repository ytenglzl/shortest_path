cmake_minimum_required(VERSION 2.8)
project(shortest_path)

## Project specific configurations go here

# enable c++11 support
set (CMAKE_CXX_FLAGS "-std=c++11 -Wall ${CMAKE_CXX_FLAGS}")
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/run)

# create the main executable
## add additional .cpp files if needed
add_executable(path_gener path_gener.cpp )
add_executable(constructor constructor.cpp )
add_executable(rgen rgen.cpp )

# create the executable for tests
add_executable(tests test.cpp)
install(PROGRAMS constructor.py DESTINATION bin)
install(TARGETS rgen path_gener master RUNTIME DESTINATION bin)
