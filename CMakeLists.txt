# CMakeList.txt : CMake project for CS715, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project (CS715)

# Add source to this project's executable.
add_executable (CS715 "CS715.cpp" "CS715.h" "VectorField.h" "VectorField.cpp" "Vector.cpp" "Vector.h")

find_library(FREEGLUT_LIB freeglut "C:\\Program Files\\Common Files\\freeglut\\lib\\x64")

target_include_directories(CS715 PUBLIC "C:\\Program Files\\Common Files\\freeglut\\include")
target_link_libraries(CS715 ${FREEGLUT_LIB})

# TODO: Add tests and install targets if needed.
