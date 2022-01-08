# CMakeList.txt : CMake project for RayTracer, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project(RayTracer)

set( DCMAKE_TOOLCHAIN_FILE "D:/Code/vcpkg/Microsoft/vcpkg/scripts/buildsystems/vcpkg.cmake")

find_package(OpenGL)
find_package(GLUT)

set( CXX_SRCS
	main.cpp
	RayTracer.cpp
)

include_directories(
	"D:/Code/vcpkg/Microsoft/vcpkg/installed/x64-windows/include/GL"
)


# Add source to this project's executable.
add_executable (RayTracer ${CXX_SRCS})

# TODO: Add tests and install targets if needed.
target_link_libraries( RayTracer
	OpenGL32
	GLUT::GLUT
)