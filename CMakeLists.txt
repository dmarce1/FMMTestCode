
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -march=native")

add_executable( fmm_test
	src/main.cpp
)

