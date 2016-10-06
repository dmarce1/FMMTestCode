
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAG} -march=native -std=c++11 -O3")

add_executable( fmm_test
	src/main.cpp
)

