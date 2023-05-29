cmake_minimum_required(VERSION 3.19)
project(final_compress)

set(CMAKE_CXX_STANDARD 14)

add_executable(compress src/basic.h src/basic.cpp src/Binary_Search.cpp src/Binary_Search.h src/BplusTreeBit.h src/BplusTreeBit.cpp src/compressing.cpp src/compressing.h
        src/exactMatchFMindex.cpp src/exactMatchFMindex.h src/Hash.h src/Hash.cpp src/inputFMindex.cpp src/inputFMindex.h src/inputRef.cpp src/inputRef.h src/load_DBG_full.h
        src/load_DBG_full.cpp src/Main.cpp src/utils.cpp src/utils.h src/split.cpp src/split.h)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")