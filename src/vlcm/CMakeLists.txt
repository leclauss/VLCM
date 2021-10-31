cmake_minimum_required(VERSION 3.16)
include(CheckCXXCompilerFlag)
# ----------------------------------------------------------------------------------------
#                              ===== Project Setup =====
project(VLCM)

# Release build by default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ----------------------------------------------------------------------------------------
#                         ===== Compiler Configuration =====

set(CMAKE_CXX_STANDARD 17)

CHECK_CXX_COMPILER_FLAG("-fPIC" COMPILER_OPT_PIC_SUPPORTED)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
CHECK_CXX_COMPILER_FLAG("-O3" COMPILER_OPT_O3_SUPPORTED)
CHECK_CXX_COMPILER_FLAG("-funroll-loops" COMPILER_OPT_UNROLL_LOOPS_SUPPORTED)
CHECK_CXX_COMPILER_FLAG("-ffp-contract=fast" COMPILER_OPT_FPCONTRACT_FAST_SUPPORTED)
CHECK_CXX_COMPILER_FLAG("-Wall" COMPILER_OPT_WARN_ALL_SUPPORTED)

if (COMPILER_OPT_PIC_SUPPORTED)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fPIC")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fPIC")
endif()

if (COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
endif()

if (COMPILER_OPT_O3_SUPPORTED)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif()

if (COMPILER_OPT_UNROLL_LOOPS_SUPPORTED)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funroll-loops")
endif()

if (COMPILER_OPT_FPCONTRACT_FAST_SUPPORTED)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffp-contract=fast")
endif()

if (COMPILER_OPT_WARN_ALL_SUPPORTED)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
endif()

CHECK_CXX_COMPILER_FLAG("-fsanitize=address" COMPILER_OPT_SANITIZE_ADDRESS_SUPPORTED)
CHECK_CXX_COMPILER_FLAG("-fno-omit-frame-pointer" COMPILER_OPT_NO_OMIT_FP_SUPPORTED)

if (COMPILER_OPT_SANITIZE_ADDRESS_SUPPORTED)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
endif()

if (COMPILER_OPT_NO_OMIT_FP_SUPPORTED)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")
endif()

# ----------------------------------------------------------------------------------------
#                              ===== Build targets =====

set(SOURCES ${SOURCES} utils.cpp utils.h exception.h precompute.cpp precompute.h lmc.cpp lmc.h graph_utils.h distance_graph.cpp distance_graph.h graph_utils.cpp fast_cliques.cpp fast_cliques.h)

add_executable(VLCM main.cpp ${SOURCES})