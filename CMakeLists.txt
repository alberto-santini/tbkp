cmake_minimum_required(VERSION 3.1)
project(tbkp)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

add_library(pdqsort src/utils/pdqsort.h src/utils/pdqsort_c.cpp)
add_library(combo src/utils/combo.h src/utils/combo.c src/utils/pdqsort_c.cpp src/utils/pdqsort_c.h)
add_library(argparse src/utils/argparse.h src/utils/argparse.c)
add_executable(tbkp
        src/main.c
        src/tbkp_instance.h
        src/tbkp_instance.c
        src/tbkp_de_sol.h
        src/tbkp_de_sol.c
        src/tbkp_boole_sol.h
        src/tbkp_boole_sol.c
        src/tbkp_bb.c
        src/tbkp_bb.h src/tbkp_bb_stats.h src/tbkp_bb_stats.c src/tbkp_solution.h src/tbkp_solution.c src/tbkp_bb_fixed_status.h src/tbkp_bb_defs.h)

set(WARNING_OPTIONS             -Wall -Wextra -Werror)
set(GCC_WARNING_OPTIONS         -Wcast-align -Wconversion -Wmisleading-indentation
                                -Wduplicated-branches -Wduplicated-cond -Wlogical-op
                                -Wnull-dereference -Wformat=2 -Wno-missing-field-initializers)
set(ARCH_OPTIONS                -march=native)
set(RELEASE_OPTIONS             -Ofast)
set(DEBUG_OPTIONS               -O0 -g)
set(GDB_DEBUG_OPTIONS           -ggdb)

find_package(GUROBI REQUIRED)

target_compile_options(     pdqsort PRIVATE ${ARCH_OPTIONS})
target_compile_options(     pdqsort PRIVATE $<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>)
target_compile_options(     pdqsort PRIVATE $<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>)
target_compile_options(     pdqsort PRIVATE $<$<AND:$<CONFIG:DEBUG>,$<CXX_COMPILER_ID:GNU>>:${GDB_DEBUG_OPTIONS}>)

target_compile_options(     combo PRIVATE ${ARCH_OPTIONS})
target_compile_options(     combo PRIVATE $<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>)
target_compile_options(     combo PRIVATE $<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>)
target_compile_options(     combo PRIVATE $<$<AND:$<CONFIG:DEBUG>,$<CXX_COMPILER_ID:GNU>>:${GDB_DEBUG_OPTIONS}>)

target_compile_options(     argparse PRIVATE ${ARCH_OPTIONS})
target_compile_options(     argparse PRIVATE $<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>)
target_compile_options(     argparse PRIVATE $<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>)
target_compile_options(     argparse PRIVATE $<$<AND:$<CONFIG:DEBUG>,$<CXX_COMPILER_ID:GNU>>:${GDB_DEBUG_OPTIONS}>)

target_include_directories( tbkp SYSTEM PRIVATE "src/utils")
target_include_directories( tbkp SYSTEM PRIVATE ${GUROBI_INCLUDE_DIRS})
target_compile_options(     tbkp PRIVATE ${WARNING_OPTIONS} ${ARCH_OPTIONS})
target_compile_options(     tbkp PRIVATE $<$<CXX_COMPILER_ID:GNU>:${GCC_WARNING_OPTIONS}>)
target_compile_options(     tbkp PRIVATE $<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>)
target_compile_options(     tbkp PRIVATE $<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>)
target_compile_options(     tbkp PRIVATE $<$<AND:$<CONFIG:DEBUG>,$<CXX_COMPILER_ID:GNU>>:${GDB_DEBUG_OPTIONS}>)
target_compile_definitions( tbkp PRIVATE $<$<CONFIG:RELEASE>:BB_VERBOSITY_CURRENT=0>)
target_link_libraries(      tbkp PRIVATE combo)
target_link_libraries(      tbkp PRIVATE pdqsort)
target_link_libraries(      tbkp PRIVATE argparse)
target_link_libraries(      tbkp PRIVATE ${GUROBI_LIBRARY})
