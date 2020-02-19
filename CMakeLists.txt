cmake_minimum_required(VERSION 3.15)
project(tbkp C)

set(CMAKE_C_STANDARD 11)

add_library(combo src/utils/combo.h src/utils/combo.c)
add_executable(tbkp src/main.c src/tbkp_instance.h src/tbkp_instance.c src/tbkp_ub.h src/tbkp_ub.c)

set(WARNING_OPTIONS             -Wall -Wextra -Werror)
set(GCC_WARNING_OPTIONS         -Wcast-align -Wconversion -Wmisleading-indentation
                                -Wduplicated-branches -Wduplicated-cond -Wlogical-op
                                -Wnull-dereference -Wdouble-promotion -Wformat=2)
set(ARCH_OPTIONS                -march=native)
set(RELEASE_OPTIONS             -Ofast)
set(DEBUG_OPTIONS               -O0 -g)
set(GDB_DEBUG_OPTIONS           -ggdb)

target_compile_options(     combo PRIVATE ${ARCH_OPTIONS})
target_compile_options(     combo PRIVATE $<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>)
target_compile_options(     combo PRIVATE $<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>)
target_compile_options(     combo PRIVATE $<$<AND:$<CONFIG:DEBUG>,$<CXX_COMPILER_ID:GNU>>:${GDB_DEBUG_OPTIONS}>)

target_include_directories( tbkp SYSTEM PRIVATE "src/utils")
target_compile_options(     tbkp PRIVATE ${WARNING_OPTIONS} ${ARCH_OPTIONS})
target_compile_options(     tbkp PRIVATE $<$<CXX_COMPILER_ID:GNU>:${GCC_WARNING_OPTIONS}>)
target_compile_options(     tbkp PRIVATE $<$<CONFIG:RELEASE>:${RELEASE_OPTIONS}>)
target_compile_options(     tbkp PRIVATE $<$<CONFIG:DEBUG>:${DEBUG_OPTIONS}>)
target_compile_options(     tbkp PRIVATE $<$<AND:$<CONFIG:DEBUG>,$<CXX_COMPILER_ID:GNU>>:${GDB_DEBUG_OPTIONS}>)
target_link_libraries(      tbkp PRIVATE combo)
