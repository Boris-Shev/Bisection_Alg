project(BisectionAlgorithm CXX)

cmake_minimum_required(VERSION 3.2)

set(BISEC_INCLUDE_DIRS "include")

if (NOT MSVC)
  set(SOME_LIBS "pthread")
endif ()

include_directories(${BISEC_INCLUDE_DIRS})

add_compile_options(
      -o3
      -mfpmath=sse
      -fstack-protector-all
      -g
      -W
      -Wall
      -Wextra
      -Wunused
      -Wcast-align
      -Werror -pedantic
      -pedantic-errors
      -Wfloat-equal
      -Wpointer-arith
      -Wformat-security
      -Wmissing-format-attribute
      -Wformat=1
      -Wwrite-strings
      -Wcast-align
      -Wno-long-long
      -Woverloaded-virtual
      -Wnon-virtual-dtor
      -Wcast-qual
      -Wno-suggest-attribute=format
)

add_executable(
  bisec_alg
  include/header.hpp
  source/main.cpp
  source/matrx_funcs.cpp
  source/bisection_alg.cpp
)

target_link_libraries(bisec_alg ${SOME_LIBS})
