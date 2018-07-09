cmake_minimum_required(VERSION 3.9)
project(mcf)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -O3")

set(source_files
        defines.h
        implicit.c
        implicit.h
        mcf.c
        mcf.h
        mcflimit.h
        mcfutil.c
        mcfutil.h
        output.c
        output.h
        pbeampp.c
        pbeampp.h
        pbla.c
        pbla.h
        pflowup.c
        pflowup.h
        prototyp.h
        psimplex.c
        psimplex.h
        pstart.c
        pstart.h
        readmin.c
        readmin.h
        treeup.c
        treeup.h)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)
add_executable(mcf
        ${source_files})
