cmake_minimum_required(VERSION 3.9)
project(mcf)

#set(CMAKE_CXX_STANDARD 11)

find_package(OpenMP REQUIRED)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -pg")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

set(C_FLAGS "-fopenmp -O3")
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


