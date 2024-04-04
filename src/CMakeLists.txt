# Copyright (c) 2016-2017 Hong Xu

# This file is part of WCSPLift.

# WCSPLift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# WCSPLift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with WCSPLift.  If not, see <http://www.gnu.org/licenses/>.

# Use C++ 11 standard
set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

set(GECODE_LIB
        gecodeminimodel
        gecodeint
        gecodesearch
        gecodekernel
        gecodesupport
)

set(CXX_SRCS
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/myDom.cpp
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/QSpace.cpp ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/QSpace.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/OptVar.cpp ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/OptVar.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/QProblem.cpp ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/QProblem.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/QCSPPlusUnblockable.cpp ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/QCSPPlusUnblockable.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/qecode.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/qsolver_parallel.cpp ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/qsolver_parallel.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/qsolver_qcop.cpp ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/qsolver_qcop.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/QCSPSolver.cpp ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/QCSPSolver.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/qsolver_unblockable.cpp ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/qsolver_unblockable.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/Strategy.cpp ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/Strategy.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/StrategyNode.cpp ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/StrategyNode.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/UnblockableBranching.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/UnblockableViewValBranching.cpp ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/UnblockableViewValBranching.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/vartype.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/Work.cpp ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/Work.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/WorkComparators.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/Worker.cpp ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/Worker.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/WorkManager.cpp ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/WorkManager.h
        ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode/WorkComparators.cpp
)


# Boost
find_package(Boost REQUIRED COMPONENTS graph program_options system timer iostreams filesystem)
set(all_depends ${all_depends} ${Boost_LIBRARIES})
#message(Boost_INCLUDE_DIRS ="${Boost_INCLUDE_DIRS}")
include_directories(${Boost_INCLUDE_DIRS})



# BLAS
find_package(BLAS REQUIRED)
set(all_depends ${all_depends} ${BLAS_LIBRARIES})
#message(BLAS_INCLUDE_DIRS="${BLAS_INCLUDE_DIRS}")
include_directories("/opt/local/include/openblas")

# GUROBI
#find_package(GUROBI REQUIRED)
#set(all_depends ${all_depends} ${GUROBI_LIBRARIES})
#include_directories(${GUROBI_INCLUDE_DIRS})

# ZLIB
find_package(ZLIB)
set (all_depends ${all_depends} ${ZLIB_LIBRARIES}) 

# LZMA
find_package(LibLZMA)
if(LIBLZMA_FOUND)
    set (all_depends ${all_depends} ${LIBLZMA_LIBRARIES}) 
else(LIBLZMA_FOUND)
    add_compile_definitions("-DNO_LZMA")
endif(LIBLZMA_FOUND)


add_compile_options(-g -O2 )

#link_directories(${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode)
#foreach(lib IN
add_library(qeCodeLib SHARED ${CXX_SRCS})
target_link_libraries(qeCodeLib ${GECODE_LIB})

#set_target_properties(qeCodeLib PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode)

#message(make source dir ="${CMAKE_SOURCE_DIR}")

set (all_depends ${all_depends} ${qeCodeLib})
set (all_depends ${all_depends} ${GECODE_LIB})
#message(all depends fin ="${all_depends}")

find_library(qeCodeLibDY libQecode.dylib PATHS ${CMAKE_SOURCE_DIR}/gecode-release-6.2.0/contribs/qecode)
message(${qeCodeLibDY})
set (all_depends ${all_depends} ${qeCodeLibDY})

# For now, we assume Gurobi is always available
#add_definitions(-DHAVE_GUROBI)
#set(GUROBI_VERSION "95" CACHE STRING "Gurobi version.")
#set(GUROBI_CXX_LIBRARY "libgurobi_c++" CACHE STRING "You may need to adapt the value of the C++ library according to your system. Please read into your Gurobi library directory (usually ${GUROBI_HOME}/lib) for available files.")
#set(GUROBI_LIBRARIES "${GUROBI_CXX_LIBRARY}" "libgurobi${GUROBI_VERSION}")

set(wcsp_BINSRCS
    main.cpp)
#LinearProgramSolver.cpp
#LinearProgramSolverGurobi.cpp

add_executable(wcsp ${wcsp_BINSRCS})
target_link_libraries(wcsp ${all_depends})

set_target_properties(wcsp PROPERTIES
    OUTPUT_NAME wcsp
    LINKER_FLAG "${BLAS_LINKER_FLAGS}")

install(TARGETS wcsp RUNTIME DESTINATION bin)
