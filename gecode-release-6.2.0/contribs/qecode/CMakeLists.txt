cmake_minimum_required(VERSION 3.13)
project(Qecode)

set(CMAKE_CXX_STANDARD 11)

set(GECODE_LIB
        gecodeminimodel
        gecodeint
        gecodesearch
        gecodekernel
        gecodesupport
)

set(CXX_SRCS
        myDom.cpp
        QSpace.cpp QSpace.h
        OptVar.cpp OptVar.h
        QProblem.cpp QProblem.h
        QCSPPlusUnblockable.cpp QCSPPlusUnblockable.h
        qecode.h
        qsolver_parallel.cpp qsolver_parallel.h
        qsolver_qcop.cpp qsolver_qcop.h
        QCSPSolver.cpp QCSPSolver.h
        qsolver_unblockable.cpp qsolver_unblockable.h
        Strategy.cpp Strategy.h
        StrategyNode.cpp StrategyNode.h
        UnblockableBranching.h
        UnblockableViewValBranching.cpp UnblockableViewValBranching.h
        vartype.h
        Work.cpp Work.h
        WorkComparators.h
        Worker.cpp Worker.h
        WorkManager.cpp WorkManager.h
        WorkComparators.cpp
)
add_compile_options(-g -O2 )
add_library(Qecode SHARED ${CXX_SRCS})
# Link Gecode to your library
target_link_libraries(Qecode ${GECODE_LIB})

#Examples
function(addQecodeExample EXAMPLE_NAME)
    add_executable(${EXAMPLE_NAME} examples/${EXAMPLE_NAME}.cpp)
    target_link_libraries(${EXAMPLE_NAME} Qecode ${GECODE_LIB})
endfunction()

#addQecodeExample(MatrixGame)
#addQecodeExample(NetworkPricing1)
#addQecodeExample(NetworkPricing2)
#addQecodeExample(NimFibo)
#addQecodeExample(OldUnitTests)
#addQecodeExample(TestQCOP)

