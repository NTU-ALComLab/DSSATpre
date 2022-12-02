# Create the DepQBF library with integrated preprocessor Bloqqer
add_library(minisat STATIC
        core/Dimacs.h
        core/Solver.h
        core/SolverTypes.h
        mtl/Sort.h
        mtl/Vec.h
        mtl/Heap.h
        mtl/Alg.h
        mtl/Alloc.h
        mtl/IntTypes.h
        mtl/Map.h
        mtl/Queue.h
        mtl/XAlloc.h
        simp/SimpSolver.h
        utils/Options.h
        utils/ParseUtils.h
        utils/System.h
        core/Solver.cc
        utils/Options.cc
        utils/System.cc
        simp/SimpSolver.cc
)

#target_compile_definitions(minisat PRIVATE SOLVER NLOG LIBBLOQQER)
target_include_directories(minisat PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

install(TARGETS minisat
        RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin
        LIBRARY DESTINATION ${CMAKE_BINARY_DIR}/bin
        ARCHIVE DESTINATION ${CMAKE_BINARY_DIR}/bin
       )
