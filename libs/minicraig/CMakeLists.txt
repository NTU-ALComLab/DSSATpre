# Create the DepQBF library with integrated preprocessor Bloqqer
add_library(minicraig STATIC
	Alg.h
	BasicHeap.h
	BoxedVec.h
	Heap.h
	Map.h
	Queue.h 
	SimpSolver.h
	Solver.h
	SolverTypes.h
	Sort.h
	Vec.h
	craig_define.hpp
	minisat_bmc.hpp
	minisat_craig.hpp
	parallel_define.hpp
	proof_define.hpp
	simple_aig.hpp
	simple_aig_inline.hpp
	time_check.hpp
	types.hpp
	SimpSolver.C
	Solver.C
	minisat_craig.cpp
	simple_aig.cpp 
	time_check.cpp
)
#target_compile_definitions(minicraig PUBLIC __STDC_LIMIT_MACROS __STDC_FORMAT_MACROS)
get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
target_include_directories(minicraig PUBLIC ${PARENT_DIR})

install(TARGETS minicraig
        RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin
        LIBRARY DESTINATION ${CMAKE_BINARY_DIR}/bin
        ARCHIVE DESTINATION ${CMAKE_BINARY_DIR}/bin
       )

