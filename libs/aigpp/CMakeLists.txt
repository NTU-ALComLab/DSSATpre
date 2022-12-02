# Create the aigpp library
add_library(aigpp STATIC
	AigerParser.hh
	ComputedTable.hh
	CuddTools.hh
	DynamicRewritingTable.hh
	Edge.hh
	Edge.icc
	EdgeRef.hh
	EdgeRef.icc
	EdgeRef_operations.icc
	EquivalenceChecker.hh
	ExtRefTable.hh
	FourInputFunction.hh
	FRAIGManager.hh
	InternalEdgeRef.hh
	InternalEdgeRef.icc
	Manager_export.icc
	Manager.hh
	Manager.icc
	NodeHashMap.hh
	NodeHashMap.icc
	Node.hh
	Node.icc
	NPNClass.hh
	Primes.hh
	RewritingManager.hh
	RewritingManager.icc
	Settings.hh
	Settings.icc
	SimulationTable.hh
	SimVector.hh
	SimVector.icc
	StatManager.hh
	UniqueTable.hh
	VariableAssignment.hh
	VariableAssignment.icc
	VarSupport.hh
	VarSupport.icc
	AigerParser.cc
	ComputedTable.cc
	CuddTools.cc
	DynamicRewritingTable.cc
	Edge.cc
	EdgeRef.cc
	EdgeRef_operations.cc
	EquivalenceChecker.cc
	ExtRefTable.cc
	FourInputFunction.cc
	FRAIGManager.cc
	InternalEdgeRef.cc
	InternalEdgeRef_decompose.cc
	InternalEdgeRef_polarity.cc
	Manager_and.cc
	Manager_assignment.cc
	Manager_BDDQuantify.cc
	Manager_BDDSweeping.cc
	Manager.cc
	Manager_export.cc
	Manager_dcmin.cc
	Manager_fraig.cc
	Manager_garbage.cc
	Manager_generateSim.cc
	Manager_quantify.cc
	Manager_quantify_craig.cc
	Manager_recursive.cc
	Manager_ref.cc
	Manager_sat.cc
	Manager_sim.cc
	Node.cc
	NPNClass.cc
	RewritingManager.cc
	Settings.cc
	SimulationTable.cc
	SimVector.cc
	StaticRewritingTable.cc
	StatManager.cc
	UniqueTable.cc
	VariableAssignment.cc
	VarSupport.cc
           )

target_link_libraries(aigpp PUBLIC lrabsutil)
target_link_libraries(aigpp PUBLIC minisat)
target_link_libraries(aigpp PUBLIC minicraig)
target_link_libraries(aigpp PUBLIC cudd)
target_compile_definitions(aigpp PUBLIC USE_MINICRAIG)
get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
target_include_directories(aigpp PUBLIC ${PARENT_DIR})


install(TARGETS aigpp
        RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin
        LIBRARY DESTINATION ${CMAKE_BINARY_DIR}/bin
        ARCHIVE DESTINATION ${CMAKE_BINARY_DIR}/bin
       )
