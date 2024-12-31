add_executable(fold_codon_graph src/fold_codon_graph.cpp)

set(targets fold_codon_graph)

foreach(target IN LISTS targets)
  target_include_directories(${target} PUBLIC ${CMAKE_SOURCE_DIR}/mwmrna_lib/src)
  target_link_libraries(${target} mwmrna_lib)
endforeach()

