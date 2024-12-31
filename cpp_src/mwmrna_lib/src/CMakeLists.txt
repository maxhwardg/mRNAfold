set(FILES
    vienna_params.hpp
    rna.hpp rna.cpp
    energy_model.hpp
    energy_calculator.hpp
    structure.hpp structure.cpp
    energy_semiring.hpp
    logging.cpp logging.hpp
    protein.hpp protein.cpp
    vienna_param_strings.hpp
    codon_table_strings.hpp
    parallel.hpp
    codon_fold.hpp
    codon_trace.hpp
    md_array.hpp
    codon_graph.hpp
    codon_graph_fold.hpp
    codon_graph_trace.hpp)

find_package(TBB REQUIRED)

add_library(mwmrna_lib STATIC ${FILES})
target_link_libraries(mwmrna_lib PRIVATE TBB::tbb)

