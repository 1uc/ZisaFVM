target_sources(zisa_generic_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/down_sample_reference.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/numerical_experiment.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/numerical_experiment_factory.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/polytrope.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/rayleigh_taylor.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/smooth_bubble.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/stellar_convection.cpp
)

if(ZISA_HAS_MPI)

endif()
