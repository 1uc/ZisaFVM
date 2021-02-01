target_sources(zisa_generic_obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/down_sample_reference.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/numerical_experiment.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/numerical_experiment_factory.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/polytrope.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/rayleigh_taylor.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/smooth_bubble.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/stellar_convection.cpp
)

if(ZISA_HAS_MPI)

endif()
