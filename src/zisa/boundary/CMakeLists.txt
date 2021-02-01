target_sources(zisa_generic
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/boundary_condition_factory.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/frozen_boundary_condition.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/halo_exchange_bc.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/no_boundary_condition.cpp
)

if(ZISA_HAS_MPI)

endif()
