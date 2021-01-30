target_sources(zisa_generic_obj
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/rate_of_change.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/runge_kutta.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/step_rejection.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/time_integration_factory.cpp
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}/time_keeper_factory.cpp
)

if(ZISA_HAS_MPI)

endif()
