target_sources(zisa_obj
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/plotting_steps.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/rate_of_change.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/runge_kutta.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/simulation_clock.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/step_rejection.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/time_integration_factory.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/time_keeper.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/time_keeper_factory.cpp
)

