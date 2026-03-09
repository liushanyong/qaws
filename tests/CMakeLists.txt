add_executable(
  qaws_tests
  test_runner.c
  # Test sources will be added here as they are implemented.
)

target_link_libraries(qaws_tests PRIVATE Qaws::qaws)
target_include_directories(qaws_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src/qaws")
set_target_properties(qaws_tests PROPERTIES C_STANDARD 11 C_STANDARD_REQUIRED YES C_EXTENSIONS NO)

if(MSVC)
  target_compile_definitions(qaws_tests PRIVATE _CRT_SECURE_NO_WARNINGS NOMINMAX)
endif()

add_test(NAME qaws_tests COMMAND qaws_tests)
