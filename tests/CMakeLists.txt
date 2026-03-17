# tests/CMakeLists.txt
#

# --------------------------------------------------------------------------- #
# Compile and add C++ (GTest) tests
# --------------------------------------------------------------------------- #
include(FetchContent)

set(BUILD_GMOCK
    OFF
    CACHE BOOL "" FORCE)
set(INSTALL_GTEST
    OFF
    CACHE BOOL "" FORCE)
set(gtest_force_shared_crt
    ON
    CACHE BOOL "" FORCE)

# NOTE: Commit hash needs to updated regularly!
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/5a9c3f9e8d9b90bbbe8feb32902146cb8f7c1757.zip
)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

file(GLOB_RECURSE TEST_SRC "*.cc")
add_executable(ackroyd_riyait_tests ${TEST_SRC})
target_include_directories(ackroyd_riyait_tests PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Link with gtest_main to get an auto-generated main()
target_link_libraries(ackroyd_riyait_tests GTest::gtest_main ar_lib)

gtest_discover_tests(ackroyd_riyait_tests)

