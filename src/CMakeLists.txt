# gather all .cc files from src directory
file(GLOB_RECURSE SOURCES "*.cc")

# create RCANA library
add_library(ar_lib)
target_compile_features(ar_lib PUBLIC cxx_std_23)
target_include_directories(ar_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_sources(ar_lib PRIVATE ${SOURCES})

# create executable and link to library
add_executable(ackroyd_riyait main.cc)
target_link_libraries(ackroyd_riyait PRIVATE ar_lib)
