# cmake file for configuring global variables used throughout build process
#

# ----------------------------------------------------------------------------#
# Set up RCANA cmake paths
# ----------------------------------------------------------------------------#
get_filename_component(ACKROYD_RIYAIT_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}" DIRECTORY)

cmake_path(APPEND ACKROYD_RIYAIT_ROOT_DIR "source" OUTPUT_VARIABLE ACKROYD_RIYAIT_SOURCE_DIR)
cmake_path(APPEND ACKROYD_RIYAIT_ROOT_DIR "include" OUTPUT_VARIABLE ACKROYD_RIYAIT_INCLUDE_DIR)
