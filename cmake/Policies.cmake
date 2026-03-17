# File holds global policies for building
#

function(require_out_of_source_build)
  if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(
      FATAL_ERROR
        "In-source builds are not allowed. Please create a separate build directory."
    )
  endif()
endfunction(require_out_of_source_build)
