file(GLOB find_use_files "*.cmake")
install(
  FILES ${find_use_files} cmake_uninstall.cmake.in
  DESTINATION share/ackroyd-riyait/cmake
  COMPONENT core)
