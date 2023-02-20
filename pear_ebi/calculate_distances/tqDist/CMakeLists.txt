cmake_minimum_required(VERSION 2.8)

set(CMAKE_LEGACY_CYGWIN_WIN32 0)

project(tqDist CXX C)

enable_testing()

add_subdirectory(tqDist)

set(CPACK_PACKAGE_NAME "tqDist")
set(CPACK_PACKAGE_VENDOR "BiRC - Bioinformatics Research Center, Aarhus University, Denmark")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

if (CMAKE_SIZEOF_VOID_P MATCHES "8")
  set(CPACK_SYSTEM_NAME x86_64)
  set(CPACK_TOPLEVEL_TAG x86_64)
else (CMAKE_SIZEOF_VOID_P MATCHES "8")
  set(CPACK_SYSTEM_NAME x86_32)
  set(CPACK_TOPLEVEL_TAG x86_32)
endif (CMAKE_SIZEOF_VOID_P MATCHES "8")

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "tqDist")
set(CPACK_RESOURCE_FILE_LICENSE "COPYING")
set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_IGNORE_FILES
  "\\\\./Makefile$"
  "tqDist/Makefile$"
  "install_manifest\\\\.txt$"
  "CMakeFiles"
  "CTestTestfile\\\\.cmake$"
  "cmake_install\\\\.cmake$"
  "svn_ignore\\\\.txt$"
  "_CPack_Packages"
  "CPackSourceConfig.cmake"
  "CPackConfig.cmake"
  "CMakeCache\\\\.txt$"
  "\\\\.svn"
  "\\\\.gz$"
  "\\\\.zip$"
  "\\\\.pyc$"
  "\\\\.a$"
  "\\\\.so$"
  "\\\\.dylib$"
  "_wrap\\\\.cxx$"
  "callgrind\\\\.out"
  "~$"
  "cmake_install\\\\.cmake"
  "bin"
  "compile\\\\ R\\\\.txt"
  "CMakeClean\\\\.sh"
  "Testing"
  "\\\\.o$"
  "tqDist/test_quartet$"
  "tqDist/test_triplet$"
  "tqDist/triplet_dist$"
  "tqDist/quartet_dist$"
  "tqDist/all_pairs_triplet_dist$"
  "tqDist/all_pairs_quartet_dist$"
  "tqDist/tqDist.egg-info"
  "tqDist/build"
  "tqDist/dist"

  ${CPACK_SOURCE_IGNORE_FILES}
  )

set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

set(CPACK_RESOURCE_FILE_README  "${CMAKE_SOURCE_DIR}/README")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")

SET(CPACK_NSIS_DISPLAY_NAME "tqDist")
SET(CPACK_NSIS_HELP_LINK "http://birc.au.dk/software/tqDist")
SET(CPACK_NSIS_URL_INFO_ABOUT "http://birc.au.dk/software/tqDist")
SET(CPACK_NSIS_CONTACT "cstorm@birc.au.dk")

set(CPACK_NSIS_MODIFY_PATH ON)
set(CPACK_NSIS_EXTRA_INSTALL_COMMANDS "\${EnvVarUpdate} \"$0\" \"PATH\"  \"A\" \"HKCU\" \"${INSTDIR}\\\\bin\"")

include(CPack)

#This adds the 'dist' target
add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
