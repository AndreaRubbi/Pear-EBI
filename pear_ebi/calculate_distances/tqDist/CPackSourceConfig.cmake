# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. The list of available CPACK_xxx variables and their associated
# documentation may be obtained using
#  cpack --help-variable-list
#
# Some variables are common to all generators (e.g. CPACK_PACKAGE_NAME)
# and some are specific to a generator
# (e.g. CPACK_NSIS_EXTRA_INSTALL_COMMANDS). The generator specific variables
# usually begin with CPACK_<GENNAME>_xxxx.


set(CPACK_BINARY_DEB "OFF")
set(CPACK_BINARY_FREEBSD "OFF")
set(CPACK_BINARY_IFW "OFF")
set(CPACK_BINARY_NSIS "OFF")
set(CPACK_BINARY_RPM "OFF")
set(CPACK_BINARY_STGZ "ON")
set(CPACK_BINARY_TBZ2 "OFF")
set(CPACK_BINARY_TGZ "ON")
set(CPACK_BINARY_TXZ "OFF")
set(CPACK_BINARY_TZ "ON")
set(CPACK_BUILD_SOURCE_DIRS "/home/andrear/tqDist;/home/andrear/tqDist")
set(CPACK_CMAKE_GENERATOR "Unix Makefiles")
set(CPACK_COMPONENT_UNSPECIFIED_HIDDEN "TRUE")
set(CPACK_COMPONENT_UNSPECIFIED_REQUIRED "TRUE")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_FILE "/usr/share/cmake-3.22/Templates/CPack.GenericDescription.txt")
set(CPACK_DEFAULT_PACKAGE_DESCRIPTION_SUMMARY "tqDist built using CMake")
set(CPACK_GENERATOR "TGZ;ZIP")
set(CPACK_IGNORE_FILES "\\./Makefile$;tqDist/Makefile$;install_manifest\\.txt$;CMakeFiles;CTestTestfile\\.cmake$;cmake_install\\.cmake$;svn_ignore\\.txt$;_CPack_Packages;CPackSourceConfig.cmake;CPackConfig.cmake;CMakeCache\\.txt$;\\.svn;\\.gz$;\\.zip$;\\.pyc$;\\.a$;\\.so$;\\.dylib$;_wrap\\.cxx$;callgrind\\.out;~$;cmake_install\\.cmake;bin;compile\\ R\\.txt;CMakeClean\\.sh;Testing;\\.o$;tqDist/test_quartet$;tqDist/test_triplet$;tqDist/triplet_dist$;tqDist/quartet_dist$;tqDist/all_pairs_triplet_dist$;tqDist/all_pairs_quartet_dist$;tqDist/tqDist.egg-info;tqDist/build;tqDist/dist")
set(CPACK_INSTALLED_DIRECTORIES "/home/andrear/tqDist;/")
set(CPACK_INSTALL_CMAKE_PROJECTS "")
set(CPACK_INSTALL_PREFIX "/usr/local")
set(CPACK_MODULE_PATH "")
set(CPACK_NSIS_CONTACT "cstorm@birc.au.dk")
set(CPACK_NSIS_DISPLAY_NAME "tqDist")
set(CPACK_NSIS_DISPLAY_NAME_SET "TRUE")
set(CPACK_NSIS_EXTRA_INSTALL_COMMANDS "${EnvVarUpdate} "$0" "PATH"  "A" "HKCU" "\\bin"")
set(CPACK_NSIS_HELP_LINK "http://birc.au.dk/software/tqDist")
set(CPACK_NSIS_INSTALLER_ICON_CODE "")
set(CPACK_NSIS_INSTALLER_MUI_ICON_CODE "")
set(CPACK_NSIS_INSTALL_ROOT "$PROGRAMFILES")
set(CPACK_NSIS_MODIFY_PATH "ON")
set(CPACK_NSIS_PACKAGE_NAME "tqDist")
set(CPACK_NSIS_UNINSTALL_NAME "Uninstall")
set(CPACK_NSIS_URL_INFO_ABOUT "http://birc.au.dk/software/tqDist")
set(CPACK_OUTPUT_CONFIG_FILE "/home/andrear/tqDist/CPackConfig.cmake")
set(CPACK_PACKAGE_DEFAULT_LOCATION "/")
set(CPACK_PACKAGE_DESCRIPTION_FILE "/usr/share/cmake-3.22/Templates/CPack.GenericDescription.txt")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "tqDist")
set(CPACK_PACKAGE_FILE_NAME "tqDist-1.0.1")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "tqDist 1.0.1")
set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "tqDist 1.0.1")
set(CPACK_PACKAGE_NAME "tqDist")
set(CPACK_PACKAGE_RELOCATABLE "true")
set(CPACK_PACKAGE_VENDOR "BiRC - Bioinformatics Research Center, Aarhus University, Denmark")
set(CPACK_PACKAGE_VERSION "1.0.1")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_RESOURCE_FILE_LICENSE "/home/andrear/tqDist/COPYING")
set(CPACK_RESOURCE_FILE_README "/home/andrear/tqDist/README")
set(CPACK_RESOURCE_FILE_WELCOME "/usr/share/cmake-3.22/Templates/CPack.GenericWelcome.txt")
set(CPACK_RPM_PACKAGE_SOURCES "ON")
set(CPACK_SET_DESTDIR "OFF")
set(CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set(CPACK_SOURCE_IGNORE_FILES "\\./Makefile$;tqDist/Makefile$;install_manifest\\.txt$;CMakeFiles;CTestTestfile\\.cmake$;cmake_install\\.cmake$;svn_ignore\\.txt$;_CPack_Packages;CPackSourceConfig.cmake;CPackConfig.cmake;CMakeCache\\.txt$;\\.svn;\\.gz$;\\.zip$;\\.pyc$;\\.a$;\\.so$;\\.dylib$;_wrap\\.cxx$;callgrind\\.out;~$;cmake_install\\.cmake;bin;compile\\ R\\.txt;CMakeClean\\.sh;Testing;\\.o$;tqDist/test_quartet$;tqDist/test_triplet$;tqDist/triplet_dist$;tqDist/quartet_dist$;tqDist/all_pairs_triplet_dist$;tqDist/all_pairs_quartet_dist$;tqDist/tqDist.egg-info;tqDist/build;tqDist/dist")
set(CPACK_SOURCE_INSTALLED_DIRECTORIES "/home/andrear/tqDist;/")
set(CPACK_SOURCE_OUTPUT_CONFIG_FILE "/home/andrear/tqDist/CPackSourceConfig.cmake")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "tqDist-1.0.1")
set(CPACK_SOURCE_TOPLEVEL_TAG "x86_64-Source")
set(CPACK_STRIP_FILES "")
set(CPACK_SYSTEM_NAME "x86_64")
set(CPACK_THREADS "1")
set(CPACK_TOPLEVEL_TAG "x86_64-Source")
set(CPACK_WIX_SIZEOF_VOID_P "8")

if(NOT CPACK_PROPERTIES_FILE)
  set(CPACK_PROPERTIES_FILE "/home/andrear/tqDist/CPackProperties.cmake")
endif()

if(EXISTS ${CPACK_PROPERTIES_FILE})
  include(${CPACK_PROPERTIES_FILE})
endif()
