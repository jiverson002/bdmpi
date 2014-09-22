# Add this directory for internal users.
include_directories(.)

# Find sources.
file(GLOB libbdmp_sources *.c)

# Build libbdmp.
add_library(bdmpi ${BDMP_LIBRARY_TYPE} ${common_sources} ${libbdmp_sources})
add_dependencies(bdmpi ${GKlib_includes})

if(BDMP_INSTALL)
  install(TARGETS bdmpi
    LIBRARY DESTINATION lib/${LINSTALL_PATH}
    RUNTIME DESTINATION lib/${LINSTALL_PATH}
    ARCHIVE DESTINATION lib/${LINSTALL_PATH})
endif()
