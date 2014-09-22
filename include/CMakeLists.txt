# Install the header file as both mpi.h and bdmpi.h

if(BDMP_INSTALL)
  install(FILES bdmpi.h DESTINATION include/${HINSTALL_PATH})
  install(FILES bdmpi.h DESTINATION include/${HINSTALL_PATH} RENAME mpi.h)
#  install(FILES bdmpif.h DESTINATION include)
endif()
