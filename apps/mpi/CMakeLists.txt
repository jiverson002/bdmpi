# Build program.
add_executable(mpi_samplesort samplesort.c)
add_executable(mpi_sphrbkmeans sphrbkmeans.c)
add_executable(mpi_sphkmeans sphkmeans.c)
add_executable(mpi_pr1d pr1d.c)
add_executable(mpi_mf2d mf2d.c)
add_executable(mpi_mf1d2 mf1d2.c)
#add_executable(mpi_rw mpi_rw.c)
#foreach(prog mpi_samplesort mpi_sphkmeans mpi_rw)
foreach(prog mpi_samplesort 
             mpi_sphkmeans 
             mpi_sphrbkmeans 
             mpi_pr1d 
             mpi_mf2d 
             mpi_mf1d2)
  target_link_libraries(${prog} GKlib m)
  #add_dependencies(${prog} ${GKlib_includes})
endforeach(prog)

set(CMAKE_C_COMPILER mpicc)

if(BDMP_INSTALL)
#  install(TARGETS mpi_samplesort mpi_sphkmeans mpi_rw
  install(TARGETS mpi_samplesort 
                  mpi_sphrbkmeans 
                  mpi_sphkmeans 
                  mpi_pr1d 
                  mpi_mf2d 
                  mpi_mf1d2
    RUNTIME DESTINATION bin)
endif()
