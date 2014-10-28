# Build program.
add_executable(helloworld helloworld.c)
add_executable(p2p p2p.c)
add_executable(pscan pscan.c)
add_executable(bcast bcast.c)
add_executable(allgather allgather.c)
add_executable(reduce reduce.c)
add_executable(commtest commtest.c)
add_executable(alltoall alltoall.c)
add_executable(gather gather.c)
add_executable(scatter scatter.c)
add_executable(test_MPI_ANY_SOURCE test_MPI_ANY_SOURCE.c)
foreach(prog helloworld
             p2p
             pscan
             bcast
             allgather
             reduce
             commtest
             scatter
             gather
             alltoall
             test_MPI_ANY_SOURCE
             )
  target_link_libraries(${prog} bdmpi GKlib rt dl)
  add_dependencies(${prog} ${GKlib_includes})
endforeach(prog)

# if(BDMP_INSTALL)
#   install(TARGETS helloworld p2p
#     RUNTIME DESTINATION bin)
# endif()
