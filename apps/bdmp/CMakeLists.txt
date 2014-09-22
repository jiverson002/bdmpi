# Build program.
add_executable(bdmp_samplesort samplesort.c)
add_executable(bdmp_sphrbkmeans sphrbkmeans.c)
add_executable(bdmp_sphkmeans sphkmeans.c)
add_executable(bdmp_pr1d pr1d.c)
add_executable(bdmp_pr1dmt pr1dmt.c)
add_executable(bdmp_pr1drmt pr1drmt.c)
add_executable(bdmp_mf2d mf2d.c)
add_executable(bdmp_mf2dmt mf2dmt.c)
add_executable(bdmp_mf1d mf1d.c)
add_executable(bdmp_mf1d2 mf1d2.c)
add_executable(bdmp_mf1d2-mpi mf1d2-mpi.c)
add_executable(sb_sphkmeans sb_sphkmeans.c)
add_executable(sb_pr1d sb_pr1d.c)
foreach(prog bdmp_samplesort 
             bdmp_sphrbkmeans
             bdmp_sphkmeans 
             bdmp_pr1d 
             bdmp_pr1dmt
             bdmp_pr1drmt
             bdmp_mf2d
             bdmp_mf2dmt
             bdmp_mf1d
             bdmp_mf1d2
             bdmp_mf1d2-mpi
             sb_sphkmeans 
             sb_pr1d 
       )
  target_link_libraries(${prog} bdmpi GKlib rt m dl)
  add_dependencies(${prog} ${GKlib_includes})
endforeach(prog)

if(BDMP_INSTALL)
  install(TARGETS bdmp_samplesort 
                  bdmp_sphrbkmeans 
                  bdmp_sphkmeans
                  bdmp_pr1d 
                  bdmp_pr1dmt 
                  bdmp_pr1drmt
                  bdmp_mf2d
                  bdmp_mf2dmt
                  bdmp_mf1d
                  bdmp_mf1d2
                  bdmp_mf1d2-mpi
                  sb_sphkmeans
                  sb_pr1d
    RUNTIME DESTINATION bin)
endif()
