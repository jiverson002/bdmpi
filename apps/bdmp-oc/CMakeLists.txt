# Build program.
add_executable(bdmp-oc_sphkmeans sphkmeans.c)
add_executable(bdmp-oc_pr1d pr1d.c)
add_executable(bdmp-oc_pr1d-dv pr1d-dv.c)
add_executable(bdmp-oc_mf2d mf2d.c)
add_executable(bdmp-oc_mf1d mf1d.c)
foreach(prog 
            bdmp-oc_sphkmeans
            bdmp-oc_pr1d 
            bdmp-oc_pr1d-dv
            bdmp-oc_mf2d
            bdmp-oc_mf1d
       )
  target_link_libraries(${prog} bdmpi GKlib rt m dl)
  #add_dependencies(${prog} ${GKlib_includes})
endforeach(prog)

if(BDMP_INSTALL)
  install(TARGETS bdmp-oc_sphkmeans 
                  bdmp-oc_pr1d 
                  bdmp-oc_pr1d-dv 
                  bdmp-oc_mf2d
                  bdmp-oc_mf1d
    RUNTIME DESTINATION bin)
endif()
