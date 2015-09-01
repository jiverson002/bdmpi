# Build program.
add_executable(oc_samplesort samplesort.c)
add_executable(oc_sphkmeans sphkmeans.c)
add_executable(oc_pr pr.c)
add_executable(oc_mf2d mf2d.c)
foreach(prog oc_samplesort oc_sphkmeans oc_pr oc_mf2d)
  target_link_libraries(${prog} bdmpi GKlib rt m dl)
  #add_dependencies(${prog} ${GKlib_includes})
endforeach(prog)

if(BDMP_INSTALL)
  install(TARGETS oc_samplesort oc_sphkmeans oc_pr oc_mf2d
    RUNTIME DESTINATION bin)
endif()

