# Build program.
add_executable(merge merge.c)
add_executable(vmqsort vmqsort.c)
add_executable(lockmem lockmem.c)
foreach(prog merge vmqsort lockmem) 
  target_link_libraries(${prog} bdmpi GKlib rt dl)
  #add_dependencies(${prog} ${GKlib_includes})
endforeach(prog)

# if(BDMP_INSTALL)
#   install(TARGETS merge vmqsort lockmem 
#     RUNTIME DESTINATION bin)
# endif()
