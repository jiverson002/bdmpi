add_executable ( graphcoloring
      graphcoloring.c
)

foreach(prog graphcoloring
             )
  target_link_libraries(${prog} bdmpi GKlib rt dl)
  add_dependencies(${prog} ${GKlib_includes})
endforeach(prog)
