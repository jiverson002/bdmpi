# Helper modules
include(FindPerl)

set(bdmpicc_in "${CMAKE_CURRENT_SOURCE_DIR}/bdmpicc.pl.in")
set(bdmpicc_out "${CMAKE_CURRENT_BINARY_DIR}/bdmpicc")

set(bdmpicxx_in "${CMAKE_CURRENT_SOURCE_DIR}/bdmpic++.pl.in")
set(bdmpicxx_out "${CMAKE_CURRENT_BINARY_DIR}/bdmpic++")

configure_file(${bdmpicc_in} ${bdmpicc_out})
configure_file(${bdmpicxx_in} ${bdmpicxx_out})

add_custom_target(
  bdmpicc 
  COMMAND touch ${bdmpicc_out}
  DEPENDS ${bdmpicc_in}
)

add_custom_target(
  bdmpic++ 
  COMMAND touch ${bdmpicxx_out}
  DEPENDS ${bdmpicxx_in}
)


if(BDMP_INSTALL)
  install(PROGRAMS ${bdmpicc_out} DESTINATION bin)
  install(PROGRAMS ${bdmpicxx_out} DESTINATION bin)
endif()
