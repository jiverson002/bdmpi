gdem [![Build Status](https://travis-ci.org/Cecca/gdem.png)](https://travis-ci.org/Cecca/gdem)
==============================================================================================

Graph Diameter Estimator on MPI

How to build
------------

From the root directory of the project

    mkdir build
    cd build
    cmake ..
    make

How to run tests
----------------

In the build directory

    make
    make test

How to generate documentation
-----------------------------

In the build directory

    make doc

The generated documentation will be available in `doc/html`

The documentation of the script `scripts/pair2adjacency.py` is available 
[here](http://www.dei.unipd.it/~ceccarel/gdem-docs/pair2adjacency.html).

Use of `scan-build` to detect potential bugs
--------------------------------------------

The `clang` static analizer can be used to detect potential bugs
in the code. First install the latest version of the `clang` compiler.
Then, from the project root directory

    mkdir check
    cd check
    scan-build cmake ..
    scan-build make

The static analyzer will output something similar to this

    scan-build: Using '/usr/bin/clang' for static analysis
    [ 12%] Building C object src/CMakeFiles/mpiHello.dir/mpi_hello.c.o
    Linking C executable ../bin/mpiHello
    [ 12%] Built target mpiHello
    [ 25%] Building C object src/CMakeFiles/serial.dir/hll_counter.c.o
    [ 37%] Building C object src/CMakeFiles/serial.dir/serial.c.o
    /home/matteo/Development/c/gdem/src/serial.c:127:3: warning: Memory is never released; potential leak of memory pointed to by 'f'
      char *line = strtok(f, "\n");
        ^
    /home/matteo/Development/c/gdem/src/serial.c:172:10: warning: Memory is never released; potential leak of memory pointed to by 'n.neighbours'
      return 0;
             ^
    2 warnings generated.
    Linking C executable ../bin/serial
    [ 37%] Built target serial
    [ 50%] Building C object src/CMakeFiles/testHllCounter.dir/hll_counter.c.o
    [ 62%] Building C object src/CMakeFiles/testHllCounter.dir/test_hll_counter.c.o
    Linking C executable ../bin/testHllCounter
    [ 62%] Built target testHllCounter
    [ 75%] Building C object src/CMakeFiles/testParser.dir/parser.c.o
    [ 87%] Building C object src/CMakeFiles/testParser.dir/graph.c.o
    [100%] Building C object src/CMakeFiles/testParser.dir/test_parser.c.o
    Linking C executable ../bin/testParser
    [100%] Built target testParser
    scan-build: 2 bugs found.
    scan-build: Run 'scan-view /tmp/scan-build-2013-05-22-1' to examine bug reports.

To view the report generated by the tool, run the command suggested in the output, 
in this case

    scan-view /tmp/scan-build-2013-05-22-1

How to configure Eclipse
------------------------

To generate project files for the Eclipse IDE, run the following from 
the project root directory

    cd ..
    mkdir gdem_eclipse
    cd gdem_eclipse
    cmake -G "Eclipse CDT4 - Unix Makefiles" ../gdem

Then you can import the project into Eclipse from the directory
`gdem_eclipse`. Please note that due to the way Eclipse works (at least
Eclipse Indigo) the Eclipse directory _must_ be a sibiling of the project's
root directory in order to perform out of source builds (that is highly desirable).