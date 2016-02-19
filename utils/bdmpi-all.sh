#!/bin/bash

SRCDIR=/scratch/jeremy/src
V1=0.1.7
V2=0.2.3
BITMASK=31
NT=2

# v0.1.x
cd ${SRCDIR}/sbma ; checkout v${V1} ; make clean ; make install
cd ${SRCDIR}/bdmpi ; make clean ; make install
cd ${SRCDIR}/parmetis-4.0.3-bdmpi ; make clean ; make install
cd ${SRCDIR}/csplatt ; make clean ; make install

./bdmpi.sh -b ${BITMASK} -r 1 -p pr1d -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p pr1d -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p pr1d -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p parmetis -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p parmetis -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p parmetis -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p splatt -m 786432 -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p splatt -m 786432 -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p splatt -m 786432 -t ${NT}

# v0.2.x
cd ${SRCDIR}/sbma ; checkout v${V2} ; make clean ; make install
cd ${SRCDIR}/bdmpi ; make clean ; make install
cd ${SRCDIR}/parmetis-4.0.3-bdmpi ; make clean ; make install
cd ${SRCDIR}/csplatt ; make clean ; make install

./bdmpi.sh -b ${BITMASK} -r 1 -p pr1d -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p pr1d -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p pr1d -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p parmetis -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p parmetis -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p parmetis -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p splatt -m 786432 -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p splatt -m 786432 -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 1 -p splatt -m 786432 -t ${NT}

./bdmpi.sh -b ${BITMASK} -r 7 -p pr1d -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 7 -p pr1d -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 7 -p pr1d -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 7 -p parmetis -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 7 -p parmetis -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 7 -p parmetis -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 7 -p splatt -m 786432 -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 7 -p splatt -m 786432 -t ${NT}
./bdmpi.sh -b ${BITMASK} -r 7 -p splatt -m 786432 -t ${NT}
