rm -rf $1
svn up

# Get sources 
# svn export http://dminers.dtc.umn.edu/svn/programs/karypis/bdmpi/pola/trunk $1
# svn export http://dminers.dtc.umn.edu/svn/libs/GKlib/trunk $1/GKlib
svn export http://dminers.dtc.umn.edu/svn/programs/karypis/bdmpi/pola/branches/jeremy $1
svn export http://dminers.dtc.umn.edu/svn/libs/GKlib/trunk $1/GKlib

# Build documentation
cd $1/doxygen
doxygen html.doxyfile
mv html ../doc
doxygen latex.doxyfile
cd latex
make
cp refman.pdf ../../doc
cd ../../../

# Remove build instructions that are not distributed
mv $1/CMakeLists.txt $1/xxx
grep -v DISTRM $1/xxx > $1/CMakeLists.txt
rm -rf $1/xxx

# Remove files/directories that are not distributed
rm -rf $1/apps
rm -rf $1/config
rm -rf $1/doxygen
rm -rf $1/misc
rm -rf $1/results
rm -rf $1/TODO
rm -rf $1/utils

# Create tar file 
tar -czf $1.tar.gz $1

# Remove distribution directory 
rm -rf $1
