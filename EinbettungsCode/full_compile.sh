#!/bin/sh
#rm -rf *.o *.eps *.png *.pdf *.bak
#rm -rf *.lst *.map
rm -f polygon_path_fail*.eps ./failed*.eps
cd CGAL_skeleton
cmake . -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=/home/ankylotech/Downloads/CGAL-6.0
make
cd ..
aloft -L-lstdc++ -g tests embedding linesegmentintersection -L./CGAL_skeleton/CMakeFiles/skeleton_construction.dir/skeleton_construction.cpp.o -L/usr/lib/x86_64-linux-gnu/libgmpxx.so -L/usr/lib/x86_64-linux-gnu/libmpfr.so -L/usr/lib/x86_64-linux-gnu/libgmp.so
aloft -L-lstdc++ -g geometric embedding linesegmentintersection -L./CGAL_skeleton/CMakeFiles/skeleton_construction.dir/skeleton_construction.cpp.o -L/usr/lib/x86_64-linux-gnu/libgmpxx.so -L/usr/lib/x86_64-linux-gnu/libmpfr.so -L/usr/lib/x86_64-linux-gnu/libgmp.so
aloft -L-lstdc++ -g combinatorial embedding linesegmentintersection -L./CGAL_skeleton/CMakeFiles/skeleton_construction.dir/skeleton_construction.cpp.o -L/usr/lib/x86_64-linux-gnu/libgmpxx.so -L/usr/lib/x86_64-linux-gnu/libmpfr.so -L/usr/lib/x86_64-linux-gnu/libgmp.so
aloft -L-lstdc++ -g graphics embedding linesegmentintersection -L./CGAL_skeleton/CMakeFiles/skeleton_construction.dir/skeleton_construction.cpp.o -L/usr/lib/x86_64-linux-gnu/libgmpxx.so -L/usr/lib/x86_64-linux-gnu/libmpfr.so -L/usr/lib/x86_64-linux-gnu/libgmp.so
