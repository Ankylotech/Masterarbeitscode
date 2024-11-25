#!/bin/sh
aloft -L-lstdc++ -g graphics embedding linesegmentintersection -L./CGAL_skeleton/CMakeFiles/skeleton_construction.dir/skeleton_construction.cpp.o -L/usr/lib/x86_64-linux-gnu/libgmpxx.so -L/usr/lib/x86_64-linux-gnu/libmpfr.so -L/usr/lib/x86_64-linux-gnu/libgmp.so
./graphics
cp ./polygon_path_fail_total2.eps ../Arbeit/img/DemoSkeleton.eps 
cp ./polygon_path_fail_total4.eps ../Arbeit/img/DemoLinePartition.eps 
cp ./polygon_path_fail_total10.eps ../Arbeit/img/DemoCircleLine.eps 

