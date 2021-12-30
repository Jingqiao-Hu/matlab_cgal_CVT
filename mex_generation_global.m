%% Build Cpp MEX File with CGAL library
% Build a single Cpp program into a MEX file. 
clc; clear;
opt_load = 1; 
%% specify path of headers, libs
if opt_load ==0
    include1 = ['-I', 'C:\dev\vcpkg\installed\x64-windows\include'];
    include2 = ['-I', 'C:\dev\vcpkg\installed\x64-windows\include\eigen3'];
    include3 = ['-I', 'F:\library\libigl\include'];

    lib1 = ['-L', 'C:\dev\vcpkg\installed\x64-windows\lib'];
    lib11 = ['-l', 'mpfr.lib'];
    lib12 = ['-l', 'mpir.lib'];

    lib2 = ['-L', 'C:\dev\vcpkg\installed\x64-windows\lib\auxiliary\gmp\lib'];
    lib21 = ['-l', 'libgmp-10.lib'];
    lib22 = ['-l', 'libmpfr-4.lib'];

    mex('-v', include1, include2, lib1,lib11,lib1,lib12, 'voronoi_global.cpp');
end
%% test
load bnodes
[edges_cell, nodes_cell] = voronoi_global(bnodes, 256, 128);

figure; voronoi(bnodes(:,1),bnodes(:,2)); hold on;
for ele = 1:size(bnodes, 1)
    nodes = nodes_cell{ele};
    scatter(nodes(:,1),nodes(:,2)); hold on;
end