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

    % mex -v CXXFLAGS='$CXXFLAGS -Wall' '-IC:\dev\vcpkg\installed\x64-windows\include' COMPFLAGS='$COMPFLAGS /openmp' 'mexCVT.cpp';
    % mex('-v', include1, include2, include3, lib1,lib11,lib1,lib12,lib2,lib21,lib2,lib22, 'voronoi_bbx.cpp');

    mex('-v', include1, include2, lib1,lib11,lib1,lib12, 'voronoi_local.cpp');
end
%% test
load bnodes

figure
for i = 1:size(bnodes, 1)
    [edges, nodes] = voronoi_local(bnodes, i, 256, 128);
    scatter(nodes(:,1),nodes(:,2)); hold on;
end