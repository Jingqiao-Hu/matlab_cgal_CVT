# matlab invoke CGAL to generate CVT

generate CVT based on input seeds

return each nodes and their connections (edges) for the required seed

openmp has been used to speed up the for iteration

### input

1. seeds

2. the selected indices of seeds
3. the x-length bounding box
4. the y-length bounding box

### output

matrix of voronoi edges of the selected seeds

## library

CGAL

## usage

1. generate .mex using mex_generation_local.m

2. the corresponding examples are also given

