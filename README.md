# computational-geometry-cpp

C++ port of the computational-geometry examples. Each algorithm is implemented directly and checked against a second implementation or an independent geometric property.

| Example | What it computes | Implementation |
| --- | --- | --- |
| `closest_pair` | Minimum Euclidean pair | divide-and-conquer, brute force, KD-tree |
| `convex_hull` | Convex hull vertices | monotone chain in 2D, incremental hull in 3D |
| `laplacian_smoothing` | Mesh Laplacian smoothing | uniform or inverse-distance weights |
| `segment_intersection` | Intersections of line segments | Bentley-Ottmann sweep line |
| `minimum_bounding_box` | Oriented minimum bounding box | rotating calipers, hull-face search in 3D |
| `nearest_neighbors` | k nearest neighbors | brute force, KD-tree |
| `shoelace` | Area of a polygon ring | shoelace sum |
| `delaunay` | Delaunay triangulation | Bowyer-Watson with an exact in-circle test |

The 2D minimum-area rectangle is exact: one side lies on a hull edge. The 3D box is the minimum over frames built from hull triangles. Laplacian smoothing of a 2D point set with no faces uses the Delaunay mesh. `weight="inverse_distance"` uses `1 / (distance + 1e-12)`.

## Build

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

## Use

```bash
./build/geometry closest-pair --n 400 --dim 2 --out closest.svg
./build/geometry convex-hull --n 200 --out hull.svg
./build/geometry obb --n 200 --angle 40 --out box.svg
./build/geometry intersections --n 30 --out segments.svg
./build/geometry neighbors --n 80 --k 2 --out neighbors.svg
./build/geometry smooth --n 10 --iterations 8 --out smooth.svg
./build/geometry shoelace --n 9 --out shoelace.svg
./build/geometry delaunay --n 40 --out delaunay.svg
```

`./build/geometry --n 400` is the closest-pair demo. Figures are SVG. Demo point sets use `std::mt19937`, so they are not the same clouds as the Python package's NumPy generator.

The headers under `include/cg` are the library API: `closest_pair`, `convex_hull`, `laplacian_smooth`, `segment_intersections`, `minimum_bounding_box`, `nearest_neighbors`, `shoelace`, and `delaunay`.
