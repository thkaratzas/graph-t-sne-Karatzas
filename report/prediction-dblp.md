# Prediction: coPapersDBLP on Apple M1 (written before the M1 run)

Machine: Apple M1, 4 Performance + 4 Efficiency cores, 8 GB, macOS 26.6.2,
gcc 16.2.0 / OpenCilk 3.0, -O2, m = 16, theta = 0.5.

Basis: measured minimum times on com-Amazon on the same machine
(results/times-MacBook-Air-Thanasis.csv), scaled by the cost model of each phase.
Note: coPapersDBLP was run sequentially earlier on the 2-vCPU Codespace while
developing; those times were not used for this prediction.

## Sizes

| graph        | n       | nnz        | mean degree |
|--------------|---------|------------|-------------|
| com-Amazon   | 334 863 | 1 851 744  | 5.5         |
| coPapersDBLP | 540 486 | 30 491 458 | 56.4        |

Scale factors from Amazon to DBLP:
- nnz: 16.47
- n: 1.61
- n log(n/m): 1.69

## Cost model per phase

| phase | model        | why                                                        |
|-------|--------------|------------------------------------------------------------|
| P     | ~ nnz        | unweighted, so uniform rows (no bisection); one search per edge |
| attr  | ~ nnz        | a few flops and one random read of y_j per edge            |
| build | ~ n log(n/m) | each tree level touches every point                        |
| com   | ~ n          | one pass over points and cells                             |
| walk  | ~ n log(n/m) | about log n cells opened per point                         |
| grad  | ~ n          | one pass over 2n numbers                                   |

Parallel times: the DBLP sequential prediction divided by the speedup
measured on com-Amazon (OpenMP, dynamic,256) at the same thread count.

## Predicted times (seconds)

| phase | sequential | 4 threads | 8 threads |
|-------|------------|-----------|-----------|
| P     | 1.31       | 0.37      | 0.30      |
| attr  | 0.061      | 0.022     | 0.022     |
| build | 0.022      | 0.010     | 0.009     |
| com   | 0.0013     |           |           |
| walk  | 1.81       | 0.55      | 0.45      |
| grad  | 0.0010     |           |           |
| total | 3.2        | 0.96      | 0.78      |

## Measured (to be filled in after the run)

| phase | sequential | 4 threads | 8 threads | difference |
|-------|------------|-----------|-----------|------------|
| P     |            |           |           |            |
| attr  |            |           |           |            |
| build |            |           |           |            |
| com   |            |           |           |            |
| walk  |            |           |           |            |
| grad  |            |           |           |            |
