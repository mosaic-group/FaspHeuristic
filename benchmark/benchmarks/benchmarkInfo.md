# BENCHMARKS
Description of what is where in directories.

### FASP-benchmarks
Generated graphs used in benchmarks. FASP of this graphs is precomputed using exact ILP solver.
 
### benchILPvsHEURISTIC
Benchmarks TIGHT-CUT and GR heuristics using precomputed graphs from FASP-benchmarks directory.

### benchGraphsFromPaper1
Benchmarks for graphs from:
"An exact method for the minimum feedback arc set problem" Ali Baharev et al.

### benchWeights
Results of TIGHT-CUT and GR for weighted graphs generated randomly with known FASP size.

### benchXL
Some benchmarks for (extra)large Erdos-Renyi graphs run for:
* large FASP size of 200 edges
* large graphs with density 5% and |FASP| = 20 with number of vertices in range 100-1000 and number of edges from 495-499950

### benchScalability
Some results showing how time needed by solver as:
- density (|e|/|v|) vs FASP size
- density (|e|/|v|) vs number of vertices

### benchIsoCut
Benchmarks ability of pure ISO-CUT to solve the FASP problem.
It shows percentage of cut FASP edges vs FASP size. 
(NOTE: generated benchmarks have old naming, previous name of ISO-CUT was... SuperALgorithm)

