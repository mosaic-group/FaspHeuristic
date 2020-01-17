# Summaries

This directory contains summarized information about the dataset and benchmark results, in machine readable formats.

Feel free to add new files to this directory.

### `summary.json`

This is a JSON file organized hierarchically in the same way as the `data` directory, containing the following information for each graph:

 - `VertexCount`: Vertex count.
 - `EdgeCount`: Edge count.
 - `CycleRank`: Cycle rank. This is computed as (# vertices) - (# edges in strongly connected components) + (# strongly connected components).
 - `MinFASSize`: Exact minimum feedback arc set size.
 - `ILPTiming`: Time taken to compute the minimum FAS size by SageMath's ILP-based FAS solver on `Particulator`.
 - `EadesLinSmythFASSize`: FAS size computed by the Eades–Lin–Smyth heuristic.
 - `EdgeConnectivityVector`: A vector of integers where the 1st element is the number of vertices, and the `k`th element is the total number of vertices in `k`-connected components, excluding size-1 components. The length of this vector minus 1 is the connectivity of the most-connected subgraph.

See the `plots` directory for various visualizations of these quantities.

### `summary-superalg.json`

- `VertexCount`
- `EdgeCount`
- `MinFASSize`
- `SARemovableCount`: number of arcs that could be immediately removed with the super-algorithm.

The super-algorithm iteratively removes the first edge it finds that has some isolated cycles.

See the `plots` directory for various visualizations of these quantities.
