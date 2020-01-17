# Benchmark graphs for FAS solvers

This directory contains graphs for benchmarking feedback arc set solvers.

Each folder contains a set of graphs from a specific class, and if it could be computed, one of its minimal feedback arc sets (or its size):

 - `.al` files contain the graphs in adjacency list format. The first number on a line is a vertex, the rest are its neighbours.

 - `.mfas` files contain a minimum feedback arc set for the corresponding graph. Each line contains two vertices, representing an arc.

 - `.mfsize` files contain the size of the minimum feedback arc sets (i.e. a single integer). These files are used only in special cases, when the size could be computed (e.g. based on analytical results), but an actual feedback arc set could not.

 - `.timing` files contain the approximate time (in seconds) needed by SageMath's ILP-based solver to solve this instance on the machine `Particulator`.

 - `.timeout` files are present if the ILP-based solver was unable to finish within the time constraint. There is always a corresponding `.timing` file which contains the amount of time after which the solver was cancelled.

File names typically have the form `Name-ID-params-VertexCount-EdgeCount`, where `ID` is a numeric ID unique to each graph within its set (but not across sets). `params` is only present for some types of graphs (e.g. de Bruijn and Kautz graphs), and represents a set of parameters that uniquely defines the graph.

### random

Random directed graphs sampled from Erdős-Rényi $G(n,m)$ models.

Solver time limit: 30 min.

### random-oriented

Random graphs sampled from undirected Erdős-Rényi $G(n,m)$ models, with each of their edges randomly oriented as `<-` or `->`.

Solver time limit: 30 min.

### tournaments

Random tournaments, i.e. randomly oriented complete undirected graphs.

Solver time limit: 30 min.

### planar-triangulations

Planar triangulations, i.e. planar graphs with the maximum possible number of edges. Such a graph on $V$ vertices has exactly $3V-6$ edges.

These are generated by:

 - Generating uniformly distributed points on the surface of a sphere.
 - Computing their convex hull.
 - Using the edges of the convex hull as the underlying undirected graph.

Then edges are oriented randomly as `<-` or `->` with probability 1/2 each.

All graphs in this directory are solved.

### planar-triangulations-2

These are like the graphs in `planar-triangulations`, but edges are oriented as `<-`, `->` or `<->` with probability 1/3 each.

All graphs in this directory are solved.

### delaunay-3d

The underlying undirected graph is a Delaunay triangulation obtained from uniformly sampled random points in a cube. Its edges are oriented randomly as `<-` or `->` with probability 1/2 each.

Solver time limit: 30 min.

### delaunay-3d-2

These are like the graphs in `delaunay-3d`, but edges are oriented as `<-`, `->` or `<->` with probability 1/3 each.

Solver time limit: 30 min.

### small-world

We take a planar triangulation without reciprocal edges and randomly rewire a fraction of its edges (while avoiding creating reciprocal edges). Parameter `-percentage-` is the fraction of rewired edges.

### de-bruijn

[De Bruijn graphs](https://en.wikipedia.org/wiki/De_Bruijn_graph) with self-loops removed. Parameters `-m-n-` denote a graph on `m` characters and string length `n`.

[There are exact results available](http://dx.doi.org/10.1016/j.camwa.2009.10.021) for the minimum feedback _vertex_ set size of de Bruijn graphs.

### kautz

[Kautz graphs](https://en.wikipedia.org/wiki/Kautz_graph). Parameters `-m-n-` denote a graph on `m+1` characters and string length `n+1`.

Notes on Kautz graphs:

 - For the sake of this discussion, $K_m^n$ denotes the Kautz graph on $m+1$ symbols and string length $n+1$. (Warning: This notation is different from the one on Wikipedia!)
 - $K_m^n$ has in- and out-degrees $m$. Therefore, it has $m$ times more edges than vertices.
 - $K_m^n$ has $(m+1)m^n$ vertices and $(m+1)m^{n+1}$ edges.
 - The line graph of $K_n^m$ is $K_n^{m+1}$. This fact connects the feedback arc set and feedback vertex set problems of this family of graphs.
 - Kautz graphs with $m>1$ have no isolated cycles (tested empirically).
 - $K_m^n$ has edge connectivity $m$ (tested empirically), which implies the previous statement.

### kautz-2

This folder contains a different set of Kautz graphs with their corresponding feedback sizes, computed based on the results of the paper ["Feedback numbers of Kautz digraphs"](https://dx.doi.org/10.1016/j.disc.2006.09.010). This paper gives exact MFAS size results for $K_m^n$ with any $m$ and $n \le 5$, and approximate one for $n > 5$. Only `.mfsize` files are available. **Warning:** Some of these graphs are very large, with up to 300,000 edges.