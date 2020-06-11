## TIGHT-CUT heuristic

Library implementing **TIGHT-CUT** heuristic for solving **FASP** (Feedback Arc Set Problem).

If you use this library in an academic context, please cite the following paper:
- *Michael Hecht, Krzysztof Gonciarz, Szabolcs Horvát* [Tight Localizations of Feedback Sets](https://arxiv.org/abs/2001.01440)

### How to obtain the code
Clone the git repository by running the following command from the command line:
```bash
git clone https://git.mpi-cbg.de/mosaic/FaspHeuristic.git
```

If you prefer to access GitLab with SSH keys, use:
```bash
git clone git@git.mpi-cbg.de:mosaic/FaspHeuristic.git
```

### How to compile
The library has these options available:
```
FASP_HEURISTIC_INSTALL (default: ON)
FASP_HEURISTIC_BUILD_TOOLS (default: ON)
FASP_HEURISTIC_BUILD_TESTS (default: OFF)
FASP_HEURISTIC_BUILD_BENCHMARKS (default: OFF)
```

To build the tools, execute the following commands (starting in the root directory after downloading):
```bash
mkdir build
cd build
cmake ..
make
```

The command-line tools will be placed in `build/tools`.

To install the library in the default location, execute:
```bash
make install
```

To install the library in a different location, provide the prefix for the installation directory when running `cmake`:
```bash
cmake -DCMAKE_INSTALL_PREFIX=/user/specified/install/directory ..
```

### Command-line tools
Two command-line tools are provided for convenience, one for unweighted and one for weighted graphs:
- ```tightCut```
- ```tightCutWeighted```

Both tools read the graph from the standard input.

Vertices are specified as non-negative integers. Any integer may be used, regardless of the total number of vertices.

The format of the input graph for the unweighted case is an adjacency list:
```plain
<src vertex> <optional dst vertex 1> <optional dst vertex 2> ...
```

For example:
```plain
1 2 3 6
2 3
3 4 6
4 1
5
6 3
```

If above definition of graph is saved in a file with name ```graph.def``` then to calculate FASP run:
```plain
$ tightCut < graph.def
# Input graph: #vertices=6 #edges=8
# FASP size (#edges to cut): 2
# Feedback arcs:
6 3
3 4
```

The format of the input graph for the weighted case is as follows. Each line must either contain a single vertex:  
```plain
<vertex>
```
Or an edge followed by the edge weight:
```plain
<src vertex> <dst vertex> <weight of edge>
```

For example:
```plain
1 2 5
2 3 4.5
8
3 1 10
```

Output with above definition of graph:
```plain
$ tightCutWeighted < graph.def
# Input graph: #vertices=4 #edges=3
# FASP size (#edges to cut): 1, capacity of edges: 4.5
# Feedback arcs:
2 3
```

### How to use in your project
TIGHT-CUT library can be used as a submodule of your project or can be used as a regular installed library.

### Simple example
Simple project using installed TIGHT-CUT library requires only two files: `main.cpp` with the code and a simple configuration `CMakeLists.txt` file.

**main.cpp**
```c++
#include <iostream>

#include "FaspTightCut/graph.h"
#include "FaspTightCut/graphFasp.h"

int main() {
    // Create graph by adding vertices and edges to the graph
    using VERTEX_TYPE = int;
    Graph::Graph<VERTEX_TYPE> g;
    for (int i = 0; i < 8; ++i) g.addVertex(i);
    g.addEdge({0, 1});
    g.addEdge({1, 2});
    g.addEdge({2, 3});
    g.addEdge({3, 4});
    g.addEdge({4, 5});
    g.addEdge({5, 0});

    g.addEdge({3, 6});
    g.addEdge({6, 2});

    g.addEdge({5, 7});
    g.addEdge({7, 4});

    // Add capacity to edges
    auto ep = Graph::Ext::getEdgeProperties(g, 1);
    ep[{2, 3}] = 3;
    ep[{3, 6}] = 1;
    ep[{6, 2}] = 1;
    ep[{4, 5}] = 3;
    ep[{5, 7}] = 1;
    ep[{7, 4}] = 1;

    // Run TIGHT-CUT heuristic
    auto [capacity, removedEdges, saEdgesCnt, saRndEdgesCnt, redRndEdgesCnt] = Graph::Fasp::tightCut<true, true, VERTEX_TYPE, int>(g, ep);

    // Print result
    std::cout << "Capacity of FASP edges: " << capacity << std::endl;

    return 0;
}
```

**CMakeLists.txt**
```cmake
cmake_minimum_required(VERSION 3.15)
project(TightCutExample LANGUAGES CXX)

find_package(FaspTightCut REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main FaspTightCut::FaspTightCut)
```
