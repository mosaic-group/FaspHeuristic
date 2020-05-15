## TIGHT-CUT heuristic

Library implementing **TIGHT-CUT** heuristic for solving **FASP** (Feedback Arc Set Problem).

If you use this library in an academic context, please cite the following paper:
- *Michael Hecht, Krzysztof Gonciarz, Szabolcs Horvát* [Tight Localizations of Feedback Sets](https://arxiv.org/abs/2001.01440)

### How to clone code
Run from command line following command:
```bash
git clone https://github.com/krzysg/FaspHeuristic.git
```

of if you prefere to access github with SSH keys:
```bash
git clone git@github.com:krzysg/FaspHeuristic.git
```

### How to compile
Library has these options available:
```
FASP_HEURISTIC_INSTALL (default: ON)
FASP_HEURISTIC_BUILD_TOOLS (default: ON)
FASP_HEURISTIC_BUILD_TESTS (default: OFF)
FASP_HEURISTIC_BUILD_BENCHMARKS (default: OFF)
```

To compile, build tools and install it in default location please execute following commands (starting in root directory after downloading):
```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
make install
```

To install library in non default location in `cmake` line provide also prefix for install direcotry:
```bash
-DCMAKE_INSTALL_PREFIX=/user/specified/install/directory
```
### Tools
Two command line tools are available one for unweighted and one for weighted graphs:
- ```tightCut```
- ```tightCutWeighted```

Both accept graph provided from standard input.

Format of input graph for unweighted case is a:  
```<src vertex> <optional dst vertex1> <optional dst vertex 2> ...```

For example:
```asm
1 3 4 5
3 6
7
6 7
7 1
```

If above definition of graph is saved in a file with name ```graph.def``` then to calculate FASP run:
```asm
$ tightCut < graph.def 
Input graph: #vertices=6 #edges=6
FASP size (#edges to cut): 1
```

Format of input graph for weighted case is:  
```<vertex>```  
```<src vertex> <dst vertex> <weight of edge>```  

For example:
```asm
1 2 5
2 3 4
8
3 1 10
```

Output with above definition of graph:
```asm
$ tightCutWeighted < graph.def 
Input graph: #vertices=4 #edges=3
FASP size (#edges to cut): 1, capacity of edges: 4
```

### How to use in your project
TIGHT-CUT library can be used as a submodule of your project or can be used as a regular installed library.

### Simple example
Simple project using installed TIGHT-CUT library requires only two files `main.cpp` with a code and simple configuration `CMakeLists.txt` file.

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

