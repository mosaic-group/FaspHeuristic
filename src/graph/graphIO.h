#ifndef GRAPHIO_H
#define GRAPHIO_H


#include "graph.h"
#include "graphExt.h"
#include <dirent.h>

namespace Graph::IO {
    /**
     * Reads data provided as adjacency list (fist number in a line is a src vertex,
     * and next number(s) are outgoing veritces). Example file content.
     *
     * 1 4
     * 3 10
     * 11 5
     *
     * <b>NOTE:</b> no error handlihng is implemented, if file is not correct this function will fail
     *
     * @param aFileName inputFileName
     * @return read graph
     */
    template<typename VERTEX_TYPE = uint16_t>
    static Graph<VERTEX_TYPE> graphFromFile(const std::string &aFileName) {
        std::ifstream infile(aFileName);
        Graph<VERTEX_TYPE> graph;

        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            std::vector<std::string> tokens;
            std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), back_inserter(tokens));

            // We expect at least vertex without any outgoing connections
            if (tokens.size() >= 1) {
                typename Graph<VERTEX_TYPE>::Vertex src = std::stoi(tokens[0]);
                graph.addVertexSafe(src);

                // Add edges to destination vertices if defined
                for (size_t i = 1; i < tokens.size(); ++i) {
                    typename Graph<VERTEX_TYPE>::Vertex dst = std::stoi(tokens[i]);
                    graph.addVertexSafe(dst);
                    graph.addEdge(src, dst);
                }
            }
        }

        return graph;
    }

    /**
     * Save provided graph as adjacency list (fist number in a line is a src vertex,
     * and next number(s) are outgoing veritces). Example file content.
     *
     * 1 4
     * 3 10
     * 11
     *
     * <b>NOTE:</b> no error handlihng is implemented, if file is not correct this function will fail
     *
     * @param aFileName outputFileName
     * @param aGraph graph to save
     */
    template<typename VERTEX_TYPE = uint16_t>
    void graphToFile(const std::string &aFileName, const Graph<VERTEX_TYPE> &aGraph) {
        std::ofstream outfile(aFileName, std::ios_base::out); // overwrite output file if exists

        // In a case of map implementation (or any other using unordered containter) output files  look
        // bad without sorting - do it so!
        auto allVertices = aGraph.getVertices();
        std::sort(allVertices.begin(), allVertices.end());

        for (auto &v : allVertices) {
            outfile << v;
            for (auto &ov : aGraph.getOutVertices(v)) {
                outfile << "\t" << ov;
            }
            outfile << "\n";
        }
    }

    /**
     * Reads solution from files to get exact capacity. Solution is in form of list FASP edges:
     *
     * 1 5
     * 4 3
     * 3 1
     *
     * @param aFileName full (with path) filename
     * @return capacity of cut edges
     */
    [[maybe_unused]]
    static int solutionFromFile(const std::string &aFileName) {
        std::ifstream infile(aFileName);

        int capacity = 0;

        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            std::vector<std::string> tokens;
            std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), back_inserter(tokens));

            // Each correct line of solution file contain src and destination vertex
            if (tokens.size() == 2) ++capacity;
        }

        return capacity;
    }

    /**
     * Reads timing from files to get time of exact solver (no check for errors done!). File should contain
     * exactly one number representing number of seconds (can be float number).
     * @param aFileName full (with path) filename
     * @return timing
     */
    [[maybe_unused]]
    static double timingFromFile(const std::string &aFileName) {
        std::ifstream infile(aFileName);

        double num = 0.0;
        infile >> num;

        return num;
    }

    enum FileType {
        FT_REGULAR_FILE = DT_REG, FT_DIRECTORY = DT_DIR
    };

    /**
     * Returns list of files in given directory
     * @param aDir directory path
     * @param aFileType type of file (regular or directory)
     * @return list of found files (empty if there is no files of given type or directory path does not exist)
     */
    [[maybe_unused]]
    auto getFilesInDir(const std::string &aDir, const FileType aFileType = FT_REGULAR_FILE) {
        std::vector<std::string> files;

        DIR *dir = opendir(aDir.c_str());
        if (dir != NULL) {
            struct dirent *ent;
            while ((ent = readdir(dir)) != NULL) {
                if (ent->d_type == aFileType) {
                    // If we are looking for directories skip "." and ".."
                    if (aFileType == FT_DIRECTORY && (!strcmp(ent->d_name, ".") || !strcmp(ent->d_name, ".."))) continue;
                    files.push_back(std::string{ent->d_name});
                }
            }
            closedir(dir);
        }
        return files;
    }

    /**
     * Reads data provided as adjacency list with weights. File should contain 3 numbers in each row
     * representing: srcVertex dstVertex weightOfEdge. For example:
     *
     * 1 4 3
     * 1 3 10
     * 4 3 2
     *
     * <b>NOTE:</b> no error handlihng is implemented, if file is not correct this function will fail
     *
     * @param aFileName inputFileName
     * @param aGraph graph to save
     * @param aWeights weights of edges in aGraph
     */
    template<typename VERTEX_TYPE = uint16_t, typename EDGE_PROP_TYPE>
    void graphWithWeightsToFile(const std::string &aFileName, const Graph<VERTEX_TYPE> &aGraph, const Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights) {
        std::ofstream outfile(aFileName, std::ios_base::out); // overwrite output file if exists

        // In a case of map implementation (or any other using unordered containter) output files  look
        // bad without sorting - do it so!
        auto allVertices = aGraph.getVertices();
        std::sort(allVertices.begin(), allVertices.end());

        for (auto &v : allVertices) {
            bool hasEdges = false;
            for (auto &ov : aGraph.getOutVertices(v)) {
                hasEdges = true;
                outfile << v;
                outfile << "\t" << ov << "\t" << aWeights.at({v, ov}) << "\n";
            }
            if (!hasEdges) {
                // vertex without edges
                outfile << v << "\n";
            }
        }
    }

    /**
     * Reads data provided as adjacency list with weights. File should contain 3 numbers in each row
     * representing: srcVertex dstVertex weightOfEdge. For example:
     *
     * 1 4 3
     * 1 3 10
     * 4 3 2
     *
     * <b>NOTE:</b> no error handlihng is implemented, if file is not correct this function will fail
     *
     * @param aFileName inputFileName
     * @return pair {created graph, edgePropertiesWithWeighs}
     */
    template<typename VERTEX_TYPE = uint16_t, typename EDGE_PROP_TYPE>
    static auto graphWithWeightsFromFile(const std::string &aFileName) {
        std::ifstream infile(aFileName);
        Graph<VERTEX_TYPE> graph;
        Ext::EdgeProperties<VERTEX_TYPE, EDGE_PROP_TYPE> weights;
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            std::vector<std::string> tokens;
            std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), back_inserter(tokens));

            // We expect 1 (just vertex) or 3 numbers in line (srcVertex dstVertex edgeWeight)
            if (tokens.size() == 1 || tokens.size() == 3) {
                // Read vertices and create edge
                typename Graph<VERTEX_TYPE>::Vertex src = std::stoi(tokens[0]);
                graph.addVertexSafe(src);

                if (tokens.size() == 3) {
                    typename Graph<VERTEX_TYPE>::Vertex dst = std::stoi(tokens[1]);
                    graph.addVertexSafe(dst);
                    graph.addEdge(src, dst);

                    // Read weight and update properites of graph
                    EDGE_PROP_TYPE w = static_cast<EDGE_PROP_TYPE>(std::stod(tokens[2]));
                    weights[{src, dst}] = w;
                }
            }
        }

        return std::pair{graph, weights};
    }

}


#endif
