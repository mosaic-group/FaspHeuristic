//
// Created by gonciarz on 2019-03-18.
//

#ifndef GRAPHIO_H
#define GRAPHIO_H


#include "graph.h"
#include <dirent.h>

namespace Graph::IO {
    /**
    * Reads data provided as adjacency list (fist number in a line is a src vertex,
    * and next number(s) are outgoing veritces)
    *
    * <b>NOTE:</b> no error handlihng is implemented, if file is not correct this function will fail
    *
    * @param aFileName inputFileName
    * @return created graph
    */
    template<typename VERTEX_TYPE = uint16_t, template<typename> class GRAPH_TYPE = GraphMap>
    static Graph<VERTEX_TYPE, GRAPH_TYPE> graphFromFile(const std::string &aFileName) {
        std::ifstream infile(aFileName);
        Graph<VERTEX_TYPE, GRAPH_TYPE> graph;

        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream iss(line);
            std::vector<std::string> tokens;
            std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), back_inserter(tokens));

            // We expect at least vertex without any outgoing connections
            if (tokens.size() >= 1) {
                typename Graph<VERTEX_TYPE>::VertexId src = std::stoi(tokens[0]);
                graph.addVertexSafe(src);

                // Add edges to destination vertices if defined
                for (size_t i = 1; i < tokens.size(); ++i) {
                    typename Graph<VERTEX_TYPE>::VertexId dst = std::stoi(tokens[i]);
                    graph.addVertexSafe(dst);
                    graph.addEdge(src, dst);
                }
            }
        }

        return graph;
    }

    template<typename VERTEX_TYPE = uint16_t, template<typename> class GRAPH_TYPE = GraphMap>
    void graphToFile(const std::string &aFileName, const Graph<VERTEX_TYPE, GRAPH_TYPE> &aGraph) {
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
     * Reads solution from "*.fas.txt" files to get exact capacity
     * @param aFileName full (with path) filename
     * @return capacity of cut edges
     */
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

    enum FileType {FT_REGULAR_FILE = DT_REG, FT_DIRECTORY = DT_DIR};

    /**
     * Returns list of files in given directory
     * @param aDir directory path
     * @param aFileType type of file (regular or directory)
     * @return list of found files (empty if there is no files of given type or directory path does not exist)
     */
    auto getFilesInDir(const std::string &aDir, const FileType aFileType = FT_REGULAR_FILE) {
        std::vector<std::string> files;

        DIR *dir = opendir(aDir.c_str());
        if (dir != NULL) {
            struct dirent *ent;
            while ((ent = readdir (dir)) != NULL) {
                if (ent->d_type == aFileType) {
                    // If we are looking for directories skip "." and ".."
                    if (aFileType == FT_DIRECTORY && (!strcmp(ent->d_name, ".") || !strcmp(ent->d_name, ".."))) continue;
                    files.emplace_back(ent->d_name);
                }
            }
        }

        return files;
    }

}


#endif
