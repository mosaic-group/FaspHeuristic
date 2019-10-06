//
// Created by gonciarz on 9/25/19.
//

#ifndef FASPHEURISTIC_DINIC3RDPARTY_H
#define FASPHEURISTIC_DINIC3RDPARTY_H

/* Dinic algorithm for the Max-Flow problem,
implemented for simple oriented graphs */

#include "graph/graph.h"
#include "graph/graphExt.h"

#include<stdio.h>
#include<assert.h>
#include<stdlib.h>
//#include<stdbool.h>
#include<string.h> //for memset
#include<limits.h> //for INT_MAX

namespace Graph {

    class Dinic {

        static const int uv = 0;
        static const int vu = 1;

        static const int MAXN = 1000;
        static const int MAXM = 2000;
        static const int MAXFLOW = 100000000; //max flow on one edge
        static const int DEBUG2 = 0;

//We assume that every edge also has its opposite edge. If not,
//we add it with capacity 0.

        typedef struct edge {
            int u, v; //endpoints
            int c[2]; //capacities of uv and vu
            int f[2]; //flows on uv and vu
        } edge_t;

        void sort_edges(edge_t *G, int n, int m) {
            int *bucket = (int *) calloc(n + 1, sizeof(int));
            edge_t *by2nd = (edge_t *) malloc(m * sizeof(edge_t));
            for (int i = 0; i < m; i++) {
                bucket[G[i].v]++;
            }
            for (int i = 1; i <= n; i++) {
                bucket[i] += bucket[i - 1];
            }
            for (int i = 0; i < m; i++) {
                by2nd[--bucket[G[i].v]] = G[i];
            }
            memset(bucket, 0, (n + 1) * sizeof(int));
            for (int i = 0; i < m; i++) {
                bucket[G[i].u]++;
            }
            for (int i = 1; i <= n; i++) {
                bucket[i] += bucket[i - 1];
            }
            for (int i = m - 1; i >= 0; i--) {
                G[--bucket[by2nd[i].u]] = by2nd[i];
            }
            free(bucket);
            free(by2nd);
        }

// merges all edges between a pair of vertices into one.
// assumes that all edges satisfy edge.u < edge.v
// retval: updated number of edges after merge
        int merge_edges(edge_t *G, int n, int m) {
            sort_edges(G, n, m);
            if (DEBUG2) {
                printf("before merge\n");
                for (int i = 0; i < m; i++) {
                    printf("%d %d %d %d\n", G[i].u, G[i].v, G[i].c[uv], G[i].c[vu]);
                }
            }
            int _m = 0;
            for (int i = 1; i < m; i++) {
                if (G[i].u == G[_m].u && G[i].v == G[_m].v) {
                    G[_m].c[uv] += G[i].c[uv];
                    G[_m].c[vu] += G[i].c[vu];
                } else {
                    G[++_m] = G[i];
                }
            }
            m = (m == 0) ? m : _m + 1;
            if (DEBUG2) {
                printf("after merge\n");
                for (int i = 0; i < m; i++) {
                    printf("%d %d %d %d\n", G[i].u, G[i].v, G[i].c[uv], G[i].c[vu]);
                }
            }
            return m;
        }

        typedef struct dedge { //directed edge in the reserve network
            int v; //target endpoint
            int c; //capacity
            int f; //flow
            int idx; //index in the original edge list
        } dedge_t;

        typedef struct bedge { //edge leading back in the cleaned reserve network
            int u; //target endpoint
            int idx; //index in the dedge structre of the reserve network

        } bedge_t;

        typedef struct vertex {
            int outdeg; //real out degree
            int d_idx, d_offset; //index and offset in the dedge structure of the reserve network
            int b_idx, b_offset; //index and offset in the bedge structure of the reserve network
            int lvl; //level in the reserve network
            bool alive; //flag for vertices still in the reserve network
        } vertex_t;

        int n, m; //# of vertices/edges
        int s, t; //source and sink
        int l; //distance of source and sink
        int reserve_m; //# of edges in the reserve network

        vertex_t V[MAXN];
        dedge_t E[2 * MAXM];
        bedge_t Eb[2 * MAXM];

        int path[MAXN]; //path in the cleaned reserve network

//vertex queue for various purposes
        int Q[MAXN];
        int Qsize;

        void reset_vertex_data() {
            for (int i = 0; i < n; i++) {
                V[i] = (vertex_t) {0};
                V[i].alive = true;
                V[i].lvl = -1;
            }
        }

        void create_reserve_network(edge_t *G) {
            reset_vertex_data();
            reserve_m = 0;

            for (int i = 0; i < m; i++) {
                if (G[i].c[uv] - G[i].f[uv] + G[i].f[vu] > 0) {
                    V[G[i].u].d_offset++;
                    reserve_m++;
                }
                if (G[i].c[vu] - G[i].f[vu] + G[i].f[uv] > 0) {
                    V[G[i].v].d_offset++;
                    reserve_m++;
                }
            }
            for (int i = 1; i < n; i++) {
                V[i].d_idx = V[i - 1].d_idx + V[i - 1].d_offset;
            }
            for (int i = 0; i < m; i++) {
                if (G[i].c[uv] - G[i].f[uv] + G[i].f[vu] > 0) {
                    int idx = V[G[i].u].d_idx + V[G[i].u].outdeg++;
                    E[idx].v = G[i].v;
                    E[idx].c = G[i].c[uv] - G[i].f[uv] + G[i].f[vu];
                    E[idx].f = 0;
                    // save the position of the edge in original list for final flow augmentation
                    E[idx].idx = i;
                }
                if (G[i].c[vu] - G[i].f[vu] + G[i].f[uv] > 0) {
                    int idx = V[G[i].v].d_idx + V[G[i].v].outdeg++;
                    E[idx].v = G[i].u;
                    E[idx].c = G[i].c[vu] - G[i].f[vu] + G[i].f[uv];
                    E[idx].f = 0;
                    // save the position of the edge in original list for final flow augmentation
                    E[idx].idx = i;
                }
            }
        }

        void calculate_levels(void) {
            Qsize = 0;
            Q[Qsize++] = s;
            V[s].lvl = 0;

            for (int q = 0; q < Qsize; q++) {
                assert(q < n);

                int u = Q[q];
                for (int i = V[u].d_idx; i < V[u].d_idx + V[u].d_offset; i++) {
                    int v = E[i].v;
                    if (V[v].lvl == -1) {
                        V[v].lvl = V[u].lvl + 1;
                        Q[Qsize++] = v;
                    }
                }
            }
            l = V[t].lvl;

            if (DEBUG2) {
                printf("lvls:\n");
                for (int i = 0; i < n; i++) {
                    printf("%d ", V[i].lvl);
                }
                printf("\n");
            }
        }

        void clean_dead_ends() {
            for (int q = 0; q < Qsize; q++) {
                assert(q < n);

                int u = Q[q];
                for (int i = V[u].b_idx; i < V[u].b_idx + V[u].b_offset; i++) {
                    if (E[Eb[i].idx].c > E[Eb[i].idx].f) { //if the edge is still there
                        if (--V[Eb[i].u].outdeg == 0) {
                            V[Eb[i].u].alive = false;
                            Q[Qsize++] = Eb[i].u;
                        }
                    }
                }
            }
            Qsize = 0;
        }

        bool clean_reserve_network() {
            calculate_levels();

            //if there is no path from source to sink, finish
            if (l == -1) {
                return false;
            }

            Qsize = 0;
            for (int i = 0; i < n; i++) {
                //we will recalculate the out degree
                V[i].outdeg = 0;

                //ignore vertices too far from the source
                if (V[i].lvl >= l) {
                    V[i].d_offset = 0;
                    if (i != t) {
                        V[i].alive = false;
                    }
                } else {
                    for (int j = V[i].d_idx; j < V[i].d_idx + V[i].d_offset; j++) {
                        //ignore edges that go to bad vertices or levels
                        if ((E[j].v == t || V[E[j].v].lvl < l) && V[E[j].v].lvl == V[i].lvl + 1) {
                            E[V[i].d_idx + V[i].outdeg++] = E[j];
                            //prepare back refs for dead end cleaning
                            V[E[j].v].b_offset++;
                        }
                    }

                    V[i].d_offset = V[i].outdeg; //at this point they are the same

                    if (V[i].outdeg == 0) {
                        V[i].alive = false;
                        Q[Qsize++] = i;
                    }
                }
            }

            assert(V[t].alive);

            //construct back refs for dead end cleaning
            for (int i = 1; i < n; i++) {
                V[i].b_idx = V[i - 1].b_idx + V[i - 1].b_offset;
            }
            for (int i = 0; i < n; i++) {
                V[i].b_offset = 0;
            }
            for (int i = 0; i < n; i++) {
                for (int j = V[i].d_idx; j < V[i].d_idx + V[i].d_offset; j++) {
                    //back ref for the edge i --> E[j].v
                    int idx = V[E[j].v].b_idx + V[E[j].v].b_offset++;
                    Eb[idx].u = i;
                    Eb[idx].idx = j;
                }
            }

            //clean dead ends from the queue
            clean_dead_ends();

            return true;
        }

//returns maximal augmentation that can be done on the path
int find_path() {
            if (DEBUG2) {
                printf("PHASE OF LENGTH %d\n", l);
                for (int i = 0; i < n; i++) {
                    printf("%d (%d, %d): ", i, V[i].d_offset, V[i].outdeg);
                    for (int j = V[i].d_idx; j < V[i].d_idx + V[i].d_offset; j++) {
                        if (V[E[j].v].alive) {
                            printf("%d (%d), ", E[j].v, E[j].c - E[j].f);
                        }
                    }
                    printf("\n");
                }
            }

            int u = s;
            int maxaug = MAXFLOW, aug, idx;

            for (int i = 0; i < l; i++) {

                //lazy elimination of dead vertices
                while (V[u].d_offset > 0 && !V[E[V[u].d_idx + V[u].d_offset - 1].v].alive) {
                    --V[u].d_offset;
                }
                idx = V[u].d_idx + V[u].d_offset - 1;

                //check for non-existent path
                if (V[u].outdeg == 0) {
                    return 0;
                }

                //continue along the last edge in the neighbour list
                path[i] = E[idx].v;
                u = path[i];

                //update minimal augmentation of an edge on the path
                aug = E[idx].c - E[idx].f;
                maxaug = (aug < maxaug) ? aug : maxaug;
            }

            assert(path[l - 1] == t);

            return maxaug;
        }

        void augment_blocking_flow(int aug) {
            int u = s, idx;
            for (int i = 0; i < l; i++) {
                idx = V[u].d_idx + V[u].d_offset - 1;
                E[idx].f += aug; //increase flow along the edge
                if (E[idx].f == E[idx].c) { //check if the capacity limit is reached
                    --V[u].d_offset;
                    if (--V[u].outdeg == 0) { //this is the last edge
                        V[u].alive = false; //if the vertex has no more active outgoing edges, kill it
                        Q[Qsize++] = u;
                    }
                }
                u = path[i];
            }
        }

        void find_blocking_flow() {
            int aug;
            while ((aug = find_path())) {
                augment_blocking_flow(aug);
                clean_dead_ends();
            }
        }

        void augment_flow(edge_t *G) {
            int offset;
            for (int u = 0; u < n; u++) {
                offset = (u < n - 1) ? V[u + 1].d_idx : reserve_m;
                for (int i = V[u].d_idx; i < offset; i++) {
                    if (E[i].f > 0) {
                        int idx = E[i].idx;
                        if (G[idx].u == u) {
                            G[idx].f[uv] += E[i].f;
                            if (G[idx].f[uv] > G[idx].c[uv]) {
                                G[idx].f[vu] -= G[idx].f[uv] - G[idx].c[uv];
                                G[idx].f[uv] = G[idx].c[uv];
                            }
                        } else {
                            G[idx].f[vu] += E[i].f;
                            if (G[idx].f[vu] > G[idx].c[vu]) {
                                G[idx].f[uv] -= G[idx].f[vu] - G[idx].c[vu];
                                G[idx].f[vu] = G[idx].c[vu];
                            }
                        }
                    }
                }
            }
        }
    public:
        template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE, template<typename> class GRAPH_TYPE>
        int runDinic(const Graph <VERTEX_TYPE, GRAPH_TYPE> &aGraph, const typename Graph<VERTEX_TYPE>::VertexId &aSrc,
                     const typename Graph<VERTEX_TYPE>::VertexId &aDst, const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights,
                     int maxV) {
//                     std::unordered_map<int, int> &em) {
//    scanf("%d%d%d%d", &n, &m, &s, &t);
            n = maxV; //aGraph.getNumOfVertices();
            m = aGraph.getNumOfEdges();
            s = aSrc;
            t = aDst;

            edge_t *G = (edge_t *) malloc(m * sizeof(edge_t));
            int i = 0;
            for (auto &e : aGraph.getEdges()) {
                int u = e.src;
                int v = e.dst;
                int c = 1; //aWeights.at({u, v});

                G[i].f[uv] = G[i].f[vu] = 0;
                if (u < v) {
                    G[i].u = u;
                    G[i].v = v;
                    G[i].c[uv] = c;
                    G[i].c[vu] = 0;
                } else {
                    G[i].u = v;
                    G[i].v = u;
                    G[i].c[uv] = 0;
                    G[i].c[vu] = c;
                }
                i++;
            }

//    for (int i = 0; i < m; i++) {
//        int u, v, c;
//        scanf("%d%d%d", &u, &v, &c);
//
//        if (u == v) { //don't want any loops
//            continue;
//        }
//        G[i].f[uv] = G[i].f[vu] = 0;
//        if (u < v) {
//            G[i].u = u;
//            G[i].v = v;
//            G[i].c[uv] = c;
//            G[i].c[vu] = 0;
//        }
//        else {
//            G[i].u = v;
//            G[i].v = u;
//            G[i].c[uv] = 0;
//            G[i].c[vu] = c;
//        }
//    }

            m = merge_edges(G, n, m);

            for (int phase = 0; phase < n; phase++) {
                create_reserve_network(G);
                if (!clean_reserve_network()) {
                    break;
                }
                find_blocking_flow();
                augment_flow(G);
            }

            assert(!clean_reserve_network());

            long long max_flow = 0;
            for (int i = 0; i < m; i++) {
                if (G[i].u == s) {
                    max_flow += (long long) G[i].f[uv];
                }
                if (G[i].v == s) {
                    max_flow += (long long) G[i].f[vu];
                }
            }
//        printf("%lld\n", max_flow);

//        for (int i = 0; i < m; i++) {
//            if (G[i].f[uv] > 0) {
//                printf("%d %d %d\n", G[i].u, G[i].v, G[i].f[uv]);
//            }
//            if (G[i].f[vu] > 0) {
//                printf("%d %d %d\n", G[i].v, G[i].u, G[i].f[vu]);
//            }
//        }

            free(G);
            return max_flow;
        }
    };
}

#endif //FASPHEURISTIC_DINIC3RDPARTY_H
