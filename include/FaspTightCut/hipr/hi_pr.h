/* Maximum flow - highest lavel push-relabel algorithm */
/* COPYRIGHT C 1995, 2000 by IG Systems, Inc., igsys@eclipse.net */


/*
 * NOTE:
 *
 * This file contain source code from all needed files taken from original implementation:
 * http://www.avglab.com/andrew/soft.html
 *
 * In addition there is a:
 * - all is encapsulated in one class HIPR
 * - new function 'parse' which creates datastructure for HIPR from provided graph instead of stdin (as in original impl.)
 * - new function 'runHipr' to launch all needed steps (more or less as in original implementation but ending as soon as we have flow value)
 */


#ifndef HIPR_H
#define HIPR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../graph.h"
#include "../graphExt.h"

namespace Graph {

class HIPR {
#define DEXCESS_TYPE_LONG
#ifdef EXCESS_TYPE_LONG
            typedef unsigned long excessType;
#else
            typedef unsigned long long int excessType; /* change to double if not supported */
#endif

            typedef unsigned long cType;

            struct nodeSt;

            /* arc */
            struct arcSt {
                cType resCap;          /* residual capasity */
                nodeSt *head;           /* arc head */
                arcSt *rev;            /* reverse arc */
            };
            using arc  = arcSt;

            /* node */
            struct nodeSt {
                arc *first;           /* first outgoing arc */
                arc *current;         /* current outgoing arc */
                excessType excess;           /* excess at the node
				        change to double if needed */
                long d;                /* distance label */
                nodeSt *bNext;           /* next node in bucket */
                nodeSt *bPrev;           /* previous node in bucket */
            };
            using node = nodeSt;


            typedef /* bucket */
            struct bucketSt {
                node *firstActive;      /* first node with positive excess */
                node *firstInactive;    /* first node with zero excess */
            } bucket;

/*
#define GLOB_UPDT_FREQ 0.5
*/
#define GLOB_UPDT_FREQ 0.5
#define ALPHA 6
#define BETA 12

#define WHITE 0
#define GREY 1
#define BLACK 2
#define MAXLONG 2147483647
/* global variables */

            long n;                    /* number of nodes */
            long m;                    /* number of arcs */
            long nm;                   /* n + ALPHA * m */
            long nMin;                 /* smallest node id */
            node *nodes;               /* array of nodes */
            arc *arcs;                /* array of arcs */
            bucket *buckets;             /* array of buckets */
            cType *cap;                 /* array of capacities */
            node *source;              /* source node pointer */
            node *sink;                /* sink node pointer */
//node   **queue;              /* queue for BFS */
//node   **qHead, **qTail, **qLast;     /* queue pointers */
            long dMax;                 /* maximum label */
            long aMax;                 /* maximum actie node label */
            long aMin;                 /* minimum active node label */
            double flow;                 /* flow value */
            long pushCnt = 0;           /* number of pushes */
            long relabelCnt = 0;       /* number of relabels */
            long updateCnt = 0;       /* number of updates */
            long gapCnt = 0;           /* number of gaps */
            long gNodeCnt = 0;           /* number of nodes after gap */
//            float t, t2;                 /* for saving times */
            node *sentinelNode;        /* end of the node list marker */
            arc *stopA;                  /* used in forAllArcs */
            long workSinceUpdate = 0;      /* the number of arc scans since last update */
            float globUpdtFreq;          /* global update frequency */

/* macros */

#define forAllNodes(i) for ( i = nodes; i != sentinelNode; i++ )
#define forAllArcs(i, a) for (a = i->first, stopA = (i+1)->first; a != stopA; a++)

#define nNode(i) ( (i) - nodes + nMin )
#define nArc(a)  ( ( a == NULL )? -1 : (a) - arcs )

template <typename T> T min(const T &a, const T &b) {return ( ( (a) < (b) ) ? a : b ); }

/* FIFO queue for BFS macros */
/*
#define qInit() \
{\
  qHead = qTail = queue;\
}

#define qEmpty ( qHead == qTail )

#define qEnqueue(i) \
{\
  *qTail = i;\
  if ( qTail == qLast ) qTail = queue;\
  else qTail++;\
}

#define qDequeue(i) \
{\
  i = *qHead;\
  if ( qHead == qLast ) qHead = queue;\
  else qHead++;\
}
*/

/* 
   bucket macros:
   bucket's active node list is singly-linked
     operations aAdd, aRemove (from the front)
   bucket's inactive list is doubly-linked
     operations iAdd, iDelete (from arbitrary position)
*/

            long i_dist;

#define aAdd(l, i)\
{\
  i->bNext = l->firstActive;\
  l->firstActive = i;\
  i_dist = i->d;\
  if (i_dist < aMin)\
    aMin = i_dist;\
  if (i_dist > aMax)\
    aMax = i_dist;\
  if (dMax < aMax)\
    dMax = aMax;\
}

/* i must be the first element */
#define aRemove(l, i)\
{\
  l->firstActive = i->bNext;\
}

            node *i_next, *i_prev;
#define iAdd(l, i)\
{\
  i_next = l->firstInactive;\
  i->bNext = i_next;\
  i->bPrev = sentinelNode;\
  i_next->bPrev = i;\
  l->firstInactive = i;\
}

#define iDelete(l, i)\
{\
  i_next = i->bNext;\
  if (l->firstInactive == i) {\
    l->firstInactive = i_next;\
    i_next->bPrev = sentinelNode;\
  }\
  else {\
    i_prev = i->bPrev;\
    i_prev->bNext = i_next;\
    i_next->bPrev = i_prev;\
  }\
}

/* allocate datastructures, initialize related variables */

            int allocDS() {

                nm = ALPHA * n + m;
                /*
                queue = (node**) calloc ( n, sizeof (node*) );
                if ( queue == NULL ) return ( 1 );
                qLast = queue + n - 1;
                qInit();
                */
                buckets = (bucket *) calloc(n + 2, sizeof(bucket));
                if (buckets == NULL) return (1);

                sentinelNode = nodes + n;
                sentinelNode->first = arcs + 2 * m;

                return (0);

            } /* end of allocate */


            void init() {
                node *i;        /* current node */
                int overflowDetected;
                bucket *l;
                arc *a;
#ifdef EXCESS_TYPE_LONG
                double testExcess;
#endif
#ifndef OLD_INIT
                unsigned long delta;
#endif

                // initialize excesses

                forAllNodes(i) {
                    i->excess = 0;
                    i->current = i->first;
                    forAllArcs(i, a)a->resCap = cap[a - arcs];
                }

                for (l = buckets; l <= buckets + n - 1; l++) {
                    l->firstActive = sentinelNode;
                    l->firstInactive = sentinelNode;
                }

                overflowDetected = 0;
#ifdef EXCESS_TYPE_LONG
                testExcess = 0;
                forAllArcs(source,a) {
                  if (a->head != source) {
                    testExcess += a->resCap;
                  }
                }
                if (testExcess > MAXLONG) {
                  printf("c WARNING: excess overflow. See README for details.\nc\n");
                  overflowDetected = 1;
                }
#endif
#ifdef OLD_INIT
                source -> excess = MAXLONG;
#else
                if (overflowDetected) {
                    source->excess = MAXLONG;
                } else {
                    source->excess = 0;
                    forAllArcs(source, a) {
                        if (a->head != source) {
                            pushCnt++;
                            delta = a->resCap;
                            a->resCap -= delta;
                            (a->rev)->resCap += delta;
                            a->head->excess += delta;
                        }
                    }
                }

                /*  setup labels and buckets */
                l = buckets + 1;

                aMax = 0;
                aMin = n;

                forAllNodes(i) {
                    if (i == sink) {
                        i->d = 0;
                        iAdd(buckets, i);
                        continue;
                    }
                    if ((i == source) && (!overflowDetected)) {
                        i->d = n;
                    } else
                        i->d = 1;
                    if (i->excess > 0) {
                        /* put into active list */
                        aAdd(l, i);
                    } else { /* i -> excess == 0 */
                        /* put into inactive list */
                        if (i->d < n) iAdd(l, i);
                    }
                }
                dMax = 1;
#endif

                //  dMax = n-1;
                //  flow = 0.0;

            } /* end of init */

            void checkMax() {
                bucket *l;

                for (l = buckets + dMax + 1; l < buckets + n; l++) {
                    assert(l->firstActive == sentinelNode);
                    assert(l->firstInactive == sentinelNode);
                }
            }

/* global update via backward breadth first search from the sink */

            void globalUpdate() {

                node *i, *j;       /* node pointers */
                arc *a;           /* current arc pointers  */
                bucket *l, *jL;          /* bucket */
                long curDist, jD;
                long state;


                updateCnt++;

                /* initialization */

                forAllNodes(i)i->d = n;
                sink->d = 0;

                for (l = buckets; l <= buckets + dMax; l++) {
                    l->firstActive = sentinelNode;
                    l->firstInactive = sentinelNode;
                }

                dMax = aMax = 0;
                aMin = n;

                /* breadth first search */

                // add sink to bucket zero

                iAdd(buckets, sink);
                for (curDist = 0; 1; curDist++) {

                    state = 0;
                    l = buckets + curDist;
                    jD = curDist + 1;
                    jL = l + 1;
                    /*
                    jL -> firstActive   = sentinelNode;
                    jL -> firstInactive  = sentinelNode;
                    */

                    if ((l->firstActive == sentinelNode) &&
                        (l->firstInactive == sentinelNode))
                        break;

                    while (1) {

                        switch (state) {
                            case 0:
                                i = l->firstInactive;
                                state = 1;
                                break;
                            case 1:
                                i = i->bNext;
                                break;
                            case 2:
                                i = l->firstActive;
                                state = 3;
                                break;
                            case 3:
                                i = i->bNext;
                                break;
                            default:
                                assert(0);
                                break;
                        }

                        if (i == sentinelNode) {
                            if (state == 1) {
                                state = 2;
                                continue;
                            } else {
                                assert(state == 3);
                                break;
                            }
                        }

                        /* scanning arcs incident to node i */
                        forAllArcs(i, a) {
                            if (a->rev->resCap > 0) {
                                j = a->head;
                                if (j->d == n) {
                                    j->d = jD;
                                    j->current = j->first;
                                    if (jD > dMax) dMax = jD;

                                    if (j->excess > 0) {
                                        /* put into active list */
                                        aAdd(jL, j);
                                    } else {
                                        /* put into inactive list */
                                        iAdd(jL, j);
                                    }
                                }
                            }
                        } /* node i is scanned */
                    }
                }

            } /* end of global update */

/* second stage -- preflow to flow */
            void stageTwo()
/*
   do dsf in the reverse flow graph from nodes with excess
   cancel cycles if found
   return excess flow in topological order
*/

/*
   i->d is used for dfs labels 
   i->bNext is used for topological order list
   buckets[i-nodes]->firstActive is used for DSF tree
*/

            {
                node *i, *j, *tos, *bos, *restart, *r;
                arc *a;
                cType delta;

                /* deal with self-loops */
                forAllNodes(i) {
                    forAllArcs(i, a)
                        if (a->head == i) {
                            a->resCap = cap[a - arcs];
                        }
                }

                /* initialize */
                tos = bos = NULL;
                forAllNodes(i) {
                    i->d = WHITE;
                    //    buckets[i-nodes].firstActive = NULL;
                    buckets[i - nodes].firstActive = sentinelNode;
                    i->current = i->first;
                }

                /* eliminate flow cycles, topologicaly order vertices */
                forAllNodes(i)
                    if ((i->d == WHITE) && (i->excess > 0) &&
                        (i != source) && (i != sink)) {
                        r = i;
                        r->d = GREY;
                        do {
                            for (; i->current != (i + 1)->first; i->current++) {
                                a = i->current;
                                if ((cap[a - arcs] == 0) && (a->resCap > 0)) {
                                    j = a->head;
                                    if (j->d == WHITE) {
                                        /* start scanning j */
                                        j->d = GREY;
                                        buckets[j - nodes].firstActive = i;
                                        i = j;
                                        break;
                                    } else if (j->d == GREY) {
                                        /* find minimum flow on the cycle */
                                        delta = a->resCap;
                                        while (1) {
                                            delta = min (delta, j->current->resCap);
                                            if (j == i)
                                                break;
                                            else
                                                j = j->current->head;
                                        }

                                        /* remove delta flow units */
                                        j = i;
                                        while (1) {
                                            a = j->current;
                                            a->resCap -= delta;
                                            a->rev->resCap += delta;
                                            j = a->head;
                                            if (j == i)
                                                break;
                                        }

                                        /* backup DFS to the first saturated arc */
                                        restart = i;
                                        for (j = i->current->head; j != i; j = a->head) {
                                            a = j->current;
                                            if ((j->d == WHITE) || (a->resCap == 0)) {
                                                j->current->head->d = WHITE;
                                                if (j->d != WHITE)
                                                    restart = j;
                                            }
                                        }

                                        if (restart != i) {
                                            i = restart;
                                            i->current++;
                                            break;
                                        }
                                    }
                                }
                            }

                            if (i->current == (i + 1)->first) {
                                /* scan of i complete */
                                i->d = BLACK;
                                if (i != source) {
                                    if (bos == NULL) {
                                        bos = i;
                                        tos = i;
                                    } else {
                                        i->bNext = tos;
                                        tos = i;
                                    }
                                }

                                if (i != r) {
                                    i = buckets[i - nodes].firstActive;
                                    i->current++;
                                } else
                                    break;
                            }
                        } while (1);
                    }


                /* return excesses */
                /* note that sink is not on the stack */
                if (bos != NULL) {
                    for (i = tos; i != bos; i = i->bNext) {
                        a = i->first;
                        while (i->excess > 0) {
                            if ((cap[a - arcs] == 0) && (a->resCap > 0)) {
                                if (a->resCap < i->excess)
                                    delta = a->resCap;
                                else
                                    delta = i->excess;
                                a->resCap -= delta;
                                a->rev->resCap += delta;
                                i->excess -= delta;
                                a->head->excess += delta;
                            }
                            a++;
                        }
                    }
                    /* now do the bottom */
                    i = bos;
                    a = i->first;
                    while (i->excess > 0) {
                        if ((cap[a - arcs] == 0) && (a->resCap > 0)) {
                            if (a->resCap < i->excess)
                                delta = a->resCap;
                            else
                                delta = i->excess;
                            a->resCap -= delta;
                            a->rev->resCap += delta;
                            i->excess -= delta;
                            a->head->excess += delta;
                        }
                        a++;
                    }
                }
            }


/* gap relabeling */

            int gap(bucket *emptyB)
            {

                bucket *l;
                node *i;
                long r;           /* index of the bucket before l  */
                int cc;          /* cc = 1 if no nodes with positive excess before
		      the gap */

                gapCnt++;
                r = (emptyB - buckets) - 1;

                /* set labels of nodes beyond the gap to "infinity" */
                for (l = emptyB + 1; l <= buckets + dMax; l++) {
                    /* this does nothing for high level selection
                    for (i = l -> firstActive; i != sentinelNode; i = i -> bNext) {
                      i -> d = n;
                      gNodeCnt++;
                    }
                    l -> firstActive = sentinelNode;
                    */

                    for (i = l->firstInactive; i != sentinelNode; i = i->bNext) {
                        i->d = n;
                        gNodeCnt++;
                    }

                    l->firstInactive = sentinelNode;
                }

                cc = (aMin > r) ? 1 : 0;

                dMax = r;
                aMax = r;

                return (cc);

            }

/*--- relabelling node i */

            long relabel(node *i/* node to relabel */)
            {

                node *j;
                long minD;     /* minimum d of a node reachable from i */
                arc *minA;    /* an arc which leads to the node with minimal d */
                arc *a;

                assert(i->excess > 0);

                relabelCnt++;
                workSinceUpdate += BETA;

                i->d = minD = n;
                minA = NULL;

                /* find the minimum */
                forAllArcs(i, a) {
                    workSinceUpdate++;
                    if (a->resCap > 0) {
                        j = a->head;
                        if (j->d < minD) {
                            minD = j->d;
                            minA = a;
                        }
                    }
                }

                minD++;

                if (minD < n) {

                    i->d = minD;
                    i->current = minA;

                    if (dMax < minD) dMax = minD;

                } /* end of minD < n */

                return (minD);

            } /* end of relabel */


/* discharge: push flow out of i until i becomes inactive */

            void discharge(node *i)


            {

                node *j;                 /* sucsessor of i */
                long jD;                 /* d of the next bucket */
                bucket *lj;               /* j's bucket */
                bucket *l;                /* i's bucket */
                arc *a;                 /* current arc (i,j) */
                cType delta;
                arc *stopA;

                assert(i->excess > 0);
                assert(i != sink);
                do {

                    jD = i->d - 1;
                    l = buckets + i->d;

                    /* scanning arcs outgoing from  i  */
                    for (a = i->current, stopA = (i + 1)->first; a != stopA; a++) {
                        if (a->resCap > 0) {
                            j = a->head;

                            if (j->d == jD) {
                                pushCnt++;
                                if (a->resCap < i->excess)
                                    delta = a->resCap;
                                else
                                    delta = i->excess;
                                a->resCap -= delta;
                                a->rev->resCap += delta;

                                if (j != sink) {

                                    lj = buckets + jD;

                                    if (j->excess == 0) {
                                        /* remove j from inactive list */
                                        iDelete(lj, j);
                                        /* add j to active list */
                                        aAdd(lj, j);
                                    }
                                }

                                j->excess += delta;
                                i->excess -= delta;

                                if (i->excess == 0) break;

                            } /* j belongs to the next bucket */
                        } /* a  is not saturated */
                    } /* end of scanning arcs from  i */

                    if (a == stopA) {
                        /* i must be relabeled */
                        relabel(i);

                        if (i->d == n) break;
                        if ((l->firstActive == sentinelNode) &&
                            (l->firstInactive == sentinelNode)
                                )
                            gap(l);

                        if (i->d == n) break;
                    } else {
                        /* i no longer active */
                        i->current = a;
                        /* put i on inactive list */
                        iAdd(l, i);
                        break;
                    }
                } while (1);
            }


// go from higher to lower buckets, push flow
            void wave() {

                node *i;
                bucket *l;

                for (l = buckets + aMax; l > buckets; l--) {
                    for (i = l->firstActive; i != sentinelNode; i = l->firstActive) {
                        aRemove(l, i);

                        assert(i->excess > 0);
                        discharge(i);

                    }
                }
            }


/* first stage  -- maximum preflow*/

            void stageOne() {

                node *i;
                bucket *l;             /* current bucket */


#if defined(INIT_UPDATE) || defined(OLD_INIT) || defined(WAVE_INIT)
                globalUpdate ();
#endif

                workSinceUpdate = 0;

#ifdef WAVE_INIT
                wave();
#endif

                /* main loop */
                while (aMax >= aMin) {
                    l = buckets + aMax;
                    i = l->firstActive;

                    if (i == sentinelNode)
                        aMax--;
                    else {
                        aRemove(l, i);

                        assert(i->excess > 0);
                        discharge(i);

                        if (aMax < aMin)
                            break;

                        /* is it time for global update? */
                        if (workSinceUpdate * globUpdtFreq > nm) {
                            globalUpdate();
                            workSinceUpdate = 0;
                        }

                    }

                } /* end of the main loop */

                flow = sink->excess;

            }

    public:

        template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE>
        double runHipr(const Graph <VERTEX_TYPE> &aGraph,
                       const typename Graph<VERTEX_TYPE>::Vertex &aSrc,
                       const typename Graph<VERTEX_TYPE>::Vertex &aDst,
                       const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights,
                       std::vector<VERTEX_TYPE> &mapVertices) {
#if (defined(PRINT_FLOW) || defined(CHECK_SOLUTION))
                node *i;
                arc *a;
#endif

#ifdef PRINT_FLOW
                long ni, na;
#endif
#ifdef PRINT_CUT
                node *j;
#endif
                int cc;
#ifdef CHECK_SOLUTION
                excessType sum;
                bucket *l;
#endif



                globUpdtFreq = GLOB_UPDT_FREQ;

                node *nodesPtr;
                parse(&n, &m, &nodes, &arcs, &cap, &source, &sink, &nMin, &nodesPtr, aGraph, aSrc, aDst, aWeights, mapVertices);

                cc = allocDS();
                if (cc) {
                    fprintf(stderr, "Allocation error\n");
                    exit(1);
                }

                init();
                stageOne();

                free(nodesPtr);
                free(arcs);
                free(cap);
                free(buckets);

                return flow;
            }


        template<typename EDGE_PROP_TYPE, typename VERTEX_TYPE>
        int parse(long    *n_ad, long    *m_ad, node    **nodes_ad, arc     **arcs_ad, unsigned long    **cap_ad,
                  node    **source_ad, node    **sink_ad, long    *node_min_ad, node **nodesPtr, const Graph <VERTEX_TYPE> &aGraph, const typename Graph<VERTEX_TYPE>::Vertex &aSrc,
                  const typename Graph<VERTEX_TYPE>::Vertex &aDst, const Ext::EdgeProperties <VERTEX_TYPE, EDGE_PROP_TYPE> &aWeights,
                  std::vector<VERTEX_TYPE> &mapVertices) {

            long n = aGraph.getNumOfVertices();
            long m = aGraph.getNumOfEdges();

            /* allocating memory for  'nodes', 'arcs'  and internal arrays */
            node *nodes    = (node*) calloc ( n+2, sizeof(node) );
            arc *arcs     = (arc*)  calloc ( 2*m+1, sizeof(arc) );
            long *arc_tail = (long*) calloc ( 2*m,   sizeof(long) );
            long *arc_first= (long*) calloc ( n+2, sizeof(long) );
            unsigned long *acap     = (unsigned long*) calloc ( 2*m, sizeof(long) );
            /* arc_first [ 0 .. n+1 ] = 0 - initialized by calloc */
            arc *arc_current = arcs;


            long source = aSrc;
            long sink = aDst;
            long node_max = 0;
            long node_min = n;

            long pos_current=0;
            long no_alines=0;

            long head, tail, cap;
            node    *ndp;
            long arc_num;

            for (std::size_t i = 0; i < mapVertices.size(); ++i) mapVertices[i] = static_cast<VERTEX_TYPE>(-1);

            int newIdx = 0;
            for (auto &e : aGraph.getEdges()) {
                auto capacity = aWeights.at(e);
                if (mapVertices[e.src] == static_cast<VERTEX_TYPE>(-1)) mapVertices[e.src] = newIdx++;
                if (mapVertices[e.dst] == static_cast<VERTEX_TYPE>(-1)) mapVertices[e.dst] = newIdx++;
                tail = mapVertices[e.src];
                head = mapVertices[e.dst];
                cap = capacity;

                /* no of arcs incident to node i is stored in arc_first[i+1] */
                arc_first[tail + 1] ++;
                arc_first[head + 1] ++;

                /* storing information about the arc */
                arc_tail[pos_current]        = tail;
                arc_tail[pos_current+1]      = head;
                arc_current       -> head    = nodes + head;
                arc_current       -> resCap    = cap;
                arc_current       -> rev  = arc_current + 1;
                ( arc_current + 1 ) -> head    = nodes + tail;
                ( arc_current + 1 ) -> resCap    = 0;
                ( arc_current + 1 ) -> rev  = arc_current;

                /* searching minimumu and maximum node */
                if ( head < node_min ) node_min = head;
                if ( tail < node_min ) node_min = tail;
                if ( head > node_max ) node_max = head;
                if ( tail > node_max ) node_max = tail;

                no_alines   ++;
                arc_current += 2;
                pos_current += 2;

            }
            source = mapVertices[source];
            sink = mapVertices[sink];

            /* first arc from the first node */
            ( nodes + node_min ) -> first = arcs;

            /* before below loop arc_first[i+1] is the number of arcs outgoing from i;
               after this loop arc_first[i] is the position of the first
               outgoing from node i arcs after they would be ordered;
               this value is transformed to pointer and written to node.first[i]
               */

            for (long i = node_min + 1; i <= node_max + 1; i ++ )
            {
                arc_first[i]          += arc_first[i-1];
                ( nodes + i ) -> first = arcs + arc_first[i];
            }


            for (long  i = node_min; i < node_max; i ++ ) /* scanning all the nodes
                                            exept the last*/
            {

                long last = ((nodes + i + 1)->first) - arcs;
                /* arcs outgoing from i must be cited
                 from position arc_first[i] to the position
                 equal to initial value of arc_first[i+1]-1  */

                for (long arc_num = arc_first[i]; arc_num < last; arc_num++) {
                    long tail = arc_tail[arc_num];

                    while (tail != i)
                        /* the arc no  arc_num  is not in place because arc cited here
                           must go out from i;
                           we'll put it to its place and continue this process
                           until an arc in this position would go out from i */

                    {
                        long arc_new_num = arc_first[tail];
                        arc_current = arcs + arc_num;
                        arc *arc_new = arcs + arc_new_num;

                        /* arc_current must be cited in the position arc_new
                           swapping these arcs:                                 */

                        node *head_p = arc_new->head;
                        arc_new->head = arc_current->head;
                        arc_current->head = head_p;

                        cap = arc_new->resCap;
                        arc_new->resCap = arc_current->resCap;
                        arc_current->resCap = cap;

                        if (arc_new != arc_current->rev) {
                            arc *arc_tmp = arc_new->rev;
                            arc_new->rev = arc_current->rev;
                            arc_current->rev = arc_tmp;

                            (arc_current->rev)->rev = arc_current;
                            (arc_new->rev)->rev = arc_new;
                        }

                        arc_tail[arc_num] = arc_tail[arc_new_num];
                        arc_tail[arc_new_num] = tail;

                        /* we increase arc_first[tail]  */
                        arc_first[tail]++;

                        tail = arc_tail[arc_num];
                    }
                }
                /* all arcs outgoing from  i  are in place */
            }

            /* -----------------------  arcs are ordered  ------------------------- */

            /*----------- constructing lists ---------------*/


            for ( ndp = nodes + node_min; ndp <= nodes + node_max;  ndp ++ )
                ndp -> first = (arc*) NULL;

            for ( arc_current = arcs + (2*m-1); arc_current >= arcs; arc_current -- )
            {
                arc_num = arc_current - arcs;
                tail = arc_tail [arc_num];
                ndp = nodes + tail;
                /* avg
                arc_current -> next = ndp -> first;
                */
                ndp -> first = arc_current;
            }

            /* ----------- assigning output values ------------*/
            *m_ad = m;
            *n_ad = node_max - node_min + 1;
            *source_ad = nodes + source;
            *sink_ad   = nodes + sink;
            *node_min_ad = node_min;
            *nodes_ad = nodes + node_min;
            *arcs_ad = arcs;
            *cap_ad = acap;
            *nodesPtr = nodes;

            for ( arc_current = arcs, arc_num = 0;
                  arc_num < 2*m;
                  arc_current ++, arc_num ++
                    )
                acap [ arc_num ] = arc_current -> resCap;


            /* free internal memory */
            free ( arc_first ); free ( arc_tail );

            /* Thanks God! all is done */
            return (0);
        }
};

}


#endif