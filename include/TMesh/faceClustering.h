/****************************************************************************
* TMesh                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
* Copyright(C) 2013: IMATI-GE / CNR                                         *
* All rights reserved.                                                      *
*                                                                           *
* This program is dual-licensed as follows:                                 *
*                                                                           *
* (1) You may use TMesh as free software; you can redistribute it and/or *
* modify it under the terms of the GNU General Public License as published  *
* by the Free Software Foundation; either version 3 of the License, or      *
* (at your option) any later version.                                       *
* In this case the program is distributed in the hope that it will be       *
* useful, but WITHOUT ANY WARRANTY; without even the implied warranty of    *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
* (2) You may use TMesh as part of a commercial software. In this case a *
* proper agreement must be reached with the Authors and with IMATI-GE/CNR   *
* based on a proper licensing contract.                                     *
*                                                                           *
****************************************************************************/

#ifndef FACECLUSTERING_H
#define FACECLUSTERING_H

#include "tin.h"
#include "heap.h"
#include "graph.h"

using namespace T_MESH;

class clusterEdge : public graphEdge
{
public:

	int index;
	double cost;

	clusterEdge(graphNode *n1, graphNode *n2, int i) : graphEdge(n1, n2) { index = i; cost = 0.0; }
};


class clusterHeap : abstractHeap
{
public:
	clusterEdge **edges;

	clusterHeap(int n, clusterEdge **edg) : abstractHeap(n)
	{
		positions = new int[n + 1]; edges = edg; for (int i = 0; i <= n; i++) positions[i] = 0;
	}

	~clusterHeap() { delete(positions); }

	void push(clusterEdge *e) { insert((void *)(e->index)); }
	clusterEdge *getFirst() { return (numels) ? (edges[(j_voidint)getHead()]) : (NULL); }
	clusterEdge *popHead() { return (numels) ? (edges[(j_voidint)removeHead()]) : (NULL); }
	int isEmpty() { return (numels == 0); }

	void remove(clusterEdge *e)
	{
		e->cost = -1; if (positions[e->index]) { upheap(positions[e->index]); removeHead(); }
	}

	void update(clusterEdge *e)
	{
		if (!positions[e->index]) push(e); else downheap(upheap(positions[e->index]));
	}

	int compare(const void *, const void *);
};


class clusterGraph : public Graph
{
	int curEdgeIndex, maxNumEdges;
	clusterEdge **ces;
	clusterHeap *ch;
	double(*costFunction)(const void *, const void *);

public:

	clusterGraph(int n, double(*cf)(const void *, const void *))
	{
		ces = new clusterEdge *[n];
		maxNumEdges = n;
		curEdgeIndex = 0;
		ch = new clusterHeap(n, ces);
		costFunction = cf;
	}
	~clusterGraph() { delete ces; delete ch; }

	clusterEdge *createEdge(graphNode *n1, graphNode *n2);

	clusterEdge *getFirstEdge();
	double getLowestCost() { clusterEdge *e = getFirstEdge(); return (e != NULL) ? (e->cost) : (DBL_MAX); }
	int collapseFirstEdge(void(*)(const void *, const void *) = NULL);
};

//////////////////////////////////////////////////////////////////////////
//
// Node of the cluster graph
//
//////////////////////////////////////////////////////////////////////////

class facecNode : public graphNode
{
 public:

 static double max_L_inf; // Maximum tolerable distance

 List triangles;        // All the triangles within the cluster
 List areas;            // Triangle areas
 Point sum_ctr;         // Weighted sum of barycenters
 SymMatrix3x3 Cov_v;    // Covariance matrix of cluster vertices
 double tot_area;       // Total area of the cluster

// Point fitplane_normal; // Parameters of the fitting primitive
// double fitplane_d;     // Parameters of the fitting primitive

 facecNode(Triangle *);	// Constructor (singleton)
 facecNode();	            // Constructor (empty cluster)
 ~facecNode() {areas.freeNodes();}

// bool initParameters();

 static void merge(const void *n1, const void *n2);
 static double edgeCostFunction(const void *n1, const void *n2);
 static double fittingPlaneCost(const void *, const void * =NULL);
};









//////////////////////////////////////////////////////////////////////////
//
// Cluster for Lloyd relaxation
//
//////////////////////////////////////////////////////////////////////////

class lloydCluster
{
 public:

 Point proxy_normal;
 double proxy_d;
 Triangle *seed;
 List triangles;
 int id;

 lloydCluster(facecNode *fc, int index);
 lloydCluster(Triangle *, int);     // Singleton
 lloydCluster(lloydCluster *, lloydCluster *, int);     // Merged cluster
 ~lloydCluster() {triangles.removeNodes();}     // Destructor

 void update(bool rt = true);	    // Update cluster pars
 double getCost(Triangle *); 	    // Return fitting residue
 double getTotalCost(); 	    // Return fitting residue
};


//////////////////////////////////////////////////////////////////////////
//
// Pair triangle-cluster to be prioritized in the queue
//
//////////////////////////////////////////////////////////////////////////

class clusteredTriangle
{
 public:
 Triangle *t;
 double cost;
 lloydCluster *l;

 clusteredTriangle(Triangle *f, lloydCluster *c) {t=f; l=c; cost=l->getCost(t);}
};


//////////////////////////////////////////////////////////////////////////
//
// Priority queue for triangle-cluster pairs
//
//////////////////////////////////////////////////////////////////////////

class triangleHeap : abstractHeap
{
 public:

 triangleHeap(int n) : abstractHeap(n) {}
 ~triangleHeap() {}

 inline void push(clusteredTriangle *e) {insert((void *)e);}
 inline clusteredTriangle *popHead() {return (clusteredTriangle *)removeHead();}
 inline clusteredTriangle *getHead() {return (clusteredTriangle *)abstractHeap::getHead();}
 inline int isEmpty() {return (numels==0);}

 int compare(const void *, const void *);
};

#endif // FACECLUSTERING_H
