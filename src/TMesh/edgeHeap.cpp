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

#include "edgeHeap.h"

namespace T_MESH
{

//////////////////////////////////////////////////////
//
// Implementation of class 'edgeHeap'
//
//////////////////////////////////////////////////////

edgeHeap::edgeHeap(List& edge_list, coord(*cf)(Edge *)) : abstractHeap(edge_list.numels())
{
	int n = edge_list.numels();
	costFunction = cf;
	positions = new int[n + 1];
	edges = new Edge *[n];	// Edge array for queue re-indexing
	Node *m;
	Edge *e;
	FOREACHVEEDGE((&edge_list), e, m) { e->info = NULL; push(e); }
}

edgeHeap::~edgeHeap()
{
	for (int i = 0; i < maxels; i++) clearEdgeCost(edges[i]);
	delete(positions); delete edges;
}

void edgeHeap::push(Edge *e)
{
	if (e->info == NULL)
	{
		edges[numels] = e;
		createEdgeCost(e, numels);
	}
	setEdgeCost(e, costFunction(e));
	insert((void *)(getEdgeIndex(e)));
}

void edgeHeap::remove(Edge *e)
{
	setEdgeCost(e, -1);
	int i = positions[getEdgeIndex(e)];
	if (i) { upheap(i); removeHead(); }
}

void edgeHeap::update(Edge *e)
{
	if (!positions[getEdgeIndex(e)]) push(e);
	else { setEdgeCost(e, costFunction(e)); downheap(upheap(positions[getEdgeIndex(e)])); }
}

int edgeHeap::compare(const void *e1, const void *e2)
{
	Edge *a = edges[(int)e1];
	Edge *b = edges[(int)e2];
	coord l1 = getEdgeCost(a);
	coord l2 = getEdgeCost(b);
	if (l1 < l2) return -1;
	if (l2 < l1) return 1;

	return 0;
}


//////////////////////////////////////////////////////
//
// Implementation of class 'dynamicEdgeHeap'
//
//////////////////////////////////////////////////////

void dynamicEdgeHeap::push(Edge *e)
{
	// If the heap size is exceeded, automatically double its size
	if (insert(e) == -1)
	{
		void **nheap = new void *[maxels * 2 + 1];
		for (int i = 0; i <= maxels; i++) nheap[i] = heap[i];
		delete(heap);
		heap = nheap;
		maxels *= 2;
		insert(e);
	}
}

int dynamicEdgeHeap::compare(const void *e1, const void *e2)
{
	coord l1 = costFunction((Edge *)e1);
	coord l2 = costFunction((Edge *)e2);
	if (l1 > l2) return -1;
	if (l2 > l1) return 1;

	return 0;
}

} //namespace T_MESH
