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

#ifndef _EDGEHEAP_H
#define _EDGEHEAP_H

#include "edge.h"
#include "heap.h"

namespace T_MESH
{

//////////////////////////////////////////////////////
//
// Priority queue for edges.
// This allows logN updates in the position of elements.
//
// Usage:
// edgeHeap *eh = new edgeHeap(your_edge_list, your_cost_function);
// .. do whatever with your heap (it is already filled with edges in your list)
// delete eh
//
// This class uses the 'info' field of your edges.
// Attention! Do not delete edges before having deleted the heap (you may unlink them, but not delete)
//
//////////////////////////////////////////////////////

	class edgeHeap : public abstractHeap
	{
	protected:
		coord(*costFunction)(Edge *);
		Edge **edges;

	public:

		class EC_Cost
		{
		public:
			int index;
			coord cost;

			EC_Cost(Edge *e, int i) : index(i) {}
		};

		static inline const coord getEdgeCost(Edge *e) { return ((EC_Cost *)e->info)->cost; }
		static inline void setEdgeCost(Edge *e, coord d) { ((EC_Cost *)e->info)->cost = d; }
		static inline int getEdgeIndex(Edge *e) { return ((EC_Cost *)e->info)->index; }
		static inline void setEdgeIndex(Edge *e, int i) { ((EC_Cost *)e->info)->index = i; }
		static inline void createEdgeCost(Edge *e, int i) { e->info = new EC_Cost(e, i); }
		static inline void clearEdgeCost(Edge *e) { delete((EC_Cost *)e->info); e->info = NULL; }

		edgeHeap(List&, coord(*cf)(Edge *));	// constructor
		~edgeHeap();			// destructor

		// Insertion of a new element
		inline void push(Edge *e);

		// Removal of the first element
		inline Edge *popHead() { return edges[(int)removeHead()]; }

		// Emptiness check
		inline int isEmpty() { return (numels == 0); }

		void remove(Edge *);		// Remove one element
		void update(Edge *);		// Update one element
		int compare(const void *, const void *);	// Comparison for sorting
	};


//////////////////////////////////////////////////////
//
// Priority queue for edges.
// This can host a dynamic number of elements (reallocs memory if necessary).
//
// Usage:
// edgeHeap *eh = new edgeHeap(your_edge_list, your_cost_function);
// .. do whatever with your heap (it is already filled with edges in your list)
// delete eh
//
// This class uses the 'info' field of your edges.
// Attention! Do not delete edges before having deleted the heap (you may unlink them, but not delete)
//
//////////////////////////////////////////////////////

class dynamicEdgeHeap : abstractHeap
{
	coord(*costFunction)(Edge *);

public:

	dynamicEdgeHeap(int n, coord(*cf)(Edge *)) : abstractHeap(n) { costFunction = cf; };

	void push(Edge *);
	inline Edge *popHead() { return (Edge *)removeHead(); }
	inline int isEmpty() { return (numels == 0); }
	int compare(const void *, const void *);
	int getnum() const { return numels; }
};

} //namespace T_MESH

#endif //_EDGEHEAP_H

