/****************************************************************************
* ImatiSTL                                                                  *
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
* (1) You may use ImatiSTL as free software; you can redistribute it and/or *
* modify it under the terms of the GNU General Public License as published  *
* by the Free Software Foundation; either version 3 of the License, or      *
* (at your option) any later version.                                       *
* In this case the program is distributed in the hope that it will be       *
* useful, but WITHOUT ANY WARRANTY; without even the implied warranty of    *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
* (2) You may use ImatiSTL as part of a commercial software. In this case a *
* proper agreement must be reached with the Authors and with IMATI-GE/CNR   *
* based on a proper licensing contract.                                     *
*                                                                           *
****************************************************************************/

#ifndef CUT_INTERSECTIONS_H
#define CUT_INTERSECTIONS_H

#include "imatistl.h"
#include "detectIntersections.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "jqsort.h"
#include <vector>
#include <algorithm>

#ifdef _HAVE_PPL
#include <ppl.h>

template<typename _Random_iterator, typename _Function>
inline void ImatiSTLsort(const _Random_iterator &_first, const _Random_iterator &_last, const _Function &_Pred)
{
	concurrency::parallel_sort(_first, _last, _Pred);
}

template<typename _Random_iterator>
inline void ImatiSTLsort(const _Random_iterator &_first, const _Random_iterator &_last)
{
	concurrency::parallel_sort(_first, _last);
}

#else
template<class _RanIt, class _Pr>
inline void ImatiSTLsort(_RanIt _first, _RanIt _last, _Pr _Pred)
{
	std::sort(_first, _last, _Pred);
}

template<class _RanIt>
inline void ImatiSTLsort(_RanIt _first, _RanIt _last)
{
	std::sort(_first, _last);
}
#endif

namespace IMATI_STL
{
	// Intersection tests for mesh elements
	inline bool ci_pointInEdge(const Point *p, const Edge *e) { return Point::pointInSegment(p, e->v1, e->v2); }
	inline bool ci_pointInInnerTriangle(const Point *p, const Triangle *t) { return Point::pointInInnerTriangle(p, t->v1(), t->v2(), t->v3()); }
	inline bool ci_pointInTriangle(const Point *p, const Triangle *t) { return Point::pointInTriangle(p, t->v1(), t->v2(), t->v3()); }
	inline bool ci_edgeIntersectsTriangle(const Edge *e, const Triangle *t) { return Point::segmentIntersectsTriangle(e->v1, e->v2, t->v1(), t->v2(), t->v3()); }
	inline bool ci_trianglesIntersect(const Triangle *t1, const Triangle *t2) { return t1->intersects(t2); }
} //namespace IMATI_STL

#endif // CUT_INTERSECTIONS_H
