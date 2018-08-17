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

#include "cutIntersections.h"

using namespace IMATI_STL;
bool _write_additional_files = false; // DEBUG STUFF


//////////////////////////////////////////////////////////////////
//
// Basic types to represent the intersection
// as a vector of segments
//
//////////////////////////////////////////////////////////////////

class ci_segment
{
public:
	Point p1, p2;						// The segment endpoints
	std::vector<Triangle *> triangles; // Triangles sharing this segment

	ci_segment() {}
	ci_segment(const Point& a, const Point& b) : p1(a), p2(b) {}

	void swapPoints() { Point tmp = p2; p2 = p1; p1 = tmp; }

	// TRUE iff the segments are geometrically equivalent (i.e. same endpoints)
	bool operator==(const ci_segment& s) const { return ((p1 == s.p1 && p2 == s.p2) || (p1 == s.p2 && p2 == s.p1)); }
	bool equal(const ci_segment& s) const { return operator==(s); }

	// Lexycographic comparison
	bool operator<(const ci_segment& s) const;

	// Splits the segment at p and returns the new subsegment
	ci_segment *split(const Point& p);

	// Simple split (assumes that p is not an endpoint)
	ci_segment *simpleSplit(const Point& p);

	// Splits the segment at all the ipoints.
	// The number 'n' of splits is returned.
	// The new parts can be retrieved by looking at the last 'n' segments
	// of any of this->triangles
	int splitAllIntersections(std::vector<ci_segment *>& s);

	void addIntersectionPoint(Point& ip) { ipoints.push_back(ip); }

	static bool segment_ptr_sort(ci_segment *a, ci_segment *b) { return ((*a) < (*b)); }

private:
	std::vector<Point> ipoints;		// Temporary list of points to split this segment

	static bool pointFarthestFromP2(const Point& q1, const Point& q2);
	void sortIPoints();
};

// Vector of the above-defined entity
typedef std::vector<ci_segment *> ci_segments;

// Information attached to each intersecting triangles
class triangleInfo 
{
public:
	ci_segments segments; // Segments on the triangle
	List triangles;		  // Other triangles which intersect this one

	triangleInfo(Triangle *t) // Constructor
	{ 
		triangles.joinTailList((List *)t->info);
		delete ((List *)t->info);
		t->info = this;
	}

	~triangleInfo() { }
};

// Vector of triangles
typedef std::vector<Triangle *> ci_Triangles;

// Access to information stored in the triangleInfo objects
inline ci_segments& ci_triangleSegments(Triangle *t) { return ((triangleInfo *)t->info)->segments; }
inline ci_segments *ci_triangleSegmentsPtr(Triangle *t) { return &((triangleInfo *)t->info)->segments; }
inline List *ci_triangleTriangles(Triangle *t) { return &((triangleInfo *)t->info)->triangles; }
inline void ci_destroyTriangleInfo(Triangle *t) { delete ((triangleInfo *)t->info); }

// Indexed segment - assumes the existence of a vector of points to be dereferenced
class ci_idxSegment
{
public:
	int i, j;
	ci_idxSegment(int a, int b) : i(a), j(b) {}

	bool operator<(const ci_idxSegment& s) const
	{
		if (i<s.i) return true; else if (i>s.i) return false; else return (j<s.j);
	}
};

// Vectors of indexed segments and points
typedef std::vector<ci_idxSegment> ci_idxSegments;
typedef std::vector<Point> ci_Points;


// Converts the vector of segments 's' into a pair made of a vector of points and
// a vector of indexed segments.
void ci_makeIndexedSegments(ci_segments *s, ci_Points& points3D, ci_idxSegments& si)
{
	ci_Points t;
	for (ci_segments::const_iterator i = s->begin(); i != s->end(); i++)
	{
		t.push_back((*i)->p1);
		t.push_back((*i)->p2);
	}
	ImatiSTLsort(t.begin(), t.end());

	for (ci_Points::const_iterator i = t.begin(); (i + 1) != t.end(); i++)
	if (((*i) != (*(i + 1)))) points3D.push_back((*i));

	if (points3D.empty() || points3D.back() != t.back()) points3D.push_back(t.back());

	unsigned int i, j;
	for (ci_segments::const_iterator k = s->begin(); k != s->end(); k++)
	{
		for (i = 0; i<points3D.size(); i++) if (points3D[i] == ((*k)->p1)) break;
		for (j = 0; j<points3D.size(); j++) if (points3D[j] == ((*k)->p2)) break;
		if (i != j) si.push_back(ci_idxSegment(i, j));
	}

	if (si.size()>1)
	{
		ImatiSTLsort(si.begin(), si.end());
		for (i = 0; i<(si.size() - 1);)
		if (!(si[i]<si[i + 1])) si.erase(si.begin() + i);
		else i++;
	}
}


bool ci_swappable(Edge *e)
{
	Vertex *v1 = e->v1;
	Vertex *v2 = e->t1->oppositeVertex(e);
	Vertex *v3 = e->v2;
	Vertex *v4 = e->t2->oppositeVertex(e);
	if (!v1->exactMisalignment(v2, v4)) return false;
	if (!v3->exactMisalignment(v2, v4)) return false;
	if (!v2->exactMisalignment(v1, v3)) return false;
	if (!v4->exactMisalignment(v1, v3)) return false;
	if (v1->exactSameSideOnPlane(v3, v2, v4)) return false;
	if (v2->exactSameSideOnPlane(v4, v1, v3)) return false;
	return true;
}

bool ci_swapOneVertex(Vertex *v, Edge *e0)
{
	List ve(e0), ve2;
	Edge *e, *el=NULL, *er=NULL;

	e = e0;
	do {
		e = e->leftTriangle(v)->prevEdge(e);
		if (IS_BIT(e, 5)) ve.appendHead(e);
	} while (IS_BIT(e, 5));
	el = e;

	e = e0;
	do {
		e = e->rightTriangle(v)->nextEdge(e);
		if (IS_BIT(e, 5)) ve.appendTail(e);
	} while (IS_BIT(e, 5));
	er = e;

	Vertex *vp = el->oppositeVertex(v);
	Vertex *vn = er->oppositeVertex(v);
	Vertex *vo = e0->oppositeVertex(v);

	if (vn->exactSameSideOnPlane(vo, vp, v))
	{
		Node *n;
		do
		{
			FOREACHVEEDGE((&(ve)), e, n)
			 if (ci_swappable(e)) { e->swap(true); ve2.appendHead(e); ve.removeCell(n); break; }
		} while (n!=NULL);

		if (ve.numels() == 0) return true; 
		else { FOREACHVEEDGE((&ve2), e, n) e->swap(true); return false; }
	}
	else return false;
}

// Swap all the edges within the polygon which are incident at e->v1.
Edge *ci_swapOneEdge(Edge *e)
{
	Edge *ne = e->t1->prevEdge(e);
	if (!IS_BIT(ne, 5) || ci_swapOneVertex(e->v1, ne)) { ne = e->t1->nextEdge(e); UNMARK_BIT(ne, 5); return ne; }
	else return NULL;
}


bool ci_swapStrip(std::vector<Triangle *>& strip)
{
	int i, ss = (int)strip.size(); // ss is the number of triangles in strip
	Edge *e;

	std::vector<Edge *> estrip; // estrip contains common edges of contiguous triangles in strip

	for (i = 0; i < ss - 1; i++)
	{
		e = strip[i]->commonEdge(strip[i + 1]);
		estrip.push_back(e);
		MARK_BIT(e, 5);
	}

	Vertex *v1 = strip[0]->oppositeVertex(estrip[0]); // First vertex
	Vertex *v2 = strip[ss - 1]->oppositeVertex(estrip[ss - 2]); // Last vertex
	
	// - reorient border edges
	for (i = 0; i < ss; i++)
	{
		e = strip[i]->e1; if (e->t1 != strip[i]) e->invert();
		e = strip[i]->e2; if (e->t1 != strip[i]) e->invert();
		e = strip[i]->e3; if (e->t1 != strip[i]) e->invert();
	}

	// - create bounding polygon
	List polygon;
	for (i = 0; i < ss-1; i++)
	{
		e = strip[i]->prevEdge(estrip[i]);
		if (!IS_BIT(e, 5)) polygon.appendTail(e);
	}
	e = strip[i]->nextEdge(estrip[i - 1]); polygon.appendTail(e);

	for (i = ss-1; i > 0; i--)
	{
		e = strip[i]->prevEdge(estrip[i-1]);
		if (!IS_BIT(e, 5)) polygon.appendTail(e);
	}
	e = strip[i]->nextEdge(estrip[i]); polygon.appendTail(e);


	// - swap the v1s of the oriented polygon edges
	Node *m, *n = polygon.head();
	Edge *ep, *en;
	int noswaps = 0;
	int gidx = 0;
	while (n != NULL)
	{
		if (noswaps++ > polygon.numels()) break;
		m = n->prev(); if (m == NULL) m = polygon.tail();
		e = (Edge *)n->data;
		ep = (Edge *)m->data;
		if (e->v1 != v1 && e->v1 != v2 && ep->hasVertex(e->v1))
		{
			if ((en = ci_swapOneEdge(e)) != NULL)
			{
				if (en->v2 != e->v2) en->invert();
				m->data = en;
				polygon.removeCell(n);
				n = m;
				noswaps = 0;
			}
		}
		n = n->next(); if (n == NULL) n = polygon.head();
	}

	// UNMARK EVERYTHING
	for (i = 0; i < ss-1; i++) UNMARK_BIT(estrip[i], 5);

	e = v1->getEdge(v2);
	if (e == NULL) return false;
	else TAG_SHARPEDGE(e);

	return true;
}


bool ci_insertSegment3D(ci_idxSegment& s, ci_Points& points3D)
{
	Vertex *v1 = (Vertex *)points3D[s.i].info;
	Vertex *v2 = (Vertex *)points3D[s.j].info;
	Edge *e = v1->getEdge(v2);

	// If edge exists just tag it for eventual splitting and return
	if (e != NULL) { TAG_SHARPEDGE(e); return true; }

	Edge *ne1, *ne2;
	std::vector<Triangle *> strip;

	// Locate first triangle in the strip within VT(v1)
	List *vt = v1->VT();
	Node *n;
	Triangle *t;
	FOREACHVTTRIANGLE(vt, t, n)
	{
		e = t->oppositeEdge(v1);
		if (Point::innerSegmentsCross(*v1, *v2, *e->v1, *e->v2))
		{
			if (IS_SHARPEDGE(e)) ImatiSTL::error("Intersecting segments detected\n");
			if (e->isOnBoundary()) ImatiSTL::error("Segment intersects boundary edge.\n");
			strip.push_back(t);
			break;
		}
	}
	delete vt;

	if (n == NULL) { ImatiSTL::warning("Could not start segment insertion.\n"); return false; }

	// Starting from first triangle, reconstruct the strip up to v2
	while (1)
	{
		t = e->oppositeTriangle(t);
		strip.push_back(t);
		ne1 = t->nextEdge(e); ne2 = t->prevEdge(e);
		if (ne1->commonVertex(ne2) == v2) break;
		if (Point::innerSegmentsCross(*v1, *v2, *ne1->v1, *ne1->v2)) e = ne1;
		else if (Point::innerSegmentsCross(*v1, *v2, *ne2->v1, *ne2->v2)) e = ne2;
		else ImatiSTL::error("Unexpected missing (or unproper) intersection.\n");
		if (IS_SHARPEDGE(e)) ImatiSTL::error("Intersecting segments detected (2)\n");
		if (e->isOnBoundary()) ImatiSTL::error("Segment intersects boundary edge (2)\n");
	}

	// Now strip contains the triangles from v1 to v2 that intersect the segment
	// We must simply retriangulate this strip with the constraint of having (v1-v2) as an edge
	return ci_swapStrip(strip);
	//return ci_retriangulateStrip(strip);
}


bool ci_exactPointInTriangle3D(Point& p, Triangle *t, Edge **ie)
{
	if (ci_pointInEdge(&p, t->e1)) { *ie = t->e1; return true; }
	if (ci_pointInEdge(&p, t->e2)) { *ie = t->e2; return true; }
	if (ci_pointInEdge(&p, t->e3)) { *ie = t->e3; return true; }
	return (ci_pointInInnerTriangle(&p, t));
}


void ci_locateTriangleToSplit(TriMesh *tin, Point& p, Triangle **ot, Edge **oe)
{
	Node *n;
	Triangle *t;

	*oe = NULL; *ot = NULL;
	FOREACHVTTRIANGLE((&(tin->T)), t, n) if (ci_exactPointInTriangle3D(p, t, oe))
	{
		if (*oe == NULL) *ot = t; else *ot = NULL;
		return;
	}

	ImatiSTL::warning("Point is not on any triangle!\n");
}


unsigned int ci_getCoplanarTriangles(Triangle *t0, std::vector<Triangle *>& itris, List *intersecting_triangles)
{
	Node *n;
	Triangle *t;
	Point t0v = t0->getVector();
	FOREACHVTTRIANGLE(intersecting_triangles, t, n)
	 if (!IS_SHARPEDGE(t) && ((t->getVector()&t0v).squaredLength() == 0)) itris.push_back(t);
	return itris.size();
}


bool ci_hasCoplanarTriangles(Triangle *t0)
{
	Node *n;
	Triangle *t;
	Point t0v = t0->getVector();
	List *intersecting_triangles = ci_triangleTriangles(t0);
	FOREACHVTTRIANGLE(intersecting_triangles, t, n)	if (((t->getVector()&t0v).squaredLength() == 0)) return true;
	return false;
}


inline bool ci_edgeIsOnTriangle(Edge *e, Triangle *t0)
{
	return (ci_pointInTriangle(e->v1, t0) && ci_pointInTriangle(e->v2, t0));
}


TriMesh *ci_triangulateTriangle(TriMesh *ntin, List *intersecting_triangles)
{
	Vertex *v1 = (Vertex *)ntin->V.head()->data;
	Vertex *v2 = (Vertex *)ntin->V.head()->next()->data;
	Vertex *v3 = (Vertex *)ntin->V.head()->next()->next()->data;
	Edge *e1 = (Edge *)ntin->E.head()->data;
	Edge *e2 = (Edge *)ntin->E.head()->next()->data;
	Edge *e3 = (Edge *)ntin->E.head()->next()->next()->data;

	Triangle *t0 = (Triangle *)ntin->T.head()->data;
	ci_segments *s = ci_triangleSegmentsPtr(t0);
	ci_idxSegments si;
	ci_Points points3D;
	ci_makeIndexedSegments(s, points3D, si);

	// 1) Split edges first. When inserted, each IPoint info-points to the vertex
	unsigned int i;
	std::vector<Point *> s1, s2, s3;
	for (i = 0; i<points3D.size(); i++)
	if (points3D[i] == (*v1)) points3D[i].info = v1;
	else if (points3D[i] == (*v2)) points3D[i].info = v2;
	else if (points3D[i] == (*v3)) points3D[i].info = v3;
	else if (Point::pointInInnerSegment(&(points3D[i]), e1->v1, e1->v2)) s1.push_back(&points3D[i]);
	else if (Point::pointInInnerSegment(&(points3D[i]), e2->v1, e2->v2)) s2.push_back(&points3D[i]);
	else if (Point::pointInInnerSegment(&(points3D[i]), e3->v1, e3->v2)) s3.push_back(&points3D[i]);

	// Here should sort s1, s2, s3 but this standard approach seems to
	// be too cumbersome for so few points. Just use brute force ...
	coord ad, md;
	unsigned int gi;
	Vertex *v;
	while (s1.size())
	{
		md = 0;
		for (i = 0; i<s1.size(); i++) if ((ad = s1[i]->squaredDistance(e1->v1)) >= md) { md = ad; gi = i; }
		v = ntin->splitEdge(e1, s1[gi]); s1[gi]->info = v; v->info = e1;
		s1.erase(s1.begin() + gi);
	}

	while (s2.size())
	{
		md = 0;
		for (i = 0; i<s2.size(); i++) if ((ad = s2[i]->squaredDistance(e2->v1)) >= md) { md = ad; gi = i; }
		v = ntin->splitEdge(e2, s2[gi]); s2[gi]->info = v; v->info = e2;
		s2.erase(s2.begin() + gi);
	}

	while (s3.size())
	{
		md = 0;
		for (i = 0; i<s3.size(); i++) if ((ad = s3[i]->squaredDistance(e3->v1)) >= md) { md = ad; gi = i; }
		v = ntin->splitEdge(e3, s3[gi]); s3[gi]->info = v; v->info = e3;
		s3.erase(s3.begin() + gi);
	}

	// 2) Insert inner points
	Triangle *t;
	Edge *e;
	for (i = 0; i<points3D.size(); i++) if (points3D[i].info == NULL)
	{
		ci_locateTriangleToSplit(ntin, points3D[i], &t, &e);
		if (t != NULL) { v = ntin->splitTriangle(t, &(points3D[i])); if (v) points3D[i].info = v; }
		else if (e != NULL) { v = ntin->splitEdge(e, &(points3D[i])); points3D[i].info = v; }
		else ImatiSTL::error("Triangle location failed. SHOULD NOT HAPPEN!\n");
	}

	// 3) Insert segments (and tag them as SHARP)
	for (i = 0; i<si.size(); i++)
	 if (!ci_insertSegment3D(si[i], points3D)) ImatiSTL::error("Segment insertion failed. SHOULD NOT HAPPEN!\n");

	Node *n;

	// Here we must check for coplanar triangles and add to them nonsharp edges as segments to guarantee
	// compatible triangulations.

	// 1) Determine coplanar triangles
	if (intersecting_triangles != NULL)
	{
		std::vector<Triangle *> itris;
		unsigned int nct = ci_getCoplanarTriangles(t0, itris, intersecting_triangles);
		ci_segment *seg;

		if (nct>0)
		{
			// 2) for each nonsharp edge -> for each coplanar triangle ct -> 
			FOREACHVEEDGE((&(ntin->E)), e, n) if (!IS_SHARPEDGE(e))
			{
				seg = NULL;
				for (i = 0; i < (int)nct; i++) if (ci_edgeIsOnTriangle(e, itris[i]))
				{
					if (seg == NULL) { seg = new ci_segment((*e->v1), (*e->v2)); TAG_SHARPEDGE(e); }
					ci_segments& cts = ci_triangleSegments(itris[i]); //  *((ci_segments *)itris[i]->info);
					cts.push_back(seg);
					seg->triangles.push_back(itris[i]);
				}
			}
		}
	}

	// Cut mesh along tagged segments

	Edge *ne;
	FOREACHVEEDGE((&(ntin->E)), e, n) if (IS_SHARPEDGE(e) && !e->isOnBoundary())
	{
		ne = ntin->newEdge(e->v1, e->v2);
		ne->t1 = e->t1; ne->t2 = NULL;
		e->t1 = NULL; ntin->E.appendHead(ne);
		ne->t1->replaceEdge(e, ne);
	}
	ntin->duplicateNonManifoldVertices();

	if (intersecting_triangles == NULL) ci_destroyTriangleInfo(t0);
	return ntin;
}



//////////////////////////////////////////////////////////////////
//
// Spatial subdivision
//
//////////////////////////////////////////////////////////////////

class cidi_cell : public di_cell
{
public:
//	int index, cost;

	cidi_cell() : di_cell() {}
	cidi_cell(TriMesh *tin) : di_cell(tin, true) {}

	void selectIntersections();
	inline cidi_cell *fork() { return (cidi_cell *)di_cell::fork(); }
	void computeIntersections(ci_segments&);
};


//////////////////////////////////////////////////////////////////
//
// Basic constructors for intersection points
//
//////////////////////////////////////////////////////////////////


// Computes the line-line intersection while assuming that the two lines actually intersect.
// Resulting intersection point is stored in res.
// This function temporarily switches to exact arithmetics to store the exact coordinates.
void ci_exactLineLineIntersection(const Point *a, const Point *b, const Point *c, const Point *d, Point& res)
{
	bool orat = ImatiSTL::isUsingRationals();
	ImatiSTL::useRationals(true);
	res = Point::lineLineIntersection((*a), (*b), (*c), (*d));
	ImatiSTL::useRationals(orat);
}

// Computes the intersection between the two coplanar segments (a-b) and (c-d).
// Segments are assumed to be exactly coplanar. Result is undetermined otherwise.
// Returns 0 if no intersection occurs, 1 if segments intersect at a point,
// 2 if segments overlap. In the latter case p1 and p2 are the endpoints
// of the intersection segment. If the intersection is just one
// point, then only p1 will be initialized.
int ci_segmentSegmentIntersection(const Point *a, const Point *b, const Point *c, const Point *d, Point& p1, Point& p2)
{
	int ni = 0;
	Point ints_points[2];
	if (Point::pointInSegment(a, c, d)) ints_points[ni++] = Point(a);
	if (Point::pointInSegment(b, c, d)) ints_points[ni++] = Point(b);
	if (ni == 2) { p1 = ints_points[0]; p2 = ints_points[1]; return 2; } // e1 subedge of e2
	if ((ni == 0 || (*(c)) != ints_points[0]) && Point::pointInSegment(c, a, b)) ints_points[ni++] = Point(c);
	if (ni == 2) { p1 = ints_points[0]; p2 = ints_points[1]; return 2; }
	if ((ni == 0 || (*(d)) != ints_points[0]) && Point::pointInSegment(d, a, b)) ints_points[ni++] = Point(d);
	if (ni == 2) { p1 = ints_points[0]; p2 = ints_points[1]; return 2; }

	if (ni>0) p1 = ints_points[0];
	if (ni>1) p2 = ints_points[1];

	if (ni) return ni;

	if (!a->exactSameSideOnPlane(b, c, d) && !c->exactSameSideOnPlane(d, a, b))
	{
		ci_exactLineLineIntersection(a, b, c, d, p1);
		return 1;
	}

	return 0;
}

int ci_edgeEdgeIntersection(const Edge *e1, const Edge *e2, Point& p1, Point& p2) { return ci_segmentSegmentIntersection(e1->v1, e1->v2, e2->v1, e2->v2, p1, p2); }


// Computes the intersection between an edge and the triangle.
// Returns 1 if intersection is made of one point only (=p1).
// Returns 2 if intersection is a segment (endpoints = p1 and p2).
// Returns 0 if no intersection occurs.
int ci_edgeTriangleIntersection(const Edge *e, const Triangle *t, Point& p1, Point& p2)
{
	Vertex *v1 = t->v1(), *v2 = t->v2(), *v3 = t->v3();
	Point p;

	coord d1 = v1->exactOrientation((v2), (v3), ((e->v1)));
	coord d2 = v1->exactOrientation((v2), (v3), ((e->v2)));

	if ((d1 > 0 && d2 > 0) || (d1 < 0 && d2 < 0)) return 0; // Edge completely over or below triangle

	if (d1 == 0 && d2 == 0) // Edge and triangle are coplanar
	{
		Point e11, e12, e21, e22, e31, e32;
		Vertex *cv = (e->hasVertex(v1)) ? (v1) : (e->hasVertex(v2) ? (v2) : (e->hasVertex(v3) ? (v3) : (NULL)));
		int e1i = ci_edgeEdgeIntersection(e, t->e1, e11, e12);
		int e2i = ci_edgeEdgeIntersection(e, t->e2, e21, e22);
		int e3i = ci_edgeEdgeIntersection(e, t->e3, e31, e32);
		if (e1i == 2) { p1 = e11; p2 = e12; return 2; }
		if (e2i == 2) { p1 = e21; p2 = e22; return 2; }
		if (e3i == 2) { p1 = e31; p2 = e32; return 2; }
		if (e1i == 1)
		{
			p1 = e11;
			if (e2i == 1 && e21 != p1) { p2 = e21; return 2; }
			if (e3i == 1 && e31 != p1) { p2 = e31; return 2; }
			if ((*(e->v1)) != p1 && ci_pointInInnerTriangle(e->v1, t)) { p2 = Point(e->v1); return 2; }
			if ((*(e->v2)) != p1 && ci_pointInInnerTriangle(e->v2, t)) { p2 = Point(e->v2); return 2; }
			if ((cv != NULL) && ((*cv) == p1)) return 0; else return 1;
		}
		if (e2i == 1)
		{
			p1 = e21;
			if (e3i == 1 && e31 != p1) { p2 = e31; return 2; }
			if ((*(e->v1)) != p1 && ci_pointInInnerTriangle(e->v1, t)) { p2 = Point(e->v1); return 2; }
			if ((*(e->v2)) != p1 && ci_pointInInnerTriangle(e->v2, t)) { p2 = Point(e->v2); return 2; }
			if ((cv != NULL) && ((*cv) == p1)) return 0; else return 1;
		}
		if (e3i == 1)
		{
			p1 = e31;
			if ((*(e->v1)) != p1 && ci_pointInInnerTriangle(e->v1, t)) { p2 = Point(e->v1); return 2; }
			if ((*(e->v2)) != p1 && ci_pointInInnerTriangle(e->v2, t)) { p2 = Point(e->v2); return 2; }
			if ((cv != NULL) && ((*cv) == p1)) return 0; else return 1;
		}
		if (ci_pointInInnerTriangle(e->v1, t) && ci_pointInInnerTriangle(e->v2, t))
		{
			p1 = Point(e->v1); p2 = Point(e->v2); return 2;
		}
	}

	// If not coplanar and have a common vertex cannot intersect
	if (e->hasVertex(v1) || e->hasVertex(v2) || e->hasVertex(v3)) return 0;

	coord e2 = v1->exactOrientation((v2), ((e->v1)), ((e->v2)));
	coord e3 = v2->exactOrientation((v3), ((e->v1)), ((e->v2)));
	coord e1 = v3->exactOrientation((v1), ((e->v1)), ((e->v2)));

	if (e1 == 0 && e2 == 0 && e3 == 0) return 0; // Should never happen ...

	if ((e1 >= 0 && e2 >= 0 && e3 >= 0) || (e1 <= 0 && e2 <= 0 && e3 <= 0))
	{
		if (d1 == 0) p1 = Point(e->v1);
		else if (d2 == 0) p1 = Point(e->v2);
		else if (e1 == 0 && e2 == 0) p1 = Point(v1);
		else if (e2 == 0 && e3 == 0) p1 = Point(v2);
		else if (e3 == 0 && e1 == 0) p1 = Point(v3);
		// Edge-edge intersection must be exactly commutative (maybe sorting edges lexycographically)
		else if (e1 == 0) ci_edgeEdgeIntersection(e, t->e1, p1, p2);
		else if (e2 == 0) ci_edgeEdgeIntersection(e, t->e2, p1, p2);
		else if (e3 == 0) ci_edgeEdgeIntersection(e, t->e3, p1, p2);
		else {
			bool orat = ImatiSTL::isUsingRationals();
			ImatiSTL::useRationals(true);
			p1 = Point::linePlaneIntersection(e->v1, e->v2, t->v1(), t->v2(), t->v3());
			ImatiSTL::useRationals(orat);
		}

		return 1;
	}

	return 0;
}


// Calculates the intersection between two triangles.
// If the triangles do not intersect returns 0, otherwise returns
// the intersecting point, segment or polygon.
// If only a point is shared, returns 1 and inits ints[0] with it.
// If a segment is shared, returns 2 and inits ints[0] and ints[].
// If t1 and t2 are coplanar, the intersection can be a convex
// polygon with up to six vertices. In this case returns the number
// N of vertices and inits ints[0] ... ints[N-1].
int ci_triangleTriangleIntersection(const Triangle *t1, const Triangle *t2, ci_segments& segments)
{
 Vertex *v11 = t1->v1(), *v12 = t1->v2(), *v13 = t1->v3();
 Vertex *v21 = t2->v1(), *v22 = t2->v2(), *v23 = t2->v3();

 coord mx = MIN(v11->x, MIN(v13->x, v12->x));
 if (v21->x < mx && v22->x < mx && v23->x < mx) return 0;
 mx = MAX(v11->x, MAX(v13->x, v12->x));
 if (v21->x > mx && v22->x > mx && v23->x > mx) return 0;
 mx = MIN(v11->y, MIN(v13->y, v12->y));
 if (v21->y < mx && v22->y < mx && v23->y < mx) return 0;
 mx = MAX(v11->y, MAX(v13->y, v12->y));
 if (v21->y > mx && v22->y > mx && v23->y > mx) return 0;
 mx = MIN(v11->z, MIN(v13->z, v12->z));
 if (v21->z < mx && v22->z < mx && v23->z < mx) return 0;
 mx = MAX(v11->z, MAX(v13->z, v12->z));
 if (v21->z > mx && v22->z > mx && v23->z > mx) return 0;

 // Case 1) All vertices of t1 are above/below plane of t2 and vice-versa
 coord o11 = v11->exactOrientation(v21, v22, v23);
 coord o12 = v12->exactOrientation(v21, v22, v23);
 coord o13 = v13->exactOrientation(v21, v22, v23);
 if ((o11>0 && o12>0 && o13>0) || (o11<0 && o12<0 && o13<0)) return 0;
 coord o21 = v21->exactOrientation(v11, v12, v13);
 coord o22 = v22->exactOrientation(v11, v12, v13);
 coord o23 = v23->exactOrientation(v11, v12, v13);
 if ((o21>0 && o22>0 && o23>0) || (o21<0 && o22<0 && o23<0)) return 0;

 int i, ni;
 Point p1, p2;
 // Case 2) t1 and t2 are coplanar
 if (o11==0 && o12==0 && o13==0) // Checking the o2i's should be redundant...
 {
  ci_segment polygon[6];
  Edge *t1e, *t2e;
  ni = 0;
  for (t1e = t1->e1; t1e != NULL; t1e = (t1e == t1->e3) ? (NULL) : (t1->nextEdge(t1e))) if (!t2->hasEdge(t1e))
  {
   i=ci_edgeTriangleIntersection(t1e, t2, p1, p2);
   if (i==1) polygon[ni++] = ci_segment(p1, p1);
   else if (i==2) polygon[ni++] = ci_segment(p1, p2);
  }
  for (t2e = t2->e1; t2e != NULL; t2e = (t2e == t2->e3) ? (NULL) : (t2->nextEdge(t2e))) if (!t1->hasEdge(t2e))
  {
   i=ci_edgeTriangleIntersection(t2e, t1, p1, p2);
   if (i == 1) polygon[ni++] = ci_segment(p1, p1);
   else if (i==2) polygon[ni++] = ci_segment(p1, p2);
  }

  for (i=0; i<ni; i++) segments.push_back(new ci_segment(polygon[i]));

  return ni;
 }

 // Case 3) Look for intersecting segment
 ni=0;
 Point ints[2];

 if (!t2->hasEdge(t1->e1))
 {
  i = ci_edgeTriangleIntersection(t1->e1, t2, p1, p2);
  if (i==2) {segments.push_back(new ci_segment(p1,p2)); return 1;}
  else if (i==1) ints[ni++]=p1;
 }

 if (!t2->hasEdge(t1->e2))
 {
  i = ci_edgeTriangleIntersection(t1->e2, t2, p1, p2);
  if (i==2) {segments.push_back(new ci_segment(p1,p2)); return 1;}
  else if (i==1) {if (ni==0) ints[ni++]=p1; else if (ni==1 && ints[0]!=p1) {segments.push_back(new ci_segment(ints[0],p1)); return 1;}}
 }

 if (!t2->hasEdge(t1->e3))
 {
  i = ci_edgeTriangleIntersection(t1->e3, t2, p1, p2);
  if (i==2) {segments.push_back(new ci_segment(p1,p2)); return 1;}
  else if (i==1) {if (ni==0) ints[ni++]=p1; else if (ni==1 && ints[0]!=p1) {segments.push_back(new ci_segment(ints[0],p1)); return 1;}}
 }

 if (!t1->hasEdge(t2->e1))
 {
  i = ci_edgeTriangleIntersection(t2->e1, t1, p1, p2);
  if (i == 2) { segments.push_back(new ci_segment(p1, p2)); return 1; }
  else if (i==1) {if (ni==0) ints[ni++]=p1; else if (ni==1 && ints[0]!=p1) {segments.push_back(new ci_segment(ints[0],p1)); return 1;}}
 }


 if (!t1->hasEdge(t2->e2))
 {
  i = ci_edgeTriangleIntersection(t2->e2, t1, p1, p2);
  if (i==2) {segments.push_back(new ci_segment(p1,p2)); return 1;}
  else if (i==1) {if (ni==0) ints[ni++]=p1; else if (ni==1 && ints[0]!=p1) {segments.push_back(new ci_segment(ints[0],p1)); return 1;}}
 }

 if (!t1->hasEdge(t2->e3))
 {
  i = ci_edgeTriangleIntersection(t2->e3, t1, p1, p2);
  if (i==2) {segments.push_back(new ci_segment(p1,p2)); return 1;}
  else if (i==1) {if (ni==0) ints[ni++]=p1; else if (ni==1 && ints[0]!=p1) {segments.push_back(new ci_segment(ints[0],p1)); return 1;}}
 }

 if (ni == 1)
 {
	 Vertex *cv = t1->commonVertex(t2);
	 if (cv != NULL) 
	 {
		 if ((*cv) == ints[0]) return 0;
		 else segments.push_back(new ci_segment(ints[0], *cv));
	 }
	 else segments.push_back(new ci_segment(ints[0], ints[0]));
 }

 return ni;
}


//////////////////////////////////////////////////////////////////
//
// Segments
//
//////////////////////////////////////////////////////////////////

bool ci_segment::operator<(const ci_segment& s) const
{
	const Point& a1 = (p1<p2) ? (p1) : (p2);
	const Point& a2 = (a1 == p1) ? (p2) : (p1);
	const Point& b1 = (s.p1<s.p2) ? (s.p1) : (s.p2);
	const Point& b2 = (b1 == s.p1) ? (s.p2) : (s.p1);
	if (a1<b1) return true; else if (a1 != b1) return false;
	if (a2<b2) return true; else return false;
}

ci_segment *ci_segment::split(const Point& p)
{
	if (p == p1 || p == p2) return this;
	else return simpleSplit(p);
}

ci_segment *ci_segment::simpleSplit(const Point& p)
{
	Point tmp = p2;
	p2 = p;
	ci_segment *ns = new ci_segment(p, tmp);
	Triangle *t;
	unsigned int i;
	for (i = 0; i<triangles.size(); i++)
	{
		t = triangles[i];
		ns->triangles.push_back(t);
		ci_segments *ts = ci_triangleSegmentsPtr(t); // (ci_segments *)(t->info);
		ts->push_back(ns);
	}

	return ns;
}


int ci_segment::splitAllIntersections(ci_segments& s)
{
	int ips;
	if ((ips = ipoints.size()) == 0) return 0;
	if (ips != 1) sortIPoints();
	ci_segment *ns;
	while (ipoints.size())
	{
		if ((ns = split(ipoints.back())) == this) ips--;
		else s.push_back(ns);
		ipoints.pop_back();
	}
	return ips;
}

bool ci_segment::pointFarthestFromP2(const Point& q1, const Point& q2)
{
	static Point s = q1;
	if (q2 == INFINITE_POINT) { s = q1; return false; }
	else return (q1.squaredDistance(&s)>q2.squaredDistance(&s));
}

void ci_segment::sortIPoints()
{
	pointFarthestFromP2(p2, Point(INFINITE_POINT));
	ImatiSTLsort(ipoints.begin(), ipoints.end(), &pointFarthestFromP2);
}


//////////////////////////////////////////////////////////////////
//
// Spatial subdivision
//
//////////////////////////////////////////////////////////////////


// Brute force all-with-all intersection test of the triangles in 'triangles'.
void cidi_cell::computeIntersections(ci_segments& segments)
{
	Triangle *t, *y;
	Node *n, *m;
	int i, ni;
	List *ts;

	for (n = triangles.head(); n != NULL; n = n->next())
	for (m = n->next(); m != NULL; m = m->next())
	{
		t = (Triangle *)n->data;
		y = (Triangle *)m->data; // For any pair (t,y) of triangles in the cell
		// The same triangle pair can be in different cells. The following avoids redoing the check.
		if (t->info == NULL || y->info == NULL || (((List *)t->info)->containsNode(y) == NULL))
		{
			if ((ni = ci_triangleTriangleIntersection(t, y, segments)) != 0)
			{
				MARK_VISIT(t); MARK_VISIT(y);
				ts = ((t->info != NULL) ? ((List *)t->info) : (new List)); ts->appendTail(y); t->info = ts;
				ts = ((y->info != NULL) ? ((List *)y->info) : (new List)); ts->appendTail(t); y->info = ts;
				for (i = segments.size() - ni; i<(int)segments.size(); i++)
				{
					segments[i]->triangles.push_back(t);
					segments[i]->triangles.push_back(y);
				}
			}
		}
	}
}

void ci_computeTriangleIntersections(Triangle *t, ci_segments& segments)
{
	List *ts = (List *)t->info;
	Node *n;
	Triangle *y;
	int i, ni;
	MARK_VISIT(t);

	FOREACHVTTRIANGLE(ts, y, n)	if (!IS_VISITED(y))
	{
		if ((ni = ci_triangleTriangleIntersection(t, y, segments)) != 0)
		{
			for (i = segments.size() - ni; i<(int)segments.size(); i++)
			{
				segments[i]->triangles.push_back(t);
				segments[i]->triangles.push_back(y);
			}
		}
	}
}

// Brute force all-with-all intersection test of the triangles in 'triangles'.
void cidi_cell::selectIntersections()
{
	di_cell::selectIntersections();

	Triangle *t;
	Node *n;

	for (n = triangles.head(); n!=NULL;)
	{
		t = (Triangle *)n->data;
		n = n->next();
		if (!IS_VISITED(t)) triangles.removeCell((n != NULL) ? (n->prev()) : triangles.tail());
	}
}

void ci_reassignSegments(ci_segments& segments)
{
 if (segments.size() == 0) return;
 Triangle *t;
 std::vector<Triangle *>::iterator j;
 ImatiSTLsort(segments.begin(), segments.end(), ci_segment::segment_ptr_sort);
 ci_segment *last = segments[0];
 bool same_as_prev;
 unsigned int i;
 for (i=0; i<segments.size();)
 {
  ci_segment *s = segments[i];
  same_as_prev=(i && s->equal(*last));
  if (!same_as_prev) last=s;
  for (j=s->triangles.begin(); j!=s->triangles.end(); j++)
  {
   t = (*j);
   ci_segments *ts = ci_triangleSegmentsPtr(t);
   ts->push_back(last);
   if (same_as_prev) last->triangles.push_back(t);
  }
  if (same_as_prev) segments.erase(segments.begin()+i);
  else i++;
 }
}


// Split all the segments on t0 that mutually intersect
bool lessThanOnX1s(const ci_segment *r, const ci_segment *s) { return (MIN(r->p1.x, r->p2.x)<MIN(s->p1.x, s->p2.x)); }
bool lessThanOnY1s(const ci_segment *r, const ci_segment *s) { return (MIN(r->p1.y, r->p2.y)<MIN(s->p1.y, s->p2.y)); }
bool lessThanOnZ1s(const ci_segment *r, const ci_segment *s) { return (MIN(r->p1.z, r->p2.z)<MIN(s->p1.z, s->p2.z)); }


void ci_computeSegmentIntersections(Triangle *t0, std::vector<Point *>& ips)
{
	unsigned int i, j;
	ci_segments& s = ci_triangleSegments(t0); // *((ci_segments *)t0->info);

	// Sort everything based on maximum extension

	Point mbb(DBL_MAX, DBL_MAX, DBL_MAX);
	Point Mbb(-DBL_MAX, -DBL_MAX, -DBL_MAX);
	for (i = 0; i < s.size(); i++)
	{
		if (s[i]->p1.x < mbb.x) mbb.x = s[i]->p1.x;
		if (s[i]->p1.y < mbb.y) mbb.y = s[i]->p1.y;
		if (s[i]->p1.z < mbb.z) mbb.z = s[i]->p1.z;
		if (s[i]->p2.x < mbb.x) mbb.x = s[i]->p2.x;
		if (s[i]->p2.y < mbb.y) mbb.y = s[i]->p2.y;
		if (s[i]->p2.z < mbb.z) mbb.z = s[i]->p2.z;
		if (s[i]->p1.x > Mbb.x) Mbb.x = s[i]->p1.x;
		if (s[i]->p1.y > Mbb.y) Mbb.y = s[i]->p1.y;
		if (s[i]->p1.z > Mbb.z) Mbb.z = s[i]->p1.z;
		if (s[i]->p2.x > Mbb.x) Mbb.x = s[i]->p2.x;
		if (s[i]->p2.y > Mbb.y) Mbb.y = s[i]->p2.y;
		if (s[i]->p2.z > Mbb.z) Mbb.z = s[i]->p2.z;
	}

	Mbb -= mbb;
	unsigned char wc = (Mbb.x > Mbb.y && Mbb.x > Mbb.z) ? (0) : ((Mbb.y > Mbb.x && Mbb.y > Mbb.z) ? (1) : (2));
	ImatiSTLsort(s.begin(), s.end(), (wc == 0) ? (&lessThanOnX1s) : ((wc == 1) ? (&lessThanOnY1s) : (&lessThanOnZ1s)));
	
	// Now look for intersections
	Point ip;
	for (i = 0; i < s.size(); i++)
	{
		for (j = i + 1; j < s.size(); j++)
		{
			if (MIN(s[j]->p1[wc], s[j]->p2[wc])>MAX(s[i]->p2[wc], s[i]->p1[wc])) break;
			if (Point::pointInInnerSegment(&s[i]->p1, &s[j]->p1, &s[j]->p2)) { ips.push_back(new Point(s[i]->p1)); ips.back()->info = s[j]; }
			if (Point::pointInInnerSegment(&s[i]->p2, &s[j]->p1, &s[j]->p2)) { ips.push_back(new Point(s[i]->p2)); ips.back()->info = s[j]; }
			if (Point::pointInInnerSegment(&s[j]->p1, &s[i]->p1, &s[i]->p2)) { ips.push_back(new Point(s[j]->p1)); ips.back()->info = s[i]; }
			if (Point::pointInInnerSegment(&s[j]->p2, &s[i]->p1, &s[i]->p2)) { ips.push_back(new Point(s[j]->p2)); ips.back()->info = s[i]; }
			if (Point::innerSegmentsCross(&(s[i]->p1), &(s[i]->p2), &(s[j]->p1), &(s[j]->p2)))
			{
				ci_exactLineLineIntersection(&(s[i]->p1), &(s[i]->p2), &(s[j]->p1), &(s[j]->p2), ip);
				ips.push_back(new Point(ip)); ips.back()->info = s[i];
				ips.push_back(new Point(ip)); ips.back()->info = s[j];
			}
		}
	}

}

// Split possible pairs of intersecting segments

int ci_splitIntersectingSegments(ci_segments& s, ci_Triangles& tris, clock_t timeout_secs)
{
	int totsplits = 0;
	int i = 0;
	int nt = (int)tris.size();
	std::vector<Point *> *trivecs = new std::vector<Point *>[nt];

	printf("Detecting intersecting segments...\n");
#pragma omp parallel for schedule(dynamic)
	for (i = 0; i < nt; i++)
	{
		printf("\r%d out of %d        ", totsplits, nt);
		ci_computeSegmentIntersections(tris[i], trivecs[i]);
		totsplits++;
		if (timeout_secs && !(totsplits % 1000)) ImatiSTL::exitOnTimeout();
	}

	totsplits = 0;
	unsigned int j;
	for (i = 0; i < nt; i++) for (j = 0; j < trivecs[i].size(); j++)
		((ci_segment *)trivecs[i][j]->info)->addIntersectionPoint(*(trivecs[i][j]));
	delete[] trivecs;

	printf("Splitting intersecting segments...\n");

	nt = 0;
	for (i = 0; i < (int)s.size(); i++)
	{
		printf("\r%d out of %d        ", nt, s.size());
		totsplits += s[i]->splitAllIntersections(s);
		nt++;
	}
	printf("\n");

	return totsplits;
}


int ci_saveVRML1(const char *filename, ci_segments& segments)
{
 FILE *fp;
 unsigned int i;
 if ((fp = fopen(filename,"w"))==NULL) {ImatiSTL::warning("Couldn't save graph!\n"); return 1;}

 fprintf(fp,"#VRML V1.0 ascii\nSeparator {\nMaterial { diffuseColor 1 0 0 }\nCoordinate3 {\npoint [\n");

 for (i=0; i<segments.size(); i++)
 {
  segments[i]->p1.printPoint(fp);
  segments[i]->p2.printPoint(fp);
 }
 fprintf(fp,"]\n}\nLineSet {\nnumVertices [\n");
 for (i=0; i<segments.size(); i++) fprintf(fp,"2,\n");
 fprintf(fp,"]\n}\n");

 fprintf(fp,"Material { diffuseColor 1 0 0 }\nDrawStyle { pointSize 4 }\n");
 fprintf(fp,"PointSet {}\n");

 fprintf(fp,"}\n");
 fclose(fp);

 return 0;
}


//int ci_saveCells(const char *filename, List& cells, TriMesh *tin)
//{
//	FILE *fp;
//	int i;
//	if ((fp = fopen(filename, "w")) == NULL) { ImatiSTL::warning("Couldn't save cells!\n"); return 1; }
//
//	fprintf(fp, "#VRML V1.0 ascii\nSeparator {\nMaterial { diffuseColor 1 0 0 }\n");
//
//	Node *n;
//	cidi_cell *c;
//	FOREACHNODE(cells, n)
//	{
//		c = (cidi_cell *)n->data;
//		Point center = (c->Mp + c->mp) / 2;
//		Point size = c->Mp - c->mp;
//		fprintf(fp, "Separator {\nTransform { translation %f %f %f }\n", center.x.toFloat(), center.y.toFloat(), center.z.toFloat());
//		fprintf(fp, "Cube {\n width %f\n height %f\n depth %f\n}\n}\n", size.x.toFloat(), size.y.toFloat(), size.z.toFloat());
//	}
//
//	fprintf(fp, "Separator {\n");
//	fprintf(fp, " Coordinate3 {\n  point [\n");
//
//	Vertex *v;
//	FOREACHVVVERTEX((&(tin->V)), v, n) fprintf(fp, "   %f %f %f,\n", v->x.toFloat(), v->y.toFloat(), v->z.toFloat());
//
//	fprintf(fp, "  ]\n }\n");
//
//	coord *ocds = new coord[tin->V.numels()];
//	i = 0; FOREACHVVVERTEX((&(tin->V)), v, n) { ocds[i] = v->x; v->x = i++; }
//
//	fprintf(fp, "Material {\n diffuseColor 0.6 0.6 0.6\n}\n");
//	fprintf(fp, " IndexedFaceSet {\n  coordIndex [\n");
//
//	Triangle *t;
//	FOREACHVTTRIANGLE((&(tin->T)), t, n)
//		fprintf(fp, "   %d, %d, %d, -1,\n", t->v1()->x.toInt(), t->v2()->x.toInt(), t->v3()->x.toInt());
//
//	i = 0; FOREACHVVVERTEX((&(tin->V)), v, n) v->x = ocds[i++];
//	delete[] ocds;
//
//	fprintf(fp, "  ]\n");
//
//	fprintf(fp, "}\n}\n}\n");
//	fclose(fp);
//
//	return 0;
//}


void ci_subdivideTriangles(TriMesh *tin, ci_Triangles& trilist, clock_t timeout_secs)
{
	std::vector<TriMesh *> triarray;  // Partial meshes of subdivided independent triangles (parallelizeable)
	std::vector<Triangle *> coplanar; // Triangles having coplanarities (require sequential processing)
	Triangle *t;
	int i;
	TriMesh *tri, *inters = new TriMesh;

	while (trilist.size()!=0)
	{
		t = trilist.back(); trilist.pop_back();
		if (ci_hasCoplanarTriangles(t)) coplanar.push_back(t);
		else  triarray.push_back((TriMesh *)tin->createSubMeshFromTriangle(t));
	}

	int ntris = triarray.size();
#pragma omp parallel for schedule(dynamic)
	for (i = 0; i<(int)triarray.size(); i++)
	{
		printf("\r%d remaining          ", ntris);
		ci_triangulateTriangle(triarray[i], NULL);
		if (timeout_secs && !(ntris % 1000)) ImatiSTL::exitOnTimeout();
		ntris--;
	}
	
	for (i = 0; i<(int)triarray.size(); i++) inters->moveMeshElements(triarray[i]);
	triarray.clear();

	for (i = 0; i<(int)coplanar.size(); i++)
	{
		printf("\r%d remaining          ", coplanar.size() - i);
		t = coplanar[i];
		tri = ci_triangulateTriangle((TriMesh *)tin->createSubMeshFromTriangle(t), ci_triangleTriangles(t)); TAG_SHARPEDGE(t);
		ci_destroyTriangleInfo(coplanar[i]);
		inters->moveMeshElements(tri);
	}
	coplanar.clear();

	tin->invertSelection();
	TriMesh *otg = (TriMesh *)tin->createSubMeshFromSelection();
	if (otg) inters->moveMeshElements(otg);
	if (_write_additional_files) inters->saveVRML1("splitMesh.wrl");
	tin->T.freeNodes(); tin->E.freeNodes(); tin->V.freeNodes();
	tin->moveMeshElements(inters);
}

////////////
//
// Cut Intersecting triangles
//
////////////

bool rm_cutIntersectingTriangles(TriMesh *tin, bool waf =false, int timeout_secs = 0)
{
 Triangle *t;
 Node *n;
 _write_additional_files = waf;

 // (1) remove all the degenerate triangles, by swaps if possible, otherwise by delete
 tin->deselectTriangles();
 if (tin->removeDegenerateTriangles()<0) tin->removeSelectedTriangles();

 // (2) organize triangles in a BSP to reduce complexity

 FOREACHVTTRIANGLE((&(tin->T)), t, n) t->info = NULL;

 const UINT16 tris_per_cell = 50; // This number is empirically the best tradeoff
 cidi_cell *c2, *c = new cidi_cell(tin);
 List cells;

 ImatiSTL::begin_progress();

 int i=0;
 List todo(c);
 while ((c = (cidi_cell *)todo.popHead()) != NULL)
 {
  if (i>DI_MAX_NUMBER_OF_CELLS || c->triangles.numels() <= tris_per_cell) cells.appendHead(c);
  else
  {
   if (!(i % 1000))
   {
		  ImatiSTL::report_progress(NULL);
		  if (timeout_secs) ImatiSTL::exitOnTimeout();
   }
   i++;
   c2 = c->fork();
   todo.appendTail(c); todo.appendTail(c2);
  }
 }
 //ci_saveCells("cells.wrl", cells, tin);

 // Now select the intersecting triangles
 i = 0; FOREACHNODE(cells, n)
 {
	 ((cidi_cell *)n->data)->selectIntersections();
	 if (!(i % 100))
	 {
		 ImatiSTL::report_progress("%d %% done   ", ((i)* 100) / cells.numels());
		 if (timeout_secs) ImatiSTL::exitOnTimeout();
	 }
	 i++;
 }

 // Dispose memory allocated for BSP cells
 while (cells.numels()) delete((cidi_cell *)cells.popHead());

 ci_Triangles trilist; // List of intersecting triangles
 FOREACHVTTRIANGLE((&(tin->T)), t, n) if (t->info != NULL) { trilist.push_back(t); UNMARK_VISIT(t); }

 // (3) Now actually compute the intersection segments
 ci_segments segments;
 for (i = 0; i < (int)trilist.size(); i++)
 {
	 t = trilist[i];
	 ci_computeTriangleIntersections(t, segments);
	  if (!(i % 100))
	  {
//if (segments.size()>200000) ImatiSTL::logToFileAndExit("Too many segments");
	   ImatiSTL::report_progress("%d %% done   ", ((i)* 100) / trilist.size());
	   if (timeout_secs) ImatiSTL::exitOnTimeout();
	  }
 }
 
 ImatiSTL::end_progress();
 ImatiSTL::printElapsedTime();

 // Associate a class triangleInfo to each triangle to retrieve intersecting triangles and segments
 for (i = 0; i < (int)trilist.size(); i++) new triangleInfo(trilist[i]);

 bool old_status = ImatiSTL::isUsingRationals();
 ImatiSTL::useRationals(true); // From now on exact coordinates are necessary
 ci_reassignSegments(segments); // Make triangles info-point to their segments (and remove duplications)

 if (trilist.size()) ImatiSTL::info("Intersection made of %d segments on %d triangles\n", segments.size(), trilist.size());
 else { ImatiSTL::info("No intersections detected.\n"); ImatiSTL::useRationals(old_status); return true; }

 // Subdivide coplanar segments that mutually intersect
 int sis = ci_splitIntersectingSegments(segments, trilist, timeout_secs);
 ImatiSTL::printElapsedTime();

 if (sis) ImatiSTL::info("%d intersecting segments have been split\n",sis);
 if (_write_additional_files) ci_saveVRML1("intersectionCurve.wrl", segments);

 // (4) Create copies of the intersecting triangles with segments inserted as constraints
 printf("Subdividing triangles\n");
 ci_subdivideTriangles(tin, trilist, timeout_secs);

 //ImatiSTL::printElapsedTime();
 //tin->saveOFF("cut.off");

 ImatiSTL::useRationals(old_status);

 return true;
}

bool TriMesh::cutIntersections(bool deb_info, int timeout_secs)
{
	return rm_cutIntersectingTriangles(this, deb_info, timeout_secs);
}







