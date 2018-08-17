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



// Class ExtEdge
// This represents a non-manifold 'virtual' edge that points to all the mesh triangles that share it.
// This is not a proper mesh element.
// Triangles in 'et' are sorted so that they form a CCW fan around the oriented vector ep1-ep2.
// Triangles in 'et' are expected to have a boundary edge with endpoints that coincide with ep1-ep2.

class ExtEdge
{
public:
	Point ep1, ep2;	// The extedge's endpoints
	List et;		// List of virtually incident triangles

	ExtEdge(Edge *e) { ep1 = *e->v1; ep2 = *e->v2; }


	// Sort incident triangles based on pointers of their clusters (see below).
	// This is used to sort all the extedges based on their whole et relation.
	void resort()
	{
		if (et.numels() < 2) return;							// Nothing to sort
		void *t, *t0 = ((Triangle *)et.head()->data)->info;		// Clusters are stored in the info pointers of triangles
		Node *n;
		FOREACHNODE(et, n) if ((t = ((Triangle *)n->data)->info) < t0) t0 = t;

		while (t0 != ((Triangle *)et.head()->data)->info) // Cluster with  minimum pointer is the first
		{
			t = et.head()->data;
			et.removeCell(et.head());
			et.appendTail(t);
		}

		if (et.numels() < 3) return;
		if (((Triangle *)et.head()->next()->data)->info > ((Triangle *)et.tail()->data)->info) // Possibly reverse the list
		{
			t0 = et.head()->data;
			et.removeCell(et.head());
			void **etl = et.toArray();
			int ne = et.numels();
			et.removeNodes();
			for (int i = 0; i < ne; i++) et.appendHead(etl[i]);
			et.appendHead(t0);
			free(etl);

			Point ept = ep1;
			ep1 = ep2; ep2 = ept;
		}
	}

	void addUnsortedTriangle(Triangle *t) { et.appendTail(t); }
	void removeTriangle(Triangle *t) { et.removeNode(t); }


	// We assume a circular ordering as follows:
	// The triangle which is inserted first is the 'base'.
	// All the other triangles have a circular distance from the base (i.e. an angle 0 <= a <= 2PI).
	// This predicate is TRUE if the circular distance of 't' is strictly less than that of 'b'.
	// The distance is not actually computed - only exact geometric resoning is performed to
	// evaluate the predicate.

	bool isLessThan(Triangle *t, Triangle *b)
	{
		Point v = (ep1 - ep2);
		Vertex *ovt = t->oppositeVertex(getCorrEdge(t));
		Vertex *ovb = b->oppositeVertex(getCorrEdge(b));
		Point vt = ((*ovt) - ep2)&(v);
		Point vb = ((*ovb) - ep2)&(v);
		coord dotbt = vb*vt;
		coord crossbt = ((vt&vb)*v);
		if ((dotbt > 0) && (crossbt == 0)) return false; // If co-angular, return false

		Triangle *a = (Triangle *)et.head()->data; // This is the base triangle
		Vertex *ova = a->oppositeVertex(getCorrEdge(a));
		Point va = ((*ova) - ep2)&(v);

		coord dotat = va*vt;
		coord crossat = ((vt&va)*v);
		coord dotab = va*vb;
		coord crossab = ((vb&va)*v);


		int q_t, q_b;

		// Quadrant of t wrt a
		if ((dotat > 0) && (crossat >= 0)) q_t = 1;
		else if ((dotat >= 0) && (crossat < 0)) q_t = 4;
		else if ((dotat <= 0) && (crossat > 0)) q_t = 2;
		else if ((dotat < 0) && (crossat <= 0)) q_t = 3;

		// Quadrant of b wrt a
		if ((dotab > 0) && (crossab >= 0)) q_b = 1;
		else if ((dotab >= 0) && (crossab < 0)) q_b = 4;
		else if ((dotab <= 0) && (crossab > 0)) q_b = 2;
		else if ((dotab < 0) && (crossab <= 0)) q_b = 3;

		if (q_t == q_b) return (ovt->exactOrientation(ovb, &ep1, &ep2) <= 0);
		else return (q_t<q_b);
	}


	// Here we assume that the number of triangles is rather small (usually 2-4).
	// Therefore is not really worth using an efficient sorting algorithm.
	void sortTriangles()
	{
		if (et.numels() < 2) return;
		Node *n = et.head()->next()->next();
		Node *n0 = n;
		Triangle *cur, *prev;
		while (n != NULL)
		{
			cur = (Triangle *)n->data;
			prev = (Triangle *)n->prev()->data;
			if (isLessThan(cur, prev)) { n->data = prev; n->prev()->data = cur; if (n != n0) n = n->prev(); else n = n->next(); }
			else n = n->next();
		}
	}


	// Return the edge of 't' that matches ep1-ep2
	Edge *getCorrEdge(Triangle *t) const
	{
	 Vertex v1(ep1), v2(ep2);
	 Edge e(&v1, &v2);
     if (lexEdgeCompare(&e, t->e1)==0) return t->e1;
     else if (lexEdgeCompare(&e, t->e2)==0) return t->e2;
     else if (lexEdgeCompare(&e, t->e3)==0) return t->e3;
     else return NULL;
	}


	Triangle *getNextTriangle(Triangle *t0, bool& is_coherent) const
	{
		Node *n;
		Edge *e0, *ce;
		Triangle *t;

		e0 = getCorrEdge(t0);
		FOREACHVTTRIANGLE((&et), t, n) if (t == t0) break;
		if (n == NULL) ImatiSTL::error("ExtEdge::getNextTriangle: t0 is not in et.\n");
		Vertex *nv, *nv0 = e0->commonVertex(t0->nextEdge(e0));

		Node *g0 = n;
		bool cw = ((*nv0) == ep2);
		if (cw)
		{
			// Look forward
			n = (n->next() == NULL) ? (et.head()) : (n->next());
			t = (Triangle *)n->data;
			ce = getCorrEdge(t);
			nv = ce->commonVertex(t->prevEdge(ce));
			is_coherent = ((*nv) == (*nv0));
			return t;
		} else
		{
			// Look backward
			n = (n->prev() == NULL) ? (et.tail()) : (n->prev());
			t = (Triangle *)n->data;
			ce = getCorrEdge(t);
			nv = ce->commonVertex(t->prevEdge(ce));
			is_coherent = ((*nv) == (*nv0));
			return t;
		}

		return NULL;
	}


	Triangle *getPrevTriangle(Triangle *t0, bool& is_coherent) const
	{
		Node *n;
		Edge *e0, *ce;
		Triangle *t;

		e0 = getCorrEdge(t0);
		FOREACHVTTRIANGLE((&et), t, n) if (t == t0) break;
		if (n == NULL) ImatiSTL::error("ExtEdge::getNextTriangle: t0 is not in et.\n");
		Vertex *nv, *nv0 = e0->commonVertex(t0->nextEdge(e0));

		Node *g0 = n;
		bool cw = ((*nv0) == ep1);
		if (cw)
		{
			// Look forward
			n = (n->next() == NULL) ? (et.head()) : (n->next());
			t = (Triangle *)n->data;
			ce = getCorrEdge(t);
			nv = ce->commonVertex(t->prevEdge(ce));
			is_coherent = ((*nv) == (*nv0));
			return t;
		} else
		{
			// Look backward
			n = (n->prev() == NULL) ? (et.tail()) : (n->prev());
			t = (Triangle *)n->data;
			ce = getCorrEdge(t);
			nv = ce->commonVertex(t->prevEdge(ce));
			is_coherent = ((*nv) == (*nv0));
			return t;
		}

		return NULL;
	}


	//void printList()
	//{
	// Node *n;
	// Triangle *t;
	// FOREACHVTTRIANGLE((&et), t, n) printf("%p ", t->info);
	// printf("\n");
	//}

	static bool mergeEdges(Edge *ths, Edge *e)
	{
		if (ths->t1 && ths->t2) return 0;
		if (e->t1 && e->t2) return 0;
		Triangle *ot = e->getBoundaryTriangle();
		Vertex *ov1, *ov2;

		// Verify that vertices are geometrically coincident
		if (*e->v1 == *ths->v1) ov1 = e->v1; else if (*e->v2 == *ths->v1) ov1 = e->v2; else return false;
		if (*e->v1 == *ths->v2) ov2 = e->v1; else if (*e->v2 == *ths->v2) ov2 = e->v2; else return false;

		// Verify that triangles are consistently oriented
		if (ov1 == e->v1 && ((e->t1 != NULL && ths->t1 != NULL) || (e->t1 == NULL && ths->t1 == NULL))) return false;
		if (ov2 == e->v1 && ((e->t1 != NULL && ths->t1 == NULL) || (e->t1 == NULL && ths->t1 != NULL))) return false;

		Node *n;
		Edge *f;

		// Unify first vertex (if necessary)
		if (ov1 != ths->v1)
		{
			List *ve1 = ov1->VE();
			FOREACHVEEDGE(ve1, f, n) f->replaceVertex(ov1, ths->v1);
			delete(ve1);
			ov1->e0 = NULL;
		}
		// Unify second vertex (if necessary)
		if (ov2 != ths->v2)
		{
			List *ve2 = ov2->VE();
			FOREACHVEEDGE(ve2, f, n) f->replaceVertex(ov2, ths->v2);
			delete(ve2);
			ov2->e0 = NULL;
		}
		// Update TE, ET, and VE* relations
		ot->replaceEdge(e, ths);
		((ths->t1 == NULL) ? (ths->t1) : (ths->t2)) = ot;
		ths->v1->e0 = ths->v2->e0 = ths;

		// Unlink old edge
		e->v1 = e->v2 = NULL;

		return 1;
	}


	bool resolveIfSimple()
	{
		if (ep1 == ep2) ImatiSTL::error("resolveIfSimple(): degenerated extedge.");
		if (et.numels() != 2) return false;
		Triangle *t1 = (Triangle *)et.head()->data;
		Triangle *t2 = (Triangle *)et.head()->next()->data;
		if (t1 == NULL || t2 == NULL) ImatiSTL::error("resolveIfSimple(): null triangle.");
		if (t1->info != t2->info) return false;
		Edge *e1 = getCorrEdge(t1);
		Edge *e2 = getCorrEdge(t2);
		if (e1 == e2) ImatiSTL::error("resolveIfSimple(): same edge.");
		if (e1 == NULL || e2 == NULL) ImatiSTL::error("resolveIfSimple(): null edge.");
		if (!e1->isLinked() || !e2->isLinked()) ImatiSTL::error("resolveIfSimple(): unlinked edge.");
		bool ret = false;
		ret = mergeEdges(e1, e2);
		if (ret) { et.removeNodes(); e1->info = e2->info = NULL; }
		return ret;
	}

	bool resolveInCluster(void *c)
	{
		if (ep1 == ep2) ImatiSTL::error("resolveIfSimple(): degenerated extedge.");
		Triangle *t, *t1 = NULL, *t2 = NULL;
		Node *n;
		FOREACHVTTRIANGLE((&et), t, n) if (t->info == c) {
			if (t1 == NULL) t1 = t;
			else if (t2 == NULL) t2 = t;
			else return false; // More than two triangles from the same cluster meet at this edge ...
		}

		if (t1 == NULL || t2 == NULL) return false;
		Edge *e1 = getCorrEdge(t1);
		Edge *e2 = getCorrEdge(t2);
		if (e1 == e2) ImatiSTL::error("resolveInCluster(): same edge.");
		if (e1 == NULL || e2 == NULL) ImatiSTL::error("resolveInCluster(): null edge.");
		if (!e1->isLinked() || !e2->isLinked()) ImatiSTL::error("resolveInCluster(): unlinked edge.");
		bool ret = false;
		ret = mergeEdges(e1, e2);
		if (ret) { et.removeNode(t1); et.removeNode(t2); e1->info = e2->info = NULL; }
		return ret;
	}

	bool isIsolated() const { return (et.numels() == 0); }
	bool isOnBoundary() const { return (et.numels() == 1); }
	bool isTwoConnected() const { return (et.numels() == 2); }

	static bool lessThan(ExtEdge *a, ExtEdge *b)
	{
		List &xa = (a)->et, &xb = (b)->et;

		if (xa.numels() < xb.numels()) return true;
		if (xa.numels() > xb.numels()) return false;
		if (xa.numels() == 0) return false;

		Node *n = xa.head(), *m = xb.head();
		while (n != NULL && ((Triangle *)n->data)->info == ((Triangle *)m->data)->info) { n = n->next(); m = m->next(); }
		if (n == NULL) return false;
		if (((Triangle *)n->data)->info < ((Triangle *)m->data)->info) return true;
		else return false;
	}
};

typedef std::vector<ExtEdge *> ci_ExtEdges;





///////////////////////////////////////////////////////////////////////////////////////////
//
// FACE CLUSTERS AND B-REP
//
///////////////////////////////////////////////////////////////////////////////////////////


// Opposite triangle of a 2-connected edge (i.e. an edge having exactly two incident triangles).
// Returns NULL is edge is not 2-connected.
Triangle *coh_2connEdgeOppositeTriangle(Edge *e, Triangle *t, bool& is_coherent)
{
	is_coherent = true;
	if (!e->isOnBoundary()) return e->oppositeTriangle(t);
	ExtEdge *xe = ((ExtEdge *)e->info);
	if (!xe->isTwoConnected()) return NULL;
	return xe->getNextTriangle(t, is_coherent);
}


// coh_faceCluster
//
// A faceCluster is a maximal 2-edge-connected cluster of triangles (can be unorientable).
// Any pair of triangles in a faceCluster is 2-edge-connected and, due to maximality, if
// t1 is part of the faceCluster and t2 is 2-edge-connected with t1, then t2 is part of
// the faceCluster too.
// Each edge in a 2-edge-connection has exactly two incident triangles.
//
// A non-manifold mesh can be subdivided into a set of faceClusters (in the worst case
// each faceCluster is made of a single triangle). In such a subdivided mesh, we say that
// each triangle "belongs" to exactly one faceCluster.
// If 't' belong to the faceCluster 'C', we say that C = cluster(t).
//
// Let 'e' be a singular edge, let ET(e) be the set of all its incident triangles, and
// let t0(e) be one triangle in ET(e). The "triangle fan" of 'e' from t0 is the ordered list
// ET0(e, t0) whose elements are all and only the triangles in ET(e), t0 is the first
// element, and the triangles preserve their radial order around 'e'.
// On a subdivided mesh, we say that two triangle fans are homeomorphic if their triangles
// belong to the same clusters in the same order. I.e. F1 ~ F2 iff the two lists cluster(F1)
// and cluster(F2) are equal.
//
// On a subdivided mesh, we can define the concept of "clusterWall".
// A clusterWall W of a faceCluster C is a set of pairs W = {(T1, E1), ..., (Tn, En)}, where Ti is a
// triangle of C and Ei is one of Ti's edges, that satisfies the following conditions:
// 1) All the ET0(Ei,Ti) are homeomorphic
// 2) The Ei's form a single maximal 1-manifold chain of non-2-connected edges of C (either open or closed) 
// The 1-manifold chain of the Ei's is the "side" of the clusterWall.
// Two clusterWalls are 'matching' if they have the same side.
//
// A clusterEdge is an equivalence class of matching clusterWalls.
//
//
// A faceCluster is bounded by a number n >= 0 of 'walls', each corresponding to a clusterEdge.
// Different walls in a same cluster may correspond to the same clusterEdge.
//
// DATA STRUCTURE
// Relations are: C <-> W <-> E
//
// faceCluster: list of pointers to its triangles, list of pointers to its walls
// clusterWall: pointer to faceCluster, one pair of pointers (t,e), pointer to a clusterEdge
// clusterEdge: list of pointers to extedges, radially ordered list of pointers to incident clusterWalls
//
// If a cluster C is part of the outer hull, then it must be orientable and any of its clusterWalls
// must have at least a matching clusterWall (either from C or from another faceCluster).


class coh_clusterEdge; // Forward definition for the compiler. See below for actual definition.

class coh_faceCluster
{
public:
	List triangles;						// Triangles of the cluster
	List clusterEdges;					// coh_clusterEdges of the cluster
	List represTriangles;				// Representative triangles for clusterEdges (one each)
	coh_faceCluster *under_orientation;	// This cluster is being oriented according to another cluster
	j_voidint cluster_id;				// Unique identifier
	bool oriented;						// true if oriented
	bool selected;						// true if triangles are selected
	bool linked;						// linked if currently active

	// Constructor
	// Starting from a seed triangle 't0' traverses the mesh across 2-edge-connections and
	// define the set of triangles constituting this cluster.
	// clusterEdges and represTriangles are not defined at this stage.
	coh_faceCluster(Triangle *t0, TriMesh *tin)
	{
		Triangle *t, *s;
		triangles.appendTail(t0); MARK_BIT(t0, 5);
		Node *n = triangles.head();
		Edge *e;
		bool is_coherent;
		static j_voidint cid = 0;

		cluster_id = ++cid;

		oriented = true;
		selected = false;
		linked = true;
		under_orientation = NULL;

		while (n != NULL)
		{
			t = (Triangle *)n->data;
			t->info = (void *)cluster_id;
			FOREACHTRIANGLEEDGE(t, e) if ((s = coh_2connEdgeOppositeTriangle(e, t, is_coherent))!=NULL)
			{
				if (!IS_BIT(s, 5))
				{ 
					MARK_BIT(s, 5);
					triangles.appendTail(s);
					if (!is_coherent) tin->flipNormals(s);
				}
				else if (!is_coherent) oriented=false;
			}
			n = n->next();
		}
		FOREACHVTTRIANGLE((&triangles), t, n) UNMARK_BIT(t, 5);
	}

	~coh_faceCluster() { unlink(); }

	void assignTriangleInfos()
	{
		Node *n;
		Triangle *t;
		FOREACHVTTRIANGLE((&triangles), t, n) t->info = this;
	}

	// Returns the represTriangle which corresponds to 'e' in the list
	Triangle *getCorrRepresTriangle(coh_clusterEdge *e)
	{
		Node *n = clusterEdges.head(), *m = represTriangles.head();
		while (n != NULL) if (n->data == e) break; else { n = n->next(); m = m->next(); }
		if (n == NULL) return NULL;
		return ((Triangle *)m->data);
	}

	void setUnderOrientation(coh_faceCluster *u) { under_orientation = u; }
	bool underOrientation(coh_faceCluster *u) const { return (under_orientation == u); }

	void unlink();
	bool isLinked() const { return linked; }

	void addEdge(coh_clusterEdge *e, Triangle *rt) { clusterEdges.appendTail(e); represTriangles.appendTail(rt); }

	bool isOnBoundary() const;

	bool isOriented() const { return oriented; }

	void selectTriangles()
	{
		Node *n;
		Triangle *t;
		FOREACHVTTRIANGLE((&triangles), t, n) MARK_VISIT(t);
		selected = true;
	}

	void deselectTriangles()
	{
		Node *n;
		Triangle *t;
		FOREACHVTTRIANGLE((&triangles), t, n) UNMARK_VISIT(t);
		selected = false;
	}

	bool isSelected() const { return selected; }

	List *getMeshEdges()
	{
		Node *n;
		Triangle *t;
		Edge *e;
		List *tredges = new List;
		FOREACHVTTRIANGLE((&triangles), t, n)
		{
			FOREACHTRIANGLEEDGE(t, e) UNMARK_BIT(e, 6);
		}
		FOREACHVTTRIANGLE((&triangles), t, n)
		{
			FOREACHTRIANGLEEDGE(t, e) if (!IS_BIT(e, 6)) { MARK_BIT(e, 6); tredges->appendTail(e); }
		}
		FOREACHVEEDGE((tredges), e, n) UNMARK_BIT(e, 6);

		return tredges;
	}

	bool hasMarked4vertex() const
	{
		Node *n;
		Triangle *t;
		FOREACHVTTRIANGLE((&triangles), t, n)
		{
			if (IS_BIT(t->v1(), 4) || IS_BIT(t->v2(), 4) || IS_BIT(t->v3(), 4)) return true;
		}
		return false;
	}

	void flip()
	{
		Node *n;
		Triangle *t;
		Vertex *tmp;
		Edge *e;
		List *tredges = getMeshEdges();
		FOREACHVTTRIANGLE((&triangles), t, n) t->invert();
		FOREACHVEEDGE((tredges), e, n) { tmp = e->v1; e->v1 = e->v2; e->v2 = tmp; }
		delete tredges;
	}

	void mergeEdges()
	{
		Node *n;
		Edge *e;
		List *tredges = getMeshEdges();
		FOREACHVEEDGE((tredges), e, n) if (e->isLinked() && e->isOnBoundary()) ((ExtEdge *)e->info)->resolveInCluster(this);
		delete tredges;
	}

	void saveSTL()
	{
		char fname[2048];
		sprintf(fname, "cluster_%d.stl", triangles.numels());
		FILE *fp = fopen(fname, "w");
		fprintf(fp,"solid %s", fname);
		Node *n;
		Triangle *t;
		FOREACHVTTRIANGLE((&triangles), t, n)
		{
			Point nor = t->getNormal();
			fprintf(fp, "facet normal %f %f %f\n", TMESH_TO_DOUBLE(nor.x), TMESH_TO_DOUBLE(nor.y), TMESH_TO_DOUBLE(nor.z));
			fprintf(fp, "outer loop\n");
			fprintf(fp, "vertex %f %f %f\n", TMESH_TO_DOUBLE(t->v1()->x), TMESH_TO_DOUBLE(t->v1()->y), TMESH_TO_DOUBLE(t->v1()->z));
			fprintf(fp, "vertex %f %f %f\n", TMESH_TO_DOUBLE(t->v2()->x), TMESH_TO_DOUBLE(t->v2()->y), TMESH_TO_DOUBLE(t->v2()->z));
			fprintf(fp, "vertex %f %f %f\n", TMESH_TO_DOUBLE(t->v3()->x), TMESH_TO_DOUBLE(t->v3()->y), TMESH_TO_DOUBLE(t->v3()->z));
			fprintf(fp, "endloop\n");
			fprintf(fp, "endfacet\n");
		}
		fprintf(fp, "endsolid %s", fname);
		fclose(fp);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// TEMPORARY CODE TO PRODUCE NICE COLORED IMAGES OF THE CLUSTERED MESH
	//
	/////////////////////////////////////////////////////////////////////////////////////////////////

	void printVRMLCluster(FILE *fp)
	{
		fprintf(fp, "Material {\n");
		double r, g, b, l;
		r = ((double)rand()) / RAND_MAX;
		g = ((double)rand()) / RAND_MAX;
		b = ((double)rand()) / RAND_MAX;
		l = sqrt(r*r + g*g + b*b);
		r /= l; g /= l; b /= l;
		//if (isOnBoundary()) { r = 0.3; g = 1; b = 0.3; }
		//else { r = 0.6; g = 0.6; b = 0.6; }
		fprintf(fp, " diffuseColor %f %f %f\n", r, g, b);
		fprintf(fp, "}\n");

		fprintf(fp, "IndexedFaceSet {\n coordIndex [\n");
		Triangle *t;
		Node *n;
		FOREACHVTTRIANGLE((&(triangles)), t, n)
			fprintf(fp, "%d, %d, %d, -1,\n", (int)t->v1()->info, (int)t->v2()->info, (int)t->v3()->info);

		fprintf(fp, "]\n}\n");
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// END OF TEMPORARY CODE TO PRODUCE NICE COLORED IMAGES OF THE CLUSTERED MESH
	//
	/////////////////////////////////////////////////////////////////////////////////////////////////
};


class coh_clusterEdge
{
	List extedges;
	List incidentClusters;
	List represTriangles; // Each cluster has a representative triangle (to distinguish in case of multiplicity)
	bool linked;
public:

	coh_clusterEdge(ExtEdge *e)
	{
		addEdge(e);
		linked = true;
	}

	// Total incidences: each cluster may count more than once
	int getNumClusters() const { return incidentClusters.numels(); }

	void addEdge(ExtEdge *e) { extedges.appendTail(e); }
	
	void removeCluster(coh_faceCluster *f)
	{
		Node *n = incidentClusters.head(), *m = represTriangles.head();
		while (n != NULL) if (n->data == f)
		{
			n = n->next(); m = m->next();
			incidentClusters.removeCell((n == NULL) ? (incidentClusters.tail()) : (n->prev()));
			represTriangles.removeCell((m == NULL) ? (represTriangles.tail()) : (m->prev()));
		}
		else { n = n->next(); m = m->next(); }
	}

	bool isIsolated() const { return (incidentClusters.numels() == 0); }

	void initRelations()
	{
		ExtEdge *e0 = (ExtEdge *)extedges.head()->data;
		Node *n;
		Triangle *t;
		List tmpClusters;
		FOREACHVTTRIANGLE((&(e0->et)), t, n)
		{
			coh_faceCluster *fc = (coh_faceCluster *)t->info;
			if (!incidentClusters.containsNode(fc))
			{
				tmpClusters.appendTail(fc);
				represTriangles.appendTail(t);
				fc->addEdge(this, t);
			}
		}
		incidentClusters.joinTailList(&tmpClusters);
	}

	bool isOnBoundary() const
	{
		return ((ExtEdge *)extedges.head()->data)->isOnBoundary();
	}


	coh_faceCluster *getNextCluster(coh_faceCluster *f, Triangle *t0, bool& is_coherent)
	{
		ExtEdge *e0 = (ExtEdge *)extedges.head()->data;
		Node *n;
		Triangle *t;
		FOREACHVTTRIANGLE((&(e0->et)), t, n) if (t == t0) break;
		if (n == NULL) return NULL;
		Triangle *t2 = e0->getNextTriangle(t, is_coherent);
		if (t2 == NULL) return NULL;
		return ((coh_faceCluster *)t2->info);
	}

	coh_faceCluster *getOppositeCluster(coh_faceCluster *f)
	{
		if (incidentClusters.numels() != 2) return NULL;
		if (incidentClusters.head()->data == f) return ((coh_faceCluster *)incidentClusters.head()->next()->data);
		if (incidentClusters.head()->next()->data == f) return ((coh_faceCluster *)incidentClusters.head()->data);
		return NULL;
	}

	void unlink()
	{
		incidentClusters.removeNodes();
		linked = false;
	}

	bool isLinked() const { return linked; }

	List *getIncidentClusters() { return &incidentClusters; }
};

bool coh_faceCluster::isOnBoundary() const
{
	Node *n;
	FOREACHNODE(clusterEdges, n) if (((coh_clusterEdge *)n->data)->isOnBoundary()) return true;
	return false;
}

void coh_faceCluster::unlink()
{
	coh_clusterEdge *e;
	Node *n;
	FOREACHNODE(clusterEdges, n)
	{
		e = (coh_clusterEdge *)n->data;
		e->removeCluster(this);
		if (e->getNumClusters() == 0) e->unlink();
	}
	clusterEdges.removeNodes();
	represTriangles.removeNodes();
	linked = false;
}


class coh_clusterRep
{
public:
	TriMesh *tin;			// Input mesh pieces
	List dangling_tins;	    // List of output dangling mesh pieces
	List clusters;			// Clusters (2-edge-connected mesh pieaces)
	List edges;				// Edges bounding clusters

	coh_clusterRep(TriMesh *_t, ci_ExtEdges& extedges)
	{
		tin = _t;
		Triangle *t;
		Node *n;
		FOREACHVTTRIANGLE((&(tin->T)), t, n) { t->info = NULL; UNMARK_BIT(t, 5); }

		printf("making clusters..");
		FOREACHVTTRIANGLE((&(tin->T)), t, n) if (t->info == NULL) clusters.appendTail(new coh_faceCluster(t, tin));
		printf("%d clusters made.\n", clusters.numels());

		// Make the 'et' start from the triangle corresp. to smallest cluster (id-wise)
		for (unsigned int i = 0; i < extedges.size(); i++) extedges[i]->resort();

		printf("sorting %d extedges..\n", extedges.size());
		ImatiSTLsort(extedges.begin(), extedges.end(), &ExtEdge::lessThan);

		printf("making edges..");
		ExtEdge *e, *ep = NULL;
		coh_clusterEdge *ce;
		for (unsigned int i = 0; i < extedges.size(); i++)
		{
			e = extedges[i];
			if (!e->isTwoConnected())
			{
				if (ep == NULL || ExtEdge::lessThan(ep, e)) { ce = new coh_clusterEdge(e); ep = e; edges.appendTail(ce); } else ce->addEdge(e);
			}
		}
		printf("%d edges made.\n", edges.numels());

		for (n = clusters.head(); n != NULL; n = n->next()) ((coh_faceCluster *)n->data)->assignTriangleInfos();

		printf("making relations..");
		FOREACHNODE(edges, n)
		{
			ce = (coh_clusterEdge *)n->data;
			ce->initRelations();
		}

		printf("done.\n");
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// TEMPORARY CODE TO PRODUCE NICE COLORED IMAGES OF THE CLUSTERED MESH
	//
	/////////////////////////////////////////////////////////////////////////////////////////////////

	void saveClustering(const char *filename)
	{
		FILE *fp = fopen(filename, "w");
		fprintf(fp, "#VRML V1.0 ascii\n");
		fprintf(fp, "Separator {\n");
		fprintf(fp, "ShapeHints{ vertexOrdering COUNTERCLOCKWISE }\n");

		void **ocds;
		Vertex *v;
		Node *n;
		int i;
		fprintf(fp, " Coordinate3 {\n  point [\n");
		FOREACHVVVERTEX((&(tin->V)), v, n) fprintf(fp, "%f %f %f,\n", TMESH_TO_FLOAT(v->x), TMESH_TO_FLOAT(v->y), TMESH_TO_FLOAT(v->z));
		fprintf(fp, "]\n}\n");

		ocds = new void *[tin->V.numels()];
		i = 0; FOREACHVVVERTEX((&(tin->V)), v, n) ocds[i++] = v->info;
		i = 0; FOREACHVVVERTEX((&(tin->V)), v, n) v->info = (void *)(i++);

		coh_faceCluster *f;
		for (n = clusters.head(); n != NULL; n = n->next())
		{
			f = (coh_faceCluster *)n->data;
			f->printVRMLCluster(fp);
		}

		i = 0; FOREACHVVVERTEX((&(tin->V)), v, n) v->info = ocds[i++];
		delete[] ocds;

		fprintf(fp, "}\n");
		fclose(fp);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// END OF TEMPORARY CODE TO PRODUCE NICE COLORED IMAGES OF THE CLUSTERED MESH
	//
	/////////////////////////////////////////////////////////////////////////////////////////////////


	void killInvalidClusters(coh_faceCluster *c0, bool save_to_dangling =true)
	{
		coh_faceCluster *c;
		coh_clusterEdge *e;
		Node *m, *n;
		Triangle *t;
		Edge *me;

		List origEdges(c0->clusterEdges);

		c0->unlink();
		c0->mergeEdges();
		clusters.removeNode(c0);
		FOREACHVTTRIANGLE((&(c0->triangles)), t, m)
		{
			FOREACHTRIANGLEEDGE(t, me) if (me->isOnBoundary()) ((ExtEdge *)me->info)->removeTriangle(t);
			MARK_VISIT(t);
		}
		delete c0;

		tin->removeUnlinkedElements(); // To clean dangling elements from mergeEdges()
		if (save_to_dangling)
		{
			TriMesh *ntin = (TriMesh *)tin->createSubMeshFromSelection();
			tin->removeSelectedTriangles();
			if (ntin != NULL) dangling_tins.appendTail(ntin);
		}
		else tin->removeSelectedTriangles();
		n = edges.head();

		List neighborsToKill;
		for (n = origEdges.head(); n != NULL; n = n->next())
		{
			e = (coh_clusterEdge *)n->data;
			if (e->getNumClusters() == 1)
			{
				c = (coh_faceCluster *)e->getIncidentClusters()->head()->data;
				if (!c->isSelected()) { c->selected = true; neighborsToKill.appendTail(c); }
			}
		}

		while (n != NULL)
		{
			e = (coh_clusterEdge *)n->data;
			n = n->next();
			if (e->isIsolated())
			{
				edges.removeCell((n != NULL) ? (n->prev()) : (edges.tail()));
				delete e;
			}
		}

		while ((c = (coh_faceCluster *)neighborsToKill.popHead()) != NULL) killInvalidClusters(c, save_to_dangling);
	}

	TriMesh *getAllDanglingClusters()
	{
		TriMesh *cp, *ntin = (TriMesh *)tin->newObject();
		Node *n;
		FOREACHNODE(dangling_tins, n)
		{
			cp = (TriMesh *)n->data;
			if (cp->T.numels()) ntin->moveMeshElements(cp);
		}

		if (ntin->T.numels()) ntin->rebuildConnectivity();

		return ntin;
	}
};


int ric_floodSelectionAcrossIntersections(TriMesh& tin)
{
	List todo;
	Triangle *t, *s;
	Edge *e;
	int nv = 0;
	Node *n;

	FOREACHVTTRIANGLE((&tin.T), t, n) if (IS_VISITED(t)) todo.appendTail(t);

	while ((t = (Triangle *)todo.popHead()) != NULL)
	{
	 FOREACHTRIANGLEEDGE(t, e)
	 {
	  s = e->oppositeTriangle(t);
      if (s != NULL) { if (!IS_VISITED(s)) { MARK_VISIT(s); nv++; todo.appendTail(s); } }
	  else
	  {
		ExtEdge *xe = (ExtEdge *)e->info;
		FOREACHVTTRIANGLE((&xe->et), s, n) if (!IS_VISITED(s)) { MARK_VISIT(s); nv++; todo.appendTail(s); }
	  }
     }
	}

	return nv;
}


//Triangle *ric_getExtremeTriangle2(TriMesh& tin, unsigned char which_coord, bool max_or_min)
//{
//	Triangle *t, *gt = NULL;
//	Vertex *v;
//	coord max_c;
//	Node *n;
//
//	if (max_or_min)
//	{
//		max_c = -DBL_MAX;
//		FOREACHVVVERTEX((&(tin.V)), v, n) if (v->at(which_coord) > max_c) max_c = v->at(which_coord);
//	} else
//	{
//		max_c = DBL_MAX;
//		FOREACHVVVERTEX((&(tin.V)), v, n) if (v->at(which_coord) < max_c) max_c = v->at(which_coord);
//	}
//
//
//	FOREACHVVVERTEX((&(tin.V)), v, n) if (v->at(which_coord) == max_c) MARK_BIT(v, 5); else UNMARK_BIT(v, 5);
//
//	coord aznor, mxznor = -DBL_MAX;
//	//	double aznor, mxznor = -DBL_MAX;
//	coord rc;
//	Point nor;
//	FOREACHVTTRIANGLE((&tin.T), t, n)
//	if (IS_BIT(t->v1(), 5) || IS_BIT(t->v2(), 5) || IS_BIT(t->v3(), 5))
//	{
//		nor = t->getVector();
//		rc = nor[which_coord];
//		aznor = (rc*rc) / nor.squaredLength();
//		//aznor = t->getNormal()[which_coord].TMESH_TO_DOUBLE();
//		//aznor = FABS(aznor);
//		if (aznor>mxznor)  { mxznor = aznor; gt = t; }
//	}
//
//	if (gt == NULL) return NULL;
//
//	rc = gt->getVector()[which_coord];
//	if ((max_or_min && rc<0) || (!max_or_min && rc>0))
//		//	 if ((max_or_min && gt->getNormal()[which_coord].TMESH_TO_DOUBLE()<0) || (!max_or_min && gt->getNormal()[which_coord].TMESH_TO_DOUBLE()>0))
//	{
//		((coh_faceCluster *)gt->info)->flip();
//	}
//
//	return gt;
//}
//

int ric_xCompare(const void *a, const void *b)
{
	if ((((Point *)a)->x < ((Point *)b)->x)) return -1;
	else if ((((Point *)a)->x >((Point *)b)->x)) return 1;
	else return 0;
}

bool ric_xLessThan(Vertex *a, Vertex *b)
{
	return (((a)->x < (b)->x));
}

// Besides sorting the vertices, this function marks bit4 for all of them.
// Newly created vertices will be appended at the beginning of V and can
// be identified because they do not have bit4 set.

void sortVerticesByX(TriMesh& tin)
{
	std::vector<Vertex *> vtx;
	vtx.reserve(tin.V.numels());
	Vertex *v;
	Node *n;
	FOREACHVVVERTEX((&(tin.V)), v, n) {	UNMARK_BIT(v, 5); vtx.push_back(v);	}
	ImatiSTLsort(vtx.begin(), vtx.end(), ric_xLessThan);
	n = tin.V.head();
	for (unsigned int i = 0; i < vtx.size(); i++) { n->data = vtx[i]; n = n->next(); }
}


Triangle *ric_getExtremeTriangleOnX(TriMesh& tin)
{
	Vertex *v, *ov;
	coord max_c, aznor, mxznor;
	Node *n, *m;

	// 1) Find extreme vertices
	if (tin.V.numels() == 0) return NULL;
	v = (Vertex *)tin.V.head()->data;
	max_c = v->x; // max_c is the minimum X coordinate overall

	// 2) Find extreme edges
	Point nor;
	Triangle *gt = NULL;
	List *ve, cedges, candidedges, candidates;
	Edge *e, *ge = NULL;

	FOREACHVVVERTEX((&(tin.V)), v, n) if (v->x == max_c) // for each vertex having minimum X coord
	{ // add its incident edges to 'candidedges'. 'vertical' edges are considered for one of their endpoints only
		ve = v->VE();
		FOREACHVEEDGE(ve, e, m)
		{
			ov = e->oppositeVertex(v);
			if (ov->x != max_c || ov == e->v2) candidedges.appendTail(e);
		}
		delete ve;
	}
	else break;

	mxznor = DBL_MAX; // Now select the edge in 'candidedges' having minimum X component
	FOREACHVEEDGE((&candidedges), e, n)
	{
		cedges.appendTail(e);
		nor = e->toVector();
		aznor = (nor.x*nor.x) / nor.squaredLength();
		if (aznor<mxznor)  { mxznor = aznor; ge = e; }
	}

	FOREACHVEEDGE((&cedges), e, n) // For each such edge with minimum X component
	{ // select its incident triangles and put them in 'candidates'
		nor = e->toVector();
		aznor = (nor.x*nor.x) / nor.squaredLength();
		if (aznor==mxznor)
		{
			if (e->t1) candidates.appendTail(e->t1);
			if (e->t2) candidates.appendTail(e->t2);
		}
	}

	if (!candidates.numels()) return NULL;

	// 3) Find extreme triangle
	Triangle *t;
	mxznor = -DBL_MAX;
	FOREACHVTTRIANGLE((&candidates), t, n) // Now select the triangle in 'candidates' having minimum X normal component
	{
		nor = t->getVector();
		aznor = (nor.x*nor.x) / nor.squaredLength();
		if (aznor>mxznor)  { mxznor = aznor; gt = t; }
	}

	if (gt == NULL) return NULL;

	if (gt->getVector().x>0) ((coh_faceCluster *)gt->info)->flip();

	return gt;
}


TriMesh *coh_extractComponent2(TriMesh& tin, coh_clusterRep *brep)
{
	Triangle *t;
	List todo, done;
	bool is_coherent;
	coh_faceCluster *f, *fn;
	coh_clusterEdge *ce;
	Node *n, *m;
	bool iterate = false;

	do
	{
		if (!iterate) printf("\nTracking outer hull...");
		else { printf("."); fflush(stdout); }
		iterate = false;

		t = ric_getExtremeTriangleOnX(tin);
		if (t == NULL) return NULL;
		f = (coh_faceCluster *)t->info;
		todo.appendTail(f);
		f->selectTriangles();

		while ((f = (coh_faceCluster *)todo.popHead()) != NULL)
		{
			done.appendTail(f);
			m = f->represTriangles.head();
			FOREACHNODE(f->clusterEdges, n)
			{
				ce = (coh_clusterEdge *)n->data;
				fn = ce->getNextCluster(f, ((Triangle *)m->data), is_coherent);
				m = m->next();
				if (fn == NULL) // Boundary cluster. KILL f, deselect everything and restart
				{
					fn = f;
					done.joinTailList(&todo); while ((f = (coh_faceCluster *)done.popHead()) != NULL) f->deselectTriangles();
					brep->killInvalidClusters(fn);
					iterate = true; break;
				}
				if (!is_coherent) // Double orientation. KILL fn, deselect everything and restart
				{
					if (fn->isSelected())
					{
						done.joinTailList(&todo); while ((f = (coh_faceCluster *)done.popHead()) != NULL) f->deselectTriangles();
						brep->killInvalidClusters(fn);
						iterate = true; break;
					}
					else fn->flip();
				}
				if (!fn->isSelected()) { todo.appendTail(fn); fn->selectTriangles(); }
			}
		}
	} while (iterate);

	TriMesh *comp = (TriMesh *)tin.createSubMeshFromSelection();

	ric_floodSelectionAcrossIntersections(tin);
	tin.removeSelectedTriangles();

	return comp;
}





Triangle *ric_getExtremeTriangleOnX3(TriMesh& tin)
{
	Vertex *v, *ov;
	coord max_c, aznor, mxznor;
	Node *n, *m;

	// 1) Find extreme vertices
	if (tin.V.numels() == 0) return NULL;
	v = (Vertex *)tin.V.head()->data;
	max_c = v->x; // max_c is the minimum X coordinate overall

	// 2) Find extreme edges
	Point nor;
	Triangle *gt = NULL;
	List *ve, cedges, candidedges, candidates;
	Edge *e, *ge = NULL;

	FOREACHVVVERTEX((&(tin.V)), v, n) if (v->x == max_c) // for each vertex having minimum X coord
	{ // add its incident edges to 'candidedges'. 'vertical' edges are considered for one of their endpoints only
		ve = v->VE();
		FOREACHVEEDGE(ve, e, m)
		{
			ov = e->oppositeVertex(v);
			if (ov->x != max_c || ov == e->v2) candidedges.appendTail(e);
		}
		delete ve;
	} else break;

	mxznor = DBL_MAX; // Now select the edge in 'candidedges' having minimum X component
	FOREACHVEEDGE((&candidedges), e, n)
	{
		cedges.appendTail(e);
		nor = e->toVector();
		aznor = (nor.x*nor.x) / nor.squaredLength();
		if (aznor<mxznor)  { mxznor = aznor; ge = e; }
	}

	FOREACHVEEDGE((&cedges), e, n) // For each such edge with minimum X component
	{ // select its incident triangles and put them in 'candidates'
		nor = e->toVector();
		aznor = (nor.x*nor.x) / nor.squaredLength();
		if (aznor == mxznor)
		{
			if (e->t1) candidates.appendTail(e->t1);
			if (e->t2) candidates.appendTail(e->t2);
		}
	}

	if (!candidates.numels()) return NULL;

	// 3) Find extreme triangle
	Triangle *t;
	mxznor = -DBL_MAX;
	FOREACHVTTRIANGLE((&candidates), t, n) // Now select the triangle in 'candidates' having minimum X normal component
	{
		nor = t->getVector();
		aznor = (nor.x*nor.x) / nor.squaredLength();
		if (aznor>mxznor)  { mxznor = aznor; gt = t; }
	}

	if (gt == NULL) return NULL;

	if (gt->getVector().x > 0) tin.flipNormals(gt);

	return gt;
}


#define TRI_IS_UNTRAVERSED			0
#define TRI_HAS_OUTER_TRAVERSAL		2
#define TRI_HAS_INNER_TRAVERSAL		4

void coh_invertInnerOuterTraversal(Triangle *t0)
{
	Triangle *t, *y;
	Edge *e;
	List todo(t0);
	t = t0;
	t->mask = TRI_HAS_OUTER_TRAVERSAL;
	t->invert();

	while ((t = (Triangle *)todo.popHead()) != NULL)
	{
		FOREACHTRIANGLEEDGE(t, e)
		{
			y = e->oppositeTriangle(t);
			if (y == NULL || (t < y)) p_swap((void **)(&(e->v1)), (void **)(&(e->v2)));
			if (y != NULL && (y->mask != TRI_HAS_OUTER_TRAVERSAL))
			{
				y->mask = TRI_HAS_OUTER_TRAVERSAL; todo.appendTail(y); 
				y->invert();
			}
		}
	}
}

class traversedTriangle{
public:
	Triangle *t;
	bool outer;

	traversedTriangle(Triangle *y, bool o) { t = y; outer = o; t->mask |= ((o) ? (TRI_HAS_OUTER_TRAVERSAL) : (TRI_HAS_INNER_TRAVERSAL)); }

	bool isOuter() const { return outer; }
};

// Returns the leftmost M_solid component 'comp' and removes it from 'tin'.
// It also removes triangles which are edge-connected with 'comp' without being part of it.
// Possible sheets are added to M_sheet.
TriMesh *coh_extractComponent3(TriMesh& tin, TriMesh *m_sheet)
{
	Triangle *t, *y;
	Edge *e;
	ExtEdge *xe;
	List todo;
	bool is_coherent;
	Node *n;

	FOREACHVTTRIANGLE((&tin.T), t, n) t->mask = TRI_IS_UNTRAVERSED;

	t = ric_getExtremeTriangleOnX3(tin);
	traversedTriangle *tt;
	todo.appendHead((tt = new traversedTriangle(t, true)));

	while ((tt = (traversedTriangle *)todo.popHead()) != NULL)
	{
		t = tt->t;

		if (tt->isOuter())
		{
			FOREACHTRIANGLEEDGE(t, e)
			{
				if (e->isOnBoundary())
				{
					xe = (ExtEdge *)e->info;
					y = xe->getNextTriangle(t, is_coherent);
				} else { y = e->oppositeTriangle(t); is_coherent = true; }

				if (y == NULL) { ImatiSTL::error("coh_extractComponent3: should not happen"); }

				if (is_coherent)
				{
					if (!(y->mask & TRI_HAS_OUTER_TRAVERSAL)) todo.appendTail(new traversedTriangle(y, true));
				} else
				{
					if (!(y->mask & TRI_HAS_INNER_TRAVERSAL)) todo.appendTail(new traversedTriangle(y, false));
				}
			}
		}
		else
		{
			FOREACHTRIANGLEEDGE(t, e)
			{
				if (e->isOnBoundary())
				{
					xe = (ExtEdge *)e->info;
					y = xe->getPrevTriangle(t, is_coherent);
				} else { y = e->oppositeTriangle(t); is_coherent = true; }

				if (y == NULL) { ImatiSTL::error("coh_extractComponent3: should not happen"); }

				if (is_coherent)
				{
					if (!(y->mask & TRI_HAS_INNER_TRAVERSAL)) todo.appendTail(new traversedTriangle(y, false));
				} else
				{
					if (!(y->mask & TRI_HAS_OUTER_TRAVERSAL)) todo.appendTail(new traversedTriangle(y, true));
				}
			}
		}
		delete tt;
	}


// So, select TRI_HAS_DOUBLE_TRAVERSAL triangles and create a component out of them within m_sheet
	FOREACHVTTRIANGLE((&tin.T), t, n) if ((t->mask & TRI_HAS_OUTER_TRAVERSAL) && (t->mask & TRI_HAS_INNER_TRAVERSAL)) MARK_VISIT(t);



	//List lnodes;
	//Vertex *v;
	//FOREACHVTTRIANGLE((&tin.T), t, n) if (IS_VISITED(t))
	//{
	//	MARK_VISIT(t->e1); MARK_VISIT(t->e2); MARK_VISIT(t->e3);
	//	lnodes.appendTail(n);
	//}
	//if (lnodes.numels())
	//{
	//	while ((n = (Node *)lnodes.popHead()) != NULL) tin.T.moveNodeTo(n, &(m_sheet->T));
	//	FOREACHVEEDGE((&tin.E), e, n) if (IS_VISITED(e))
	//	{
	//		MARK_VISIT(e->v1); MARK_VISIT(e->v2);
	//		lnodes.appendTail(n);
	//	}
	//	while ((n = (Node *)lnodes.popHead()) != NULL) tin.E.moveNodeTo(n, &(m_sheet->E));
	//	FOREACHVVVERTEX((&tin.V), v, n) if (IS_VISITED(v)) lnodes.appendTail(n);
	//	while ((n = (Node *)lnodes.popHead()) != NULL) tin.V.moveNodeTo(n, &(m_sheet->V));
	//}




	TriMesh *comp = (TriMesh *)tin.createSubMeshFromSelection();
	if (comp != NULL)
	{
		FOREACHVTTRIANGLE((&tin.T), t, n) if (IS_VISITED(t))
		{
			FOREACHTRIANGLEEDGE(t, e)
			{
				ExtEdge *xe = (ExtEdge *)e->info;
				if (xe != NULL) xe->et.removeNode(t);
			}
		}
		tin.removeSelectedTriangles();
		m_sheet->moveMeshElements(comp);
	}

	FOREACHVTTRIANGLE((&tin.T), t, n)
	{
		if (t->mask & TRI_HAS_INNER_TRAVERSAL)
		{
			// Invert component and revert INNER to OUTER
			coh_invertInnerOuterTraversal(t);
		}
	}

	// Now, select TRI_HAS_OUTER_TRAVERSAL triangles and create a component out of them to be returned
	FOREACHVTTRIANGLE((&tin.T), t, n) if ((t->mask & TRI_HAS_OUTER_TRAVERSAL)) MARK_VISIT(t);
	comp = (TriMesh *)tin.createSubMeshFromSelection();

	// Then, flood the selection and remove from tin both the outer hull component and all its connected triangles
	ric_floodSelectionAcrossIntersections(tin);
	tin.removeSelectedTriangles();

	return comp;
}








//int ric_ptrTriangleCompare(const void *a, const void *b)
//{
//	Vertex *va1 = (Vertex *)((Triangle *)a)->v1()->info;
//	Vertex *va2 = (Vertex *)((Triangle *)a)->v2()->info;
//	Vertex *va3 = (Vertex *)((Triangle *)a)->v3()->info;
//	Vertex *vb1 = (Vertex *)((Triangle *)b)->v1()->info;
//	Vertex *vb2 = (Vertex *)((Triangle *)b)->v2()->info;
//	Vertex *vb3 = (Vertex *)((Triangle *)b)->v3()->info;
//	Vertex *v11 = MIN(va1, MIN(va2, va3));
//	Vertex *v13 = MAX(va1, MAX(va2, va3));
//	Vertex *v12 = (v11 != va1 && v13 != va1) ? (va1) : ((v11 != va2 && v13 != va2) ? (va2) : (va3));
//	Vertex *v21 = MIN(vb1, MIN(vb2, vb3));
//	Vertex *v23 = MAX(vb1, MAX(vb2, vb3));
//	Vertex *v22 = (v21 != vb1 && v23 != vb1) ? (vb1) : ((v21 != vb2 && v23 != vb2) ? (vb2) : (vb3));
//
//	if (v11 < v21) return 1;
//	else if (v11>v21) return -1;
//	if (v12 < v22) return 1;
//	else if (v12>v22) return -1;
//	if (v13 < v23) return 1;
//	else if (v13>v23) return -1;
//
//	return 0;
//}

bool ric_ptrTriangleLessThan(Triangle *a, Triangle *b)
{
	j_voidint va1 = (j_voidint )(a)->v1()->info;
	j_voidint va2 = (j_voidint )(a)->v2()->info;
	j_voidint va3 = (j_voidint )(a)->v3()->info;
	j_voidint vb1 = (j_voidint )(b)->v1()->info;
	j_voidint vb2 = (j_voidint )(b)->v2()->info;
	j_voidint vb3 = (j_voidint )(b)->v3()->info;
	j_voidint v11 = MIN(va1, MIN(va2, va3));
	j_voidint v13 = MAX(va1, MAX(va2, va3));
	j_voidint v12 = (v11 != va1 && v13 != va1) ? (va1) : ((v11 != va2 && v13 != va2) ? (va2) : (va3));
	j_voidint v21 = MIN(vb1, MIN(vb2, vb3));
	j_voidint v23 = MAX(vb1, MAX(vb2, vb3));
	j_voidint v22 = (v21 != vb1 && v23 != vb1) ? (vb1) : ((v21 != vb2 && v23 != vb2) ? (vb2) : (vb3));

	if (v11 < v21) return true;
	else if (v11>v21) return false;
	if (v12 < v22) return true;
	else if (v12>v22) return false;
	if (v13 < v23) return true;
	else return false;
}


bool vertexPtrLessThan(Vertex *a, Vertex *b) { return ((*a) < (*b)); }

int oh_removeDoubleTriangles(TriMesh& tin)
{
 // 1) Unify vertices
 Node *n;
 Vertex *v, *pv;
 j_voidint vidx = 1;

 std::vector<Vertex *> vtx; vtx.reserve(tin.V.numels());
 FOREACHVVVERTEX((&tin.V), v, n) vtx.push_back(v);
 ImatiSTLsort(vtx.begin(), vtx.end(), &vertexPtrLessThan);
 pv = vtx[0];
 for (std::vector<Vertex *>::iterator i = vtx.begin(); i != vtx.end(); i++)
 {
	 v = (*i);
	 //if ((*v) != (*pv)) pv = v;
	 //v->info = pv;
	 if ((*v) != (*pv)) { pv = v; vidx++; }
	 v->info = (void *)vidx;
 }

 // At this point any vertex points (through 'info') to its unique representative (possibly itself)

 tin.deselectTriangles(); // Prova a toglierla. forse non serve...

 Triangle *t;
 std::vector<Triangle *> tris; tris.reserve(tin.T.numels());
 FOREACHVTTRIANGLE((&(tin.T)), t, n) tris.push_back(t);
 ImatiSTLsort(tris.begin(), tris.end(), &ric_ptrTriangleLessThan);
 n = tin.T.head();
 for (std::vector<Triangle *>::iterator i = tris.begin(); i != tris.end(); i++) { n->data = (*i); n = n->next(); }

 Triangle *lt = NULL;
 int nit = 0;
 FOREACHVTTRIANGLE((&(tin.T)), t, n)
 {
	 if (lt != NULL && !ric_ptrTriangleLessThan(lt, t)) { MARK_VISIT(t); nit++; }
	 lt = t;
 }

 //for (n = tin.T.head(); n != NULL; n = oh_markAllCopies(n)) nit++;
 //nit = tin.T.numels() - nit;

 tin.removeSelectedTriangles();
 FOREACHVVVERTEX((&tin.V), v, n) v->info = NULL;

 return nit;
}


void coh_removeUnlinkedExtEdges(List& extedges)
{
	Node *n = extedges.head();
	while (n != NULL)
	{
		ExtEdge *xe = (ExtEdge *)n->data;
		n = n->next();
		if (xe->isIsolated()) extedges.freeCell((n == NULL) ? (extedges.tail()) : (n->prev()));
	}
}

inline void preLexEdgeCompare(Edge *e) { if (xyzCompare(e->v1, e->v2) > 0) e->invert(); }

bool vbesEdgeCompare(Edge *a, Edge *b)
{
	return (((*(a->v1)) < (*(b->v1))) || ((*(a->v1)) == (*(b->v1)) && ((*(a->v2)) < (*(b->v2)))));
}

inline bool vbesDifferentEdges(Edge *a, Edge *b)
{
	return (((*(a->v1)) != (*(b->v1))) || ((*(a->v2)) != (*(b->v2))));
}

inline bool ric_bboxContainmentCheck(Point& ma, Point& Ma, Point&mb, Point& Mb)
{
	return (ma.x>mb.x && Ma.x<Mb.x && ma.y>mb.y && Ma.y<Mb.y && ma.z>mb.z && Ma.z < Mb.z);
}

int ric_removeInnerShells(TriMesh *solids, TriMesh *sheets)
{
	int removed = 0;

	// 1) Split 'solids' into an array of shells (use TriMesh::split())
	std::vector<TriMesh *> solidShells;
	while (solids->T.numels()) solidShells.push_back((TriMesh *)solids->split());

	// 2) Same for sheets
	std::vector<TriMesh *> sheetShells;
	while (sheets->T.numels()) sheetShells.push_back((TriMesh *)sheets->split());

	// 3) For each solid 'i' check inclusion within all the other solids. Mark if inclusion occurs.
	std::vector<bool> markedSolids, markedSheets;
	markedSolids.assign(solidShells.size(), false);
	markedSheets.assign(sheetShells.size(), false);
	std::vector<Point> minbboxSolids, maxbboxSolids, minbboxSheets, maxbboxSheets;
	Point dummy;
	minbboxSolids.assign(solidShells.size(), dummy);
	maxbboxSolids.assign(solidShells.size(), dummy);
	minbboxSheets.assign(sheetShells.size(), dummy);
	maxbboxSheets.assign(sheetShells.size(), dummy);

	int i, j;
	for (i = 0; i < (int)solidShells.size(); i++) solidShells[i]->getBoundingBox(minbboxSolids[i], maxbboxSolids[i]);
	for (i = 0; i < (int)sheetShells.size(); i++) sheetShells[i]->getBoundingBox(minbboxSheets[i], maxbboxSheets[i]);

#pragma omp parallel for private(j) schedule(dynamic)
	for (i = 0; i < (int)solidShells.size(); i++)
	 for (j = 0; j < (int)solidShells.size(); j++)
	  if (i != j && !markedSolids[j] && ric_bboxContainmentCheck(minbboxSolids[i], maxbboxSolids[i], minbboxSolids[j], maxbboxSolids[j]) && 
		 solidShells[j]->isInnerPoint(((Triangle *)solidShells[i]->T.head()->data)->getCenter()))
	  {
		 markedSolids[i] = true; removed++;
		 j = solidShells.size();
	  }

	// 4) For each sheet, check for inclusion within all the remaining solids and mark if inclusion occurs.
#pragma omp parallel for private(j) schedule(dynamic)
	for (i = 0; i < (int)sheetShells.size(); i++)
	 for (j = 0; j < (int)solidShells.size(); j++)
	 if (i != j && !markedSolids[j] && ric_bboxContainmentCheck(minbboxSheets[i], maxbboxSheets[i], minbboxSolids[j], maxbboxSolids[j]) && 
		 solidShells[j]->isInnerPoint(((Triangle *)sheetShells[i]->T.head()->data)->getCenter()))
	 {
		 markedSheets[i] = true; removed++;
		 j = solidShells.size();
	 }

	 // 5) Rebuild 'solids' by merging the unmarked solid shells (use TriMesh::moveMeshElements())
	for (i = 0; i < (int)solidShells.size(); i++) if (!markedSolids[i]) solids->moveMeshElements(solidShells[i]);

	// 6) Rebuild 'sheets' by merging the unmarked sheet shells (use TriMesh::moveMeshElements())
	for (i = 0; i < (int)sheetShells.size(); i++) if (!markedSheets[i]) sheets->moveMeshElements(sheetShells[i]);

	return removed;
}


TriMesh *ric_computeOuterHull(TriMesh& tin, bool waf = false, int timeout_secs = 0)
{
 Node *n;
 Edge *e, *pe;
 ExtEdge *xe;
 List bes;
 int i;

 if (timeout_secs) ImatiSTL::exitOnTimeout(timeout_secs); // Init timeout counter
 ImatiSTL::printElapsedTime(true); // Init elapsed time reporting

 // First step: resolve all the intersections
 // The new mesh will be a set of manifold patches meeting at borders in a nonmanifold manner
 if (!tin.cutIntersections(waf, timeout_secs)) return NULL;

 // Only one copy of duplicated triangles can stay in the mesh to avoid ambiguous outer hull trackings
 ImatiSTL::info("Removing double triangles...\n");
 int ntb = oh_removeDoubleTriangles(tin);
 if (ntb) ImatiSTL::info("%d double triangles removed.\n",ntb);
 else ImatiSTL::info("No double triangles detected.\n");

 // Now create Extended Edges at boundaries. These can host a complete ET relation for nonmanifold configs
 std::vector<Edge *> vbes;
 FOREACHVEEDGE((&(tin.E)), e, n) if (e->isOnBoundary()) { vbes.push_back(e); preLexEdgeCompare(e); } else e->info = NULL;
 ImatiSTLsort(vbes.begin(), vbes.end(), &vbesEdgeCompare);

 ImatiSTL::info("Creating triangle fans...\n");
 ci_ExtEdges vextedges;
 pe = NULL;
 for (i = 0; i < (int)vbes.size(); i++)
 {
	 e = vbes[i];
	 if (pe == NULL || vbesDifferentEdges(e, pe)) { xe = new ExtEdge(e); vextedges.push_back(xe); pe = e; }
	 xe->addUnsortedTriangle(e->getBoundaryTriangle());
	 e->info = xe;
 }

 bool old_status = ImatiSTL::isUsingRationals();
 ImatiSTL::useRationals(true);

#pragma omp parallel for schedule(dynamic)
 for (i = 0; i < (int)vextedges.size(); i++) vextedges[i]->sortTriangles();

 ImatiSTL::printElapsedTime();

 // Now extract all the connected components constituting the outer hull.
 // 'tin' is made empty during this operation.
 TriMesh *comp, *ttin = new TriMesh;

 // This subdivides the mesh in 2-edge-connected clusters, and removes triangles belonging to boundary clusters.

//#define USE_CLUSTER_REP

#ifdef USE_CLUSTER_REP
 coh_clusterRep *brep = new coh_clusterRep(&tin, vextedges);
 if (waf) brep->saveClustering("clusteredPolyhedron.wrl");
 sortVerticesByX(tin);

 while (tin.T.numels())
 {
  comp = coh_extractComponent2(tin, brep);
  if (comp != NULL) { printf("%d tris.\n", comp->T.numels()); ttin->moveMeshElements(comp); }
 }
 printf("\n");
 comp = brep->getAllDanglingClusters();

#else
 comp = new TriMesh;
 TriMesh *m_solid_comp;
 sortVerticesByX(tin);
 while (tin.T.numels())
 {
  m_solid_comp = coh_extractComponent3(tin, comp);
  if (m_solid_comp != NULL) { printf("%d tris.\n", m_solid_comp->T.numels()); ttin->moveMeshElements(m_solid_comp); }
  else printf("0 tris.\n");
 }
 if (comp != NULL) comp->rebuildConnectivity();
#endif

 // Move the outer hulls back in 'tin' and rebuild a proper connectivity
 if (ttin!=NULL && ttin->T.numels()) tin.moveMeshElements(ttin);

 tin.rebuildConnectivity();

 ImatiSTL::printElapsedTime();
 int r = ric_removeInnerShells(&tin, comp);

 if (r) ImatiSTL::info("%d inner shells have been removed.\n", r);
 ImatiSTL::printElapsedTime();

 ImatiSTL::useRationals(old_status);

 r = tin.removeRedundantVertices(true);
 if (r) ImatiSTL::info("%d flat vertices removed from solid.\n", r);

 if (comp != NULL)
 {
	 comp->rebuildConnectivity();
	 r = comp->removeRedundantVertices(true);
	 if (r) ImatiSTL::info("%d flat vertices removed from sheets.\n", r);
 }

 return comp;
}

TriMesh *TriMesh::computeOuterHull(bool print_debug_info, int timeout_secs)
{
	return ric_computeOuterHull(*this, print_debug_info, timeout_secs);
}
