//#include "stdafx.h"
#include <fstream>
#include <iostream>

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"
//#include <omp.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/centroid.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/lloyd_optimize_mesh_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Segment_2										Segment_2;
typedef K::Ray_2											Ray_2;
typedef K::Line_2											Line_2;
typedef K::Point_2											Point_2;
typedef CGAL::Delaunay_mesh_vertex_base_2<K>                Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K>                  Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        Tds;

typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K> Vb1;
typedef CGAL::Triangulation_data_structure_2<Vb1>			Tds1;
typedef CGAL::Delaunay_triangulation_2<K, Tds1>				DT;

typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT, AT, AP>                                  VD;
typedef AT::Site_2                    Site_2;
typedef VD::Locate_result             Locate_result;
typedef VD::Ccb_halfedge_circulator   Ccb_halfedge_circulator;

typedef std::vector<std::vector<double>> segments_lists;

class MexFunction : public matlab::mex::Function {
public:
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {

		matlab::data::TypedArray<double> pnts = std::move(inputs[0]);
		int idx		= inputs[1][0];
		double nelx = inputs[2][0];
		double nely = inputs[3][0];

		idx -= 1; // matlab start with 1

		unsigned int pnts_num = pnts.getDimensions()[0];

		// matlab::data::ArrayFactory factory;
		//matlab::data::TypedArray<double> seeds = factory.createArray<double>({ pnts_num, 2 });
		//matlab::data::TypedArray<int> edges;
		//matlab::data::TypedArray<double> nodes;
		// matlab::data::CellArray edges_cell = factory.createCellArray({ pnts_num, 1 });
		// matlab::data::CellArray nodes_cell = factory.createCellArray({ pnts_num, 1 });

		matlab::data::CellArray results = voronoi_local(pnts, idx, nelx, nely);

		outputs[0] = results[0];
		outputs[1] = results[1];
		//outputs[2] = seeds;
	}

	void DT2VD_edgs(DT dt, segments_lists & rayEdge, segments_lists & halfEdge)
	{
		// pre-process to get rayEdgeBound for voronoi diagram
		//segments_lists rayEdge, halfEdge;
		for (DT::Edge_iterator eit = dt.edges_begin(); eit != dt.edges_end(); eit++) {
			CGAL::Object o = dt.dual(eit);
			std::vector<double> epi;
			if (CGAL::object_cast<K::Segment_2>(&o)) //如果这条边是线段，则绘制线段
			{
				epi.push_back(CGAL::object_cast<K::Segment_2>(&o)->source().x());
				epi.push_back(CGAL::object_cast<K::Segment_2>(&o)->source().y());
				epi.push_back(CGAL::object_cast<K::Segment_2>(&o)->target().x());
				epi.push_back(CGAL::object_cast<K::Segment_2>(&o)->target().y());
				halfEdge.push_back(epi);

				//Eigen::RowVector2d P1(epi[0], epi[1]);
				//Eigen::RowVector2d P2(epi[2], epi[3]);
				//viewer.data().add_edges(P1, P2, Eigen::RowVector3d(1, 1, 1));
			}
			else if (CGAL::object_cast<K::Ray_2>(&o))//如果这条边是射线，则绘制射线
			{
				epi.push_back(CGAL::object_cast<K::Ray_2>(&o)->source().x());
				epi.push_back(CGAL::object_cast<K::Ray_2>(&o)->source().y());
				epi.push_back(CGAL::object_cast<K::Ray_2>(&o)->point(1).x());
				epi.push_back(CGAL::object_cast<K::Ray_2>(&o)->point(1).y());
				rayEdge.push_back(epi);

				/*		Eigen::RowVector2d P1(epi[0], epi[1]);
				Eigen::RowVector2d P2(epi[2], epi[3]);
				viewer.data().add_edges(P1, P2, Eigen::RowVector3d(1, 1, 1));*/
			}
		}
	}

	bool point_equal(Point_2 p1, Point_2 p2)
	{
		if (abs(p1.x() - p2.x()) < 1e-5 && abs(p1.y() - p2.y()) < 1e-5)
			return 1;
		else
			return 0;
	}

	//ps is the source of the ray, pt is the target of the ray
	Point_2 VoronoiD_CrossPoint(Ray_2 ray, double nelx, double nely, std::vector<Segment_2> boundarys) {

		for (int i = 0; i < 4; ++i)
		{
			const auto result = intersection(ray, boundarys[i]);
			if (result) {
				if (const Point_2* p = boost::get<Point_2 >(&*result))
					return *p;
			}
		}
		std::cout << "Error!no cross point!";
		Point_2 result(0, 0);
		return result;
	}

	// redivide the segments & rays into boundary
	std::vector<Segment_2> VoronoiD_UpdateEdge(double nelx, double nely,
		segments_lists rayEdge, segments_lists halfEdge, std::vector<Segment_2> boundarys)
	{
		segments_lists halfEdge_r, rayEdge_r;

		for (int i = 0; i < halfEdge.size(); i++) {
			double x1 = halfEdge[i][0];
			double y1 = halfEdge[i][1];
			double x2 = halfEdge[i][2];
			double y2 = halfEdge[i][3];
			std::vector<double> halfEdage_i;
			if (abs(x1) <= nelx && abs(y1) <= nely && abs(x2) <= nelx && abs(y2) <= nely) {
				halfEdage_i.insert(halfEdage_i.end(), halfEdge[i].begin(), halfEdge[i].end());
				halfEdge_r.push_back(halfEdage_i);
			}
			else if (abs(x1) <= nelx && abs(y1) <= nely) {
				halfEdage_i.push_back(x1);
				halfEdage_i.push_back(y1);
				halfEdage_i.push_back(x2);
				halfEdage_i.push_back(y2);
				rayEdge_r.push_back(halfEdage_i);
			}
			else if (abs(x2) <= nelx && abs(y2) <= nely) {
				halfEdage_i.push_back(x2);
				halfEdage_i.push_back(y2);
				halfEdage_i.push_back(x1);
				halfEdage_i.push_back(y1);
				rayEdge_r.push_back(halfEdage_i);
			}
			else {
				continue;
			}
		}

		for (int i = 0; i < rayEdge.size(); i++) {
			if (abs(rayEdge[i][0]) < nelx && abs(rayEdge[i][1]) < nely) {
				rayEdge_r.push_back(rayEdge[i]);
			}
			else {
				continue;
			}
		}

		std::vector<Segment_2> rayEdgeBound;
		for (int i = 0; i < rayEdge_r.size(); i++) {
			Point_2 ps(rayEdge_r[i][0], rayEdge_r[i][1]);
			Point_2 pt(rayEdge_r[i][2], rayEdge_r[i][3]);
			Ray_2 ray(ps, pt);
			Point_2 ptn = VoronoiD_CrossPoint(ray, nelx, nely, boundarys); // intersection of the ray with the boundary

			Segment_2 rayEdge_i(ps, ptn);
			rayEdgeBound.push_back(rayEdge_i);
		}
		return rayEdgeBound;
	}

	void add_boundary_edge(std::vector<Segment_2> added_rays, std::vector<Segment_2> & boundary_edges,
		double halfx, double halfy)
	{
		if (added_rays.size() != 2)
			throw "The size of added_rays is not 2!";

		// add the edge on the bounding box
		Segment_2 seg1 = added_rays[0];
		Segment_2 seg2 = added_rays[1];

		K::Vector_2 dir1 = seg1.point(1) - seg1.point(0);
		K::Vector_2 dir2 = seg2.point(1) - seg2.point(0);

		// two situations:
		Point_2 pf1, pf2;
		bool flag = 1;
		for (int i = 0; i < 2; ++i)
		{
			Point_2 p1 = seg1.point(i);
			for (int j = 0; j < 2; ++j)
			{
				Point_2 p2 = seg2.point(j);

				if ((abs(abs(p1.x()) - halfx) < 1e-5 && abs(abs(p2.x()) - halfx) < 1e-5) ||
					(abs(abs(p1.y()) - halfy) < 1e-5 && abs(abs(p2.y()) - halfy) < 1e-5))
				{
					Segment_2 s(p2, p1);
					boundary_edges.push_back(s);
					flag = 0;
					return;
				}
				else
				{
					if (abs(abs(p1.x()) - halfx) < 1e-5 || abs(abs(p1.y()) - halfy) < 1e-5)
						pf1 = p1;

					if (abs(abs(p2.x()) - halfx) < 1e-5 || abs(abs(p2.y()) - halfy) < 1e-5)
						pf2 = p2;
				}
			}
		}
		
		if (flag)
		{
			Point_2 cand1(pf1.x(), pf2.y());
			Segment_2 candidate11(cand1, pf1);
			Segment_2 candidate12(cand1, pf2);
			double dist1 = 0;

			Point_2 cand2(pf2.x(), pf1.y());
			Segment_2 candidate21(cand2, pf1);
			Segment_2 candidate22(cand2, pf2);
			double dist2 = 0;

			for (int i = 0; i < boundary_edges.size(); ++i)
			{
				dist1 += CGAL::squared_distance(boundary_edges[i], cand1);
				dist2 += CGAL::squared_distance(boundary_edges[i], cand2);
			}

			if (dist1 > dist2)
			{
				boundary_edges.push_back(candidate11);
				boundary_edges.push_back(candidate12);
			}
			else
			{
				boundary_edges.push_back(candidate21);
				boundary_edges.push_back(candidate22);
			}
		}
	}

	void add_ray(Point_2 p, Ccb_halfedge_circulator ec, std::vector<Segment_2> rayEdgeBound,
		std::vector<Segment_2> & boundary_edges, std::vector<Segment_2> & added_rays)
	{
		VD::Delaunay_vertex_handle v1 = ec->up();
		VD::Delaunay_vertex_handle v2 = ec->down();

		K::Vector_2 direction(v1->point().y() - v2->point().y(),
			v2->point().x() - v1->point().x());

		for (int i = 0; i < rayEdgeBound.size(); ++i)
		{
			K::Vector_2 r1(rayEdgeBound[i]);
			Point_2 p0 = rayEdgeBound[i].point(0);
			Point_2 p1 = rayEdgeBound[i].point(1);

			if ((abs(r1.y() * direction.x() - r1.x() * direction.y()) < 1e-5) &&
				(point_equal(p, p1) || point_equal(p, p0)))
			{
				boundary_edges.push_back(rayEdgeBound[i]);
				added_rays.push_back(rayEdgeBound[i]);
				break;
			}
		}
	}

	void add_segment(Ccb_halfedge_circulator ec, std::vector<Segment_2> & boundary_edges, bool & flag,
		std::vector<Segment_2> & added_rays, double halfx, double halfy, std::vector<Segment_2> boundarys)
	{
		double xs = ec->source()->point()[0];
		double ys = ec->source()->point()[1];
		double xt = ec->target()->point()[0];
		double yt = ec->target()->point()[1];
		Segment_2 s(ec->source()->point(), ec->target()->point());

		// the segment maybe out of the box
		if (abs(xs) <= halfx + 1e-5 && abs(ys) <= halfy + 1e-5 &&
			abs(xt) <= halfx + 1e-5 && abs(yt) <= halfy + 1e-5)
		{
			boundary_edges.push_back(s);
		}
		else
		{
			// find the intersection point
			Point_2 p1(0, 0);
			for (int i = 0; i < boundarys.size(); ++i)
			{
				const auto result = intersection(s, boundarys[i]);
				if (result) {
					if (const Point_2* p = boost::get<Point_2 >(&*result))
					{
						p1 = *p;
						break;
					}
				}
			}

			// find another point to build the cut segment
			if (abs(xs) <= halfx + 1e-5 && abs(ys) <= halfy + 1e-5)
			{
				Segment_2 s1(ec->source()->point(), p1);
				boundary_edges.push_back(s1);
				added_rays.push_back(s1);
			}
			else if (abs(xt) <= halfx + 1e-5 && abs(yt) <= halfy + 1e-5)
			{
				Segment_2 s1(ec->target()->point(), p1);
				boundary_edges.push_back(s1);
				added_rays.push_back(s1);
			}
			flag = 1;
		}
	}

	void reorder_edges(std::vector<Point_2> pnts, std::vector<std::vector<int>> & edges_idx)
	{
		// reorder the edge_idx of end to end
		std::vector<std::vector<int>> new_edges = edges_idx;
		edges_idx.erase(edges_idx.begin());

		for (int i = 0; i < new_edges.size() - 1; ++i)
		{
			std::vector<int> ei = new_edges[i];

			for (int j = 0; j < edges_idx.size(); ++j)
			{
				std::vector<int> ej = edges_idx[j];
				if (ei[1] == ej[0])
				{
					new_edges[i + 1] = ej;
					edges_idx.erase(edges_idx.begin() + j);
					break;
				}
				else
					if ((ei[1] == ej[1]))
					{
						new_edges[i + 1][0] = ej[1];
						new_edges[i + 1][1] = ej[0];
						edges_idx.erase(edges_idx.begin() + j);
						break;
					}
			}
		}
		edges_idx.clear();
		edges_idx = new_edges;
	}

	void preprocess_vcell(std::vector<Segment_2> edges, std::vector<Point_2> & pnts,
		std::vector<std::vector<int>> & edges_idx)
	{
		int p1_idx = 0, p2_idx = 1;
		for (int i = 0; i < edges.size(); ++i)
		{
			Point_2 p1 = edges[i].point(0);
			Point_2 p2 = edges[i].point(1);

			bool flag1 = 1, flag2 = 1;
			for (int j = 0; j < pnts.size(); ++j)
			{
				Point_2 p0 = pnts[j];
				if (point_equal(p0, p1)) {
					flag1 = 0;
					p1_idx = j;
				}
				if (point_equal(p0, p2)) {
					flag2 = 0;
					p2_idx = j;
				}
			}
			if (flag1) {
				pnts.push_back(p1);
				p1_idx = pnts.size() - 1;
			}
			if (flag2) {
				pnts.push_back(p2);
				p2_idx = pnts.size() - 1;
			}
			edges_idx[i] = std::vector<int>{ p1_idx, p2_idx };
		}
		Point_2 center = CGAL::centroid(pnts.begin(), pnts.end(), CGAL::Dimension_tag<0>());

		reorder_edges(pnts, edges_idx);
	}

	std::vector<Segment_2> voronoi_single(Point_2 p, VD vd, double halfx, double halfy, 
		std::vector<Segment_2> rayEdgeBound, std::vector<Segment_2> boundarys)
	{
		std::vector<Segment_2> boundary_edges;

		std::vector<Segment_2> added_rays;

		Locate_result lr = vd.locate(p); // the face of p located
		VD::Face_handle* f = boost::get<VD::Face_handle>(&lr);

		Ccb_halfedge_circulator ec_start = (*f)->ccb(); // traversing the halfedges on the boundary of f
		Ccb_halfedge_circulator ec = ec_start;

		bool flag = 0;
		do {
			//add_segments_and_update_bounding_box(ec);
			if (ec->is_bisector())
			{
				// TODO:
				std::cout << "bisector!" <<std::endl;
				throw "Bisector!";
			}
			else
			{
				if (ec->is_ray())
				{
					if (ec->has_source())
					{
						Point_2 p1 = ec->source()->point();
						add_ray(p1, ec, rayEdgeBound, boundary_edges, added_rays);
						flag = 1;
					}
					if (ec->has_target())
					{
						Point_2 p1 = ec->target()->point();
						add_ray(p1, ec, rayEdgeBound, boundary_edges, added_rays);
						flag = 1;
					}
				}
				else if (ec->is_segment())
				{
					add_segment(ec, boundary_edges, flag, added_rays, halfx, halfy, boundarys);
				}
			}
		} while (++ec != ec_start);

		// add the edge on the bounding box
		if (flag) // including boundary rays
			add_boundary_edge(added_rays, boundary_edges, halfx, halfy);

		return boundary_edges;
	}

	matlab::data::CellArray voronoi_local(matlab::data::TypedArray<double> seeds, int idx, double nelx, double nely)
		//matlab::data::TypedArray<int> & edges_local, matlab::data::TypedArray<double> & nodes_local)
	{
		DT dt;
		VD vd;
		std::vector<std::pair<Point_2, unsigned>> points;

		unsigned int pnts_num = seeds.getDimensions()[0];

		for (int i = 0; i < pnts_num; ++i)
		{
			Point_2 p(seeds[i][0], seeds[i][1]);
			points.push_back(std::make_pair(p, i));

			Site_2 t(p.x(), p.y());
			vd.insert(t);
		}
		dt = DT(points.begin(), points.end());
		
		// get the ray-edge and segment-edge of the voronoi diagram from the DT
		segments_lists rayEdge, halfEdge;
		DT2VD_edgs(dt, rayEdge, halfEdge);

		// cut the ray-edge under the boundary of [-halfx,-halfy] & [halfx, halfy]
		double halfx = nelx / 2.;
		double halfy = nely / 2.;

		Segment_2 up(Point_2(halfx, halfy), Point_2(-halfx, halfy));
		Segment_2 down(Point_2(-halfx, -halfy), Point_2(halfx, -halfy));
		Segment_2 left(Point_2(-halfx, -halfy), Point_2(-halfx, halfy));
		Segment_2 right(Point_2(halfx, halfy), Point_2(halfx, -halfy));

		std::vector<Segment_2> boundarys;
		boundarys.push_back(up);
		boundarys.push_back(down);
		boundarys.push_back(left);
		boundarys.push_back(right);
		std::vector<Segment_2> rayEdgeBound = VoronoiD_UpdateEdge(halfx, halfy, rayEdge, halfEdge, boundarys);

		// step.3 get voronoi		
		Point_2 p = points[idx].first;

		// get the v_edges of each cell
		std::vector<Segment_2> voronoi_edges = voronoi_single(p, vd, halfx, halfy, rayEdgeBound, boundarys);

		// the new_pnts is ordered to assemble polygon
		std::vector<Point_2> pnts;
		std::vector<std::vector<int>> edges_idx(voronoi_edges.size());
		preprocess_vcell(voronoi_edges, pnts, edges_idx);

		matlab::data::ArrayFactory factory;
		matlab::data::TypedArray<int> edges_local = factory.createArray<int>({ edges_idx.size(), 2 });
		matlab::data::TypedArray<double> nodes_local = factory.createArray<double>({ pnts.size(), 2 });
		for (int i = 0; i < edges_idx.size(); ++i)
		{
			edges_local[i][0] = edges_idx[i][0] + 1;
			edges_local[i][1] = edges_idx[i][1] + 1;

			nodes_local[i][0] = pnts[i][0];
			nodes_local[i][1] = pnts[i][1];
		}

		 matlab::data::CellArray result_cell = factory.createCellArray({ 2, 1 });
		 result_cell[ 0 ] = edges_local;
		 result_cell[ 1 ] = nodes_local;
		 return result_cell;
	}
};
