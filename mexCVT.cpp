#define CGAL_USE_BASIC_VIEWER

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Gmpq.h>
//#include <CGAL/Lazy_exact_nt.h>
//#include <CGAL/Simple_cartesian.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>

#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Polygon_2.h>

#include <iostream>
#include <cstdint>
#include <fstream>
#include <iostream>

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"
#include <omp.h>

//typedef CGAL::Lazy_exact_nt<CGAL::Gmpq> FT;
//typedef CGAL::Simple_cartesian<FT> K;
typedef CGAL::Exact_predicates_exact_constructions_kernel	K;

typedef K::Iso_rectangle_2									Iso_rectangle2;
typedef K::Segment_2										Segment2;
typedef K::Ray_2											Ray2;
typedef K::Line_2											Line2;
typedef K::Point_2											Point2;
typedef CGAL::Polygon_2<K>									Polygon2;
typedef Polygon2::Edge_const_iterator						EdgeIterator;

typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K> Vb1;
typedef CGAL::Triangulation_data_structure_2<Vb1>					 Tds1;
typedef CGAL::Delaunay_triangulation_2<K, Tds1>						 DT;

typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT, AT, AP>                                  VD;
typedef VD::Locate_result					Locate_result;
typedef VD::Face::Ccb_halfedge_circulator   Ccb_halfedge_circulator;

//Used to convert otherwise infinite rays into looooong line segments
const int RAY_LENGTH = 100;

//CGAL often returns objects that are either segments or rays. This converts
//these objects into segments. If the object would have resolved into a ray,
//that ray is intersected with the bounding box defined above and returned as
//a segment.
const auto ConvertToSeg = [&](const CGAL::Object seg_obj, bool outgoing) -> K::Segment_2
{
	//One of these will succeed and one will have a NULL pointer
	const K::Segment_2 *dseg = CGAL::object_cast<K::Segment_2>(&seg_obj);
	const K::Ray_2     *dray = CGAL::object_cast<K::Ray_2>(&seg_obj);
	if (dseg) { //Okay, we have a segment
		return *dseg;
	}
	else {    //Must be a ray
		const auto &source = dray->source();
		const auto dsx = source.x();
		const auto dsy = source.y();
		const auto &dir = dray->direction();
		const auto tpoint = K::Point_2(dsx + RAY_LENGTH * dir.dx(), dsy + RAY_LENGTH * dir.dy());
		if (outgoing)
			return K::Segment_2(
				dray->source(),
				tpoint
			);
		else
			return K::Segment_2(
				tpoint,
				dray->source()
			);
	}
};

class MexFunction : public matlab::mex::Function {
public:
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {

		matlab::data::TypedArray<double> pnts = std::move(inputs[0]);
		matlab::data::TypedArray<int> idx = std::move(inputs[1]);
		double nelx = inputs[2][0];
		double nely = inputs[3][0];

		std::vector<Segment2> edges = mexCVT(pnts, idx, nelx, nely);

		matlab::data::ArrayFactory factory;
		matlab::data::TypedArray<double> edges_mat = factory.createArray<double>({ edges.size(), 4 });

//#pragma omp parallel for
		for (int i = 0; i < edges.size(); ++i)
		{
			edges_mat[i][0] = CGAL::to_double(edges[i].point(0).x());
			edges_mat[i][1] = CGAL::to_double(edges[i].point(0).y());
			edges_mat[i][2] = CGAL::to_double(edges[i].point(1).x());
			edges_mat[i][3] = CGAL::to_double(edges[i].point(1).y());
		}

		outputs[0] = edges_mat;
		//outputs[1] = nodes_cell;
		//outputs[2] = seeds;
	}

	std::vector<Segment2> mexCVT(matlab::data::TypedArray<double> seeds,
		matlab::data::TypedArray<int> idx, double nelx, double nely)
	{
		std::vector<std::pair<Point2, unsigned>> points;

		unsigned int pnts_num = seeds.getDimensions()[0];

//#pragma omp parallel for
		for (int i = 0; i < pnts_num; ++i)
		{
			double x = seeds[i][0];
			double y = seeds[i][1];
			Point2 p(x, y);
			//Point2 p(seeds[i][0], seeds[i][1]);
			points.push_back(std::make_pair(p, i));
		}
		DT dt = DT(points.begin(), points.end());
		dt.is_valid();

		VD vd(dt);

		/// step.2 generate the outer-bounding
		// cut the ray-edge under the boundary of [-halfx,-halfy] & [halfx, halfy]
		double halfx = nelx / 2.;
		double halfy = nely / 2.;

		//In order to crop the Voronoi diagram, we need to convert the bounding box into a polygon. 
		CGAL::Polygon_2<K> bpoly;
		bpoly.push_back(Point2(-halfx, -halfy));
		bpoly.push_back(Point2(halfx, -halfy));
		bpoly.push_back(Point2(halfx, halfy));
		bpoly.push_back(Point2(-halfx, halfy));

		/// step.3 get voronoi	
		std::vector<Segment2> edges;

		unsigned int nidx = idx.getDimensions()[0];
#pragma omp parallel for
		for (int j = 0; j < nidx; ++j)
		{
			int k = idx[j] - 1;
			//std::cout << k << std::end;
			Point2 p = points[k].first;

			Locate_result lr = vd.locate(p); // the face of p located
			VD::Face_handle* f = boost::get<VD::Face_handle>(&lr);

			Ccb_halfedge_circulator ec_start = (*f)->ccb(); // traversing the halfedges on the boundary of f
			Ccb_halfedge_circulator ec = ec_start;

			CGAL::Polygon_2<K> pgon;

			do {
				//A half edge circulator representing a ray doesn't carry direction
				//information. To get it, we take the dual of the dual of the half-edge.
				//The dual of a half-edge circulator is the edge of a Delaunay triangle.
				//The dual of the edge of Delaunay triangle is either a segment or a ray.
				const CGAL::Object seg_dual = vd.dual().dual(ec->dual());

				//Convert the segment/ray into a segment
				const auto this_seg = ConvertToSeg(seg_dual, ec->has_target());

				pgon.push_back(this_seg.source());

				//If the segment has no target, it's a ray. This means that the next
				//segment will also be a ray. We need to connect those two rays with a
				//segment. The following accomplishes this.
				if (!ec->has_target()) {
					const CGAL::Object nseg_dual = vd.dual().dual(ec->next()->dual());
					const auto next_seg = ConvertToSeg(nseg_dual, ec->next()->has_target());
					pgon.push_back(next_seg.target());
				}
			} while (++ec != ec_start); //Loop until we get back to the beginning

			//if (CGAL::do_intersect(pgon, bpoly)) // in fact, not necessarily
			{
				//Perform the intersection. Since CGAL is very general, it believes the
				//result might be multiple polygons with holes.
				std::list<CGAL::Polygon_with_holes_2<K>> isect;
				CGAL::intersection(pgon, bpoly, std::back_inserter(isect));

				//But we know better. The intersection of a convex polygon and a box is
				//always a single polygon without holes. Let's assert this.
				assert(isect.size() == 1);

				//And recover the polygon of interest
				auto &poly_w_holes = isect.front();
				auto &poly_outer = poly_w_holes.outer_boundary();

				//Print the polygon as a WKT polygon
				//std::cout << fnum << ", "
				//	"\"POLYGON ((";
				for (auto v = poly_outer.vertices_begin(); v+1 != poly_outer.vertices_end(); v++)
				{
					Point2 p1(v->x(), v->y());
					Point2 p2((v+1)->x(), (v+1)->y());
					edges.push_back(Segment2(p1,p2));
				}
					//std::cout << v->x() << " " << v->y() << ", ";
				//std::cout << poly_outer.vertices_begin()->x() << " " << poly_outer.vertices_begin()->y() << "))\"\n";
			}
		}
		return edges;
	}
};