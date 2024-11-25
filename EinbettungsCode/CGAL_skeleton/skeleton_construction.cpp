#include<vector>
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_with_holes_2.h>
#include<CGAL/create_offset_polygons_2.h>
#include <CGAL/create_offset_polygons_from_polygon_with_holes_2.h>
#include<CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>
#include<CGAL/create_weighted_straight_skeleton_from_polygon_with_holes_2.h>
#include<CGAL/Polygon_2.h>
#include <CGAL/squared_distance_2.h>
#include <float.h>
#include <cmath>

int buffer = 4;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;

typedef K::Point_2                    Point_2;
typedef CGAL::Polygon_2<K>            Polygon_2;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes;
typedef CGAL::Straight_skeleton_2<K>  Ss;

typedef std::shared_ptr<Ss> SsPtr;

typedef Ss::Vertex_const_handle		Vertex_const_handle;
typedef Ss::Vertex_const_iterator		Vertex_const_iterator;
typedef Ss::Halfedge_const_handle	Halfedge_const_handle;
typedef Ss::Halfedge_const_iterator	Halfedge_const_iterator;

typedef std::shared_ptr<Polygon_2> PolygonPtr;
typedef std::vector<PolygonPtr> PolygonPtrVector;

struct OffsetPolygon {
    double* coordinates;
    int length;
};

OffsetPolygon offset(double divisor,unsigned long outer_count, double* outer_points, unsigned long hole_count, double* hole_points) {
    Polygon_2 polygon;
    Polygon_2 hole;
    for(int i = buffer; i + 1 < outer_count; i+= 2) {
        polygon.push_back(Point_2(outer_points[i], outer_points[i+1]));
    }
    for(int i = buffer; i + 1 < hole_count; i+= 2) {
        hole.push_back(Point_2(hole_points[i], hole_points[i+1]));
    }
    assert(polygon.is_counterclockwise_oriented());
    Polygon_with_holes poly( polygon ) ;
    if (hole_count > 0) {
        assert(hole.is_clockwise_oriented());
        poly.add_hole( hole ) ;
    }
    SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly);
    Ss& skeleton = *iss;
    double first = 0.0;
    //CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*iss);
    for (Halfedge_const_iterator i = skeleton.halfedges_begin(); i != skeleton.halfedges_end(); ++i)
	{
        if (!i->is_bisector() || i->is_inner_bisector()) continue;
		auto start = i->vertex()->point();
		auto end = i->opposite()->vertex()->point();
        auto d = sqrt(CGAL::squared_distance(start,end));
        if (d > first) {
            first = d;
        }
    }

    // You can pass the polygon via an iterator pair
    PolygonPtr ss = CGAL::create_offset_polygons_2<Polygon_2> (first / divisor,skeleton)[0];
    //CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*iss);
    std::vector<double> coordinates;
    coordinates.assign(buffer,0);
    for(const Point_2& p : ss->vertices()){
        coordinates.push_back(CGAL::to_double(p.x()));
        coordinates.push_back(CGAL::to_double(p.y()));
    }
    OffsetPolygon result;
    result.coordinates = &coordinates[0];
    result.length = coordinates.size();
    return result;
}

struct Point {
    double x;
    double y;
    int id;
};

struct Edge {
    int start;
    int end;
};

struct StraightSkeleton {
    Point* points;
    int num_points;
    Edge* edges;
    int num_edges;
};

StraightSkeleton unweighted_skeleton(unsigned long outer_count, double *outer_points, unsigned long hole_count, double* hole_points)
{
    Polygon_2 polygon;
    Polygon_2 hole;
    for(int i = buffer; i + 1 < outer_count; i+= 2) {
        polygon.push_back(Point_2(outer_points[i], outer_points[i+1]));
    }
    for(int i = buffer; i + 1 < hole_count; i+= 2) {
        hole.push_back(Point_2(hole_points[i], hole_points[i+1]));
    }
    assert(polygon.is_counterclockwise_oriented());
    Polygon_with_holes poly( polygon ) ;
    if (hole_count > 0) {
        assert(hole.is_clockwise_oriented());
        poly.add_hole( hole ) ;
    }
    // You can pass the polygon via an iterator pair
    SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly);
    Ss& ss = *iss;
    std::vector<Point> points;
    points.assign(buffer,Point());
    //CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*iss);
    for (Vertex_const_iterator i = ss.vertices_begin(); i != ss.vertices_end(); ++i)
	{
		auto current = i->point();
        Point p;
        p.x = current.x();
        p.y = current.y();
        p.id = i->id();
        points.push_back(p);
    }
    std::vector<Edge> edges;
    edges.assign(buffer,Edge());
    for (Halfedge_const_iterator i = ss.halfedges_begin(); i != ss.halfedges_end(); ++i)
	{
		auto start = i->vertex()->id();
		auto end = i->opposite()->vertex()->id();
        if (!i->is_bisector() || start > end) continue;
        Edge e;
        e.start = start;
        e.end = end;
        edges.push_back(e);
    }
    StraightSkeleton result;
    result.points = points.data();
    result.num_points = points.size();
    result.edges = edges.data();
    result.num_edges = edges.size();
    return result;
}

StraightSkeleton weighted_skeleton(unsigned long outer_count, double *outer_points) {
    Polygon_2 polygon;
    std::vector<std::vector<FT> > weights(1);
    for(int i = buffer; i < outer_count; i+= 3) {
        polygon.push_back(Point_2(outer_points[i], outer_points[i+1]));
        weights[0].push_back(outer_points[i+2]);
    }

    assert(polygon.is_counterclockwise_oriented());
    Polygon_with_holes poly( polygon ) ;
    // You can pass the polygon via an iterator pair
    SsPtr iss = CGAL::create_interior_weighted_straight_skeleton_2(poly,weights);
    Ss& ss = *iss;
    std::vector<Point> points;
    points.assign(buffer,Point());
    //CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*iss);
    for (Vertex_const_iterator i = ss.vertices_begin(); i != ss.vertices_end(); ++i)
	{
		auto current = i->point();
        Point p;
        p.x = current.x();
        p.y = current.y();
        p.id = i->id();
        points.push_back(p);
    }
    std::vector<Edge> edges;
    edges.assign(buffer,Edge());
    for (Halfedge_const_iterator i = ss.halfedges_begin(); i != ss.halfedges_end(); ++i)
	{
		auto start = i->vertex()->id();
		auto end = i->opposite()->vertex()->id();
        if (!i->is_bisector() || start > end) continue;
        Edge e;
        e.start = start;
        e.end = end;
        edges.push_back(e);
    }
    StraightSkeleton result;
    result.points = points.data();
    result.num_points = points.size();
    result.edges = edges.data();
    result.num_edges = edges.size();
    return result;
}