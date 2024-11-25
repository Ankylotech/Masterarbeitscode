#include<vector>
#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Polygon_with_holes_2.h>
#include<CGAL/create_offset_polygons_2.h>
#include<CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>
#include<CGAL/create_weighted_straight_skeleton_from_polygon_with_holes_2.h>
#include<CGAL/Polygon_2.h>

int buffer = 4;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;

typedef K::Point_2                    Point;
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

int main() {
    Polygon_2 polygon;
    polygon.push_back(Point(8.66667,-4.33333));
    polygon.push_back(Point(17.33333,-4.33333));
    polygon.push_back(Point(15.1667,-2.16667));
    polygon.push_back(Point(8.66667,4.33333));
    assert(polygon.is_counterclockwise_oriented());
    Polygon_with_holes poly( polygon ) ;
    // You can pass the polygon via an iterator pair
    SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly);
    Ss& ss = *iss;
    //CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*iss);
    for (Vertex_const_iterator i = ss.vertices_begin(); i != ss.vertices_end(); ++i)
	{
		auto current = i->point();
        std::cout << current.x() << "," << current.y() << " " << i->id() << std::endl;
    }
    for (Halfedge_const_iterator i = ss.halfedges_begin(); i != ss.halfedges_end(); ++i)
	{
		auto start = i->vertex()->id();
		auto end = i->opposite()->vertex()->id();
        if (!i->is_bisector() || start > end) continue;
        
        std::cout << start << "," << end  << std::endl;
    }
}