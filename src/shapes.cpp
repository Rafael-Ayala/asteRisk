// #include <Rcpp.h>
// #include <boost/lexical_cast.hpp>
// #include <list>
// #include <cassert>
// 
// #include <CGAL/Arrangement_2.h>
// #include <CGAL/Polyhedron_3.h>
// #include <CGAL/Polyhedron_traits_with_normals_3.h>
// #include <CGAL/IO/Polyhedron_iostream.h>
// #include <CGAL/IO/Gps_iostream.h>
// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
// #include <CGAL/Boolean_set_operations_2.h>
// #include <CGAL/Polygon_set_2.h>
// 
// typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
// typedef CGAL::Polyhedron_traits_with_normals_3<Kernel> Polyhedron_traits;
// typedef CGAL::Polyhedron_3<Polyhedron_traits>       Polyhedron;
// typedef Kernel::Point_3                             Point_3;
// typedef Kernel::Plane_3                             Plane_3;
// typedef Kernel::Direction_3                         Direction_3;
// typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
// typedef Kernel::Point_2                            Point;
// typedef CGAL::Polygon_2<Kernel>                    Polygon;
// typedef CGAL::Polygon_with_holes_2<Kernel>         Polygon_with_holes;
// typedef std::list<Polygon_with_holes>              Pgn_with_holes_container;
// typedef CGAL::Polygon_set_2<Kernel>             Polygon_set;
// typedef std::list<Polygon_with_holes>            Pgn_with_holes_2_container;
// 
// using namespace Rcpp;
// 
// struct Normal_equation{
//     template <typename Facet> typename Facet::Plane_3 operator()(Facet & f){
//         typename Facet::Halfedge_handle h = f.halfedge();
//         return CGAL::cross_product(h->next()->vertex()->point() - h->vertex()->point(), h->next()->next()->vertex()->point() - h->next()->vertex()->point());
//     }
// };
// 
// // [[Rcpp::plugins("cpp14")]]
// // [[Rcpp::export]]
// int testProj(int argc, std::string argv){
//     // Read the direction from the command line.
//     CGAL_assertion(argc > 3);
//     Kernel::FT x = boost::lexical_cast<double>(argv[1]);
//     Kernel::FT y = boost::lexical_cast<double>(argv[2]);
//     Kernel::FT z = boost::lexical_cast<double>(argv[3]);
//     Direction_3 direction(x, y, z);
//     
//     // Read the polyhedron from the specified input file.
//     // const char* filename = (argc > 4) ? argv[4] : "hand.off";
//     std::string filename = argv;
//     std::ifstream  in_file(filename);
//     if (!in_file.is_open()) {
//         std::cerr << "Failed to open " << filename << "!" << std::endl;
//         return -1;
//     }
//     Polyhedron polyhedron;
//     in_file >> polyhedron;
//     std::transform(polyhedron.facets_begin(), polyhedron.facets_end(),
//                    polyhedron.planes_begin(), Normal_equation());
//     
//     // Go over the polyhedron facets and project them onto the plane.
//     std::list<Polygon> polygons;
//     Kernel kernel;
//     Kernel::Compare_z_3 cmp_z = kernel.compare_z_3_object();
//     Kernel::Construct_projected_xy_point_2 proj =
//         kernel.construct_projected_xy_point_2_object();
//     Kernel::Construct_translated_point_3 translate =
//         kernel.construct_translated_point_3_object();
//     Point_3 origin = kernel.construct_point_3_object()(CGAL::ORIGIN);
//     Plane_3 plane =  kernel.construct_plane_3_object()(origin, direction);
//     Point_3 r = translate(origin, direction.vector());
//     Polyhedron::Facet_const_iterator fit;
//     for (fit = polyhedron.facets_begin(); fit != polyhedron.facets_end(); ++fit) {
//         // Discard facets facing the given direction.
//         if (CGAL::angle(translate(origin, fit->plane()), origin, r) == CGAL::OBTUSE)
//             continue;
//         
//         // Go over the facet vertices and project them.
//         Polygon polygon;
//         Polyhedron::Halfedge_around_facet_const_circulator hit = fit->facet_begin();
//         do {
//             const Point_3& point = hit->vertex()->point();
//             polygon.push_back(proj(plane, point));
//         } while (++hit != fit->facet_begin());
//         polygons.push_back(polygon);
//     }
//     polyhedron.clear();
//     Polygon_set S;
//     S.join(polygons.begin(), polygons.end());
//     std::cout << S;
//     return 0;
// };
