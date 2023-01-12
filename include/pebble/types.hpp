/*!
  \file types.hpp
  \brief Output structures of the library.

  All lines traced by the library on the surface of the input mesh is associated
to the halfedge containing the saddle point it originates from. Specific
sections of the curves are associated to the face or edge they are running on as
well. For example, a saddle contour line is a sequence of connected circular
arcs, each arc running on a single face of the mesh. Each such arc can be
associated to the halfedge containing the saddle point the contour originates
from and also to the face it is running on.

The above mentioned associations could be made bidirectional, but the library
doesn't store them as such. Each face and halfedge has a collection of contour
line sections running through them, and these sections reference the saddle the
contour line originated from. For ascending curves, the direction of the
references is opposite, the saddle knows all the curves originating from it, and
the intividual segments of these curves reference the edges or faces they are
running through.
*/

#ifndef PEBBLE_TYPES_HPP
#define PEBBLE_TYPES_HPP 1

#include <boost/graph/graph_traits.hpp>
#include <variant>
#include <iostream>

namespace pebble {
  /*!
    \class FaceContour
    \brief Portion of a contour line running through a face of the input mesh.

    A portion of a contour line running through a face of the input mesh is always
    a circular arc. The center of the arc is the centroid projected to the plane
    of the face, this point is not stored in this structure. Only the arc's
    endpoint and angle is stored.
    \tparam MeshHalfedge the halfedge descriptor of the input mesh
    \tparam Point a model of Point_3
    \tparam Scalar a floating point type
  */
  template <typename MeshHalfedge,
	    typename Point,
	    typename Scalar>
  struct FaceContour {
    //! The halfedge that has the saddle this contour belongs to.
    MeshHalfedge origin;
    //! The point the arc originates from.
    Point from;
    //! The angle of the arc.
    Scalar angle;
    //! Assigning constructor.
    FaceContour(const MeshHalfedge &origin,
		const Point &from,
		const Scalar &angle)
      : origin(origin), from(from), angle(angle) {
    }
  };
  /*!
    \class HalfedgeContour
    \brief Stores where a contour line intersects a given halfedge.
  */
  template <typename MeshHalfedge,
	    typename Point>
  struct HalfedgeContour {
    //! The halfedge that has the saddle this contour belongs to.
    MeshHalfedge origin;
    //! The point where the halfedge intersects a contour line.
    Point at;
    //! Assigning constructor.
    HalfedgeContour(const MeshHalfedge &origin,
		    const Point &at)
      : origin(origin), at(at) {
    }
  };

  /*!
    \class ascending_segment
    \brief A segment of an ascending curve.
    \tparam Mesh a model of FaceGraph
    \tparam Segment a model of Segment_3
   */
  template <typename Mesh,
	    typename Segment>
  class ascending_segment {
    using MeshHalfedge = typename boost::graph_traits<Mesh>::halfedge_descriptor;
    using MeshFace = typename boost::graph_traits<Mesh>::face_descriptor;

    //! The halfedge or face this segment is running through.
    std::variant<MeshHalfedge, MeshFace> m_simplex;
    //! The geometry of the segment.
    Segment m_segment;
  public:
    //! \brief Constructs a segment on a halfedge.
    ascending_segment(const MeshHalfedge &halfedge,
		      const Segment &segment)
      : m_simplex(halfedge), m_segment(segment) {
    }
    //! \brief Constructs a segment on a face.
    ascending_segment(const MeshFace &face,
		      const Segment &segment)
      : m_simplex(face), m_segment(segment) {
    }
    //! \return returns whether the segment runs on a face
    bool on_face() const {
      return m_simplex.index() != 0;
    }
    //! \return returns the halfedge the segment is running on
    MeshHalfedge halfedge() const {
      return get<MeshHalfedge>(m_simplex);
    }
    //! \return returns the face the segment is running on
    MeshFace face() const {
      return get<MeshFace>(m_simplex);
    }
    //! \return returns the geometry of the segment
    Segment segment() const {
      return m_segment;
    }
  };

  /*!
    \class Intersection
    \brief An intersection of a saddle contour and an ascending curve.
    \tparam MeshHalfedge the halfedge descriptor of the input mesh
    \tparam Point a model of Point_3
  */
  template <typename MeshHalfedge,
	    typename Point>
  struct Intersection {
    //! The halfedge having the saddle whose saddle contour line goes through this intersection.
    MeshHalfedge origin;
    //! Position of the intersection.
    Point at;
  };

  //! Type of a vertex in the master graph.
  enum class VertexType { MIN, MAX, SADDLE, INTERSECTION };
  // inline std::ostream& operator<<(std::ostream &stream, const VertexType &value) {
  //   switch(value) {
  //   case VertexType::MIN:
  //     stream << 'S';
  //     break;
  //   case VertexType::MAX:
  //     stream << 'U';
  //     break;
  //   case VertexType::SADDLE:
  //     stream << 'H';
  //     break;
  //   case VertexType::INTERSECTION:
  //     stream << 'I';
  //     break;
  //   }
  //   return stream;
  // }
  
  //! Type of an edge in the master graph.
  enum class EdgeType {ASCENDING, CONTOUR};
  // inline std::ostream& operator<<(std::ostream &stream, const EdgeType &value) {
  //   switch(value) {
  //   case EdgeType::ASCENDING:
  //     stream << 'a';
  //     break;
  //   case EdgeType::CONTOUR:
  //     stream << 'c';
  //     break;
  //   }
  //   return stream;
  // }
}

#endif // PEBBLE_TYPES_HPP
