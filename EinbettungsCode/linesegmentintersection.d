mixin(import("aloft.all"));
import std.stdio : writeln, write;
import std.math.algebraic: abs;
import std.algorithm.comparison : min,max;
import std.conv;
import std.math.rounding: floor,ceil, round;
import std.math.exponential: pow;
// https://github.com/williamfiset/Algorithms
auto EPS = 1e-7;
// Finds the orientation of point 'c' relative to the line segment (a, b)
// Returns  0 if all three points are collinear.
// Returns -1 if 'c' is clockwise to segment (a, b), i.e right of line formed by the segment.
// Returns +1 if 'c' is counter clockwise to segment (a, b), i.e left of line
// formed by the segment.
int orientation(Apair a, Apair b, Apair c) {
    real value = (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y);
    if (abs(value) < EPS) return 0;
    return (value > 0) ? -1 : +1;
  }
// Tests whether point 'c' is on the line segment (a, b).
// Ensure first that point c is collinear to segment (a, b) and
// then check whether c is within the rectangle formed by (a, b)
bool pointOnLine(Apair a, Apair b, Apair c) {
  return orientation(a, b, c) == 0 && 
         min(a.x, b.x) <= c.x + EPS && c.x <= max(a.x, b.x) + EPS && 
         min(a.y, b.y) <= c.y + EPS && c.y <= max(a.y, b.y) + EPS;
}
// Determines whether two segments intersect.
bool segmentsIntersect(Apair p1, Apair p2, Apair p3, Apair p4) {
  // first easy case: they cannot intersect as they are in completely different locations
  if (max(p1.x,p2.x) < min(p3.x,p4.x) || max(p3.x,p4.x) < min(p1.x,p2.x) || 
    max(p1.y,p2.y) < min(p3.y,p4.y) || max(p3.y,p4.y) < min(p1.y,p2.y)) {
    return false;
  }
  // Get the orientation of points p3 and p4 in relation
  // to the line segment (p1, p2)
  int o1 = orientation(p1, p2, p3);
  int o2 = orientation(p1, p2, p4);
  int o3 = orientation(p3, p4, p1);
  int o4 = orientation(p3, p4, p2);
  // If the points p1, p2 are on opposite sides of the infinite
  // line formed by (p3, p4) and conversly p3, p4 are on opposite
  // sides of the infinite line formed by (p1, p2) then there is
  // an intersection.
  if (o1 != o2 && o3 != o4) return true;
  // Collinear special cases (perhaps these if checks can be simplified?)
  if (o1 == 0 && pointOnLine(p1, p2, p3)) return true;
  if (o2 == 0 && pointOnLine(p1, p2, p4)) return true;
  if (o3 == 0 && pointOnLine(p3, p4, p1)) return true;
  if (o4 == 0 && pointOnLine(p3, p4, p2)) return true;
  return false;
}
Apair[] getCommonEndpoints(Apair p1, Apair p2, Apair p3, Apair p4) {
  Apair[] points = [];
  if (p1.eq(p3)) {
    points ~= p1;
    if (p2.eq(p4)) points ~= p2;
  } else if (p1.eq(p4)) {
    points ~= p1;
    if (p2.eq(p3)) points ~= p2;
  } else if (p2.eq(p3)) {
    points ~= p2;
  } else if (p2.eq(p4)) {
    points ~= p2;
  }
  return points;
}
// Finds the intersection point(s) of two line segments. Unlike regular line 
// segments, segments which are points (x1 = x2 and y1 = y2) are allowed.
Apair[] lineSegmentLineSegmentIntersection(Apair p1, Apair p2, Apair p3, Apair p4) {
  // No intersection.
  if (!segmentsIntersect(p1, p2, p3, p4)) return [];
  // Both segments are a single point.
  if (p1.eq(p2) && p2.eq(p3) && p3.eq(p4))
    return [p1];
  Apair[] endpoints = getCommonEndpoints(p1, p2, p3, p4);

  ulong n = endpoints.length;
  // One of the line segments is an intersecting single point.
  // NOTE: checking only n == 1 is insufficient to return early
  // because the solution might be a sub segment.
  bool singleton = p1.eq(p2) || p3.eq(p4);
  if ((n == 1 && singleton) || n == 2) return endpoints;

  bool collinearSegments = (orientation(p1, p2, p3) == 0) && 
                              (orientation(p1, p2, p4) == 0);
  // The intersection will be a sub-segment of the two
  // segments since they overlap each other. 
  if (collinearSegments) {
    // Segment #2 is enclosed in segment #1
    if (pointOnLine(p1, p2, p3) && pointOnLine(p1, p2, p4))
      return [p3, p4];
    // Segment #1 is enclosed in segment #2
    if (pointOnLine(p3, p4, p1) && pointOnLine(p3, p4, p2))
      return [p1, p2];
    // The subsegment is part of segment #1 and part of segment #2.
    // Find the middle points which correspond to this segment.
    Apair midPoint1 = pointOnLine(p1, p2, p3) ? p3.dup : p4.dup;
    Apair midPoint2 = pointOnLine(p3, p4, p1) ? p1.dup : p2.dup;
    // There is actually only one middle point!
    if (midPoint1.eq(midPoint2)) return [midPoint1];
    return [midPoint1, midPoint2];
  }
  if (pointOnLine(p1,p2,p3)) {
    return [p3];
  }
  if (pointOnLine(p1,p2,p4)) {
    return [p4];
  }
  if (pointOnLine(p3,p4,p1)) {
    return [p1];
  }
  if (pointOnLine(p3,p4,p2)) {
    return [p2];
  }
  /* Beyond this point there is a unique intersection point. */
  // Segment #1 is a vertical line.
  if (abs(p1.x - p2.x) < EPS) {
    real m = (p4.y - p3.y) / (p4.x - p3.x);
    real b = p3.y - m * p3.x;
    auto result = Apair(p1.x, m * p1.x + b);
    return [result];
  }
  // Segment #2 is a vertical line.
  if (abs(p3.x - p4.x) < EPS) {
    real m = (p2.y - p1.y) / (p2.x - p1.x);
    real b = p1.y - m * p1.x;
    auto result = Apair(p3.x, m * p3.x + b);
    return [result];
  }
  real m1 = (p2.y - p1.y) / (p2.x - p1.x);
  real m2 = (p4.y - p3.y) / (p4.x - p3.x);
  real b1 = p1.y - m1 * p1.x;
  real b2 = p3.y - m2 * p3.x;
  real x = (b2 - b1) / (m1 - m2);
  real y = (m1 * b2 - m2 * b1) / (m1 - m2);
  auto result = Apair(x, y);
  return [result];
}
