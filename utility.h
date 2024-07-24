#ifndef UTILITY_H
#define UTILITY_H
#include <vector>
#include <cmath>

class Utility{

public:
    typedef std::pair<float, float> Point;
    static std::vector<Point> genInterpolatedPoints(const std::vector<Point>& points,
                                             int num);
    static std::vector<Point> genParallelLine(const std::vector<Point>& points, float offset);
private:
    static Point catmullRom(Point p0, Point p1, Point p2, Point p3, float t);
    static void intersectionLineDetect(std::vector<Utility::Point>& parallel_points);
    static void adjustPoints(std::vector<Utility::Point>& parallel_points,
                                     int separated_point_num = 0);
    static inline float distance(const Point& p1, const Point& p2){
        float dx = p2.first - p1.first;
        float dy = p2.second - p2.second;
        return std::sqrt(dx * dx + dy * dy);
    }
    static inline Point parallelPoint(const Point& p, const Point& normal, float distance){
        return {p.first + normal.first * distance, p.second + normal.second * distance};
    }
    static inline Point computeNormal(const Point& p1, const Point& p2){
        float dx = p2.first - p1.first;
        float dy = p2.second - p1.second;
        float length = std::sqrt(dx * dx + dy * dy);
        return {-dy / length, dx / length };
    }
    static inline Point midpoint(const Point& p1, const Point& p2){
        return {(p1.first + p2.first) / 2, (p1.second + p2.second) / 2};
    }
    static inline float cross(const Point& p1, const Point& p2, const Point& p3){
        return (p2.first - p1.first) * (p3.second - p1.second) -
               (p2.second - p1.second) * (p3.first - p1.first);
    }
    static inline bool segmentsIntersect(const Point& p1,
                                         const Point& p2,
                                         const Point& p3,
                                         const Point& p4)
    {
        return (cross(p1, p3, p4) * cross(p2, p3, p4) <= 0 ) &&
               (cross(p3, p1, p2) * cross(p4, p1, p2) <= 0);
    }
    static inline Point computeIntersection(const Point& p1,
                                            const Point& p2,
                                            const Point& p3,
                                            const Point& p4)
    {
        float t = cross(p1, p3, p4) / (cross(p1, p3, p4) - cross(p2, p3, p4));
        return {p1.first + t * (p2.first - p1.first),
                p1.second + t * (p2.second - p1.second)};
    }
};


#endif // UTILITY_H
