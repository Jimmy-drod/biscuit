#include "utility.h"
#include <list>
#include <iostream>
#define PARALLEL_START_OFFSET (0)
#define PARALLEL_END_OFFSET PARALLEL_START_OFFSET
std::vector<Utility::Point>
Utility::genInterpolatedPoints
(const std::vector<Utility::Point>& points, int num){
    std::vector<Utility::Point> interpolated_points;
    for(int i = 0; i < (int)(points.size()) - 1; ++i){
        const auto& p0 = (i == 0) ? points[0] : points[i-1];
        const auto& p1 = points[i];
        const auto& p2 = points[i+1];
        const auto& p3 = (i == (int)(points.size()) - 2) ?
                         points[i+1] : points[i+2];

        for(int j = 0; j < num; ++j){
            float t = static_cast<float>(j) / num;
            auto interpolated_point = catmullRom(p0, p1, p2, p3, t);
            interpolated_points.push_back(interpolated_point);
        }
    }
    intersectionLineDetect(interpolated_points);
    return interpolated_points;
}

std::vector<Utility::Point>
Utility::genParallelLine
(const std::vector<Utility::Point>& points, float offset)
{
    std::vector<Utility::Point> parallelLine;
    for(size_t i = 0; i < points.size() - 1; ++i){
        const Utility::Point& p1 = points[i];
        const Utility::Point& p2 = points[i+1];
        Utility::Point normal = computeNormal(p1, p2);
        Utility::Point parallelPoint1 = parallelPoint(p1, normal, offset);
        Utility::Point parallelPoint2 = parallelPoint(p2, normal, offset);
        parallelLine.push_back(parallelPoint1);
        parallelLine.push_back(parallelPoint2);
    }

    intersectionLineDetect(parallelLine);
    return parallelLine;
}

Utility::Point Utility::catmullRom(Utility::Point p0,
                                   Utility::Point p1,
                                   Utility::Point p2,
                                   Utility::Point p3,
                                   float t)
{
    float t2 = t * t;
    float t3 = t2 * t;

    float b0 = 0.5f * (-t3 + 2.0f * t2 - t);
    float b1 = 0.5f * (3.0f * t3 - 5.0f * t2 + 2.0f);
    float b2 = 0.5f * (-3.0f * t3 + 4.0f * t2 + t);
    float b3 = 0.5f * (t3 - t2);

    float x = p0.first * b0 + p1.first *b1 + p2.first * b2 + p3.first * b3;
    float y = p0.second * b0 + p1.second *b1 + p2.second * b2 + p3.second * b3;

    return {x, y};
}

void
Utility::intersectionLineDetect
(std::vector<Utility::Point>& parallel_points)
{
    int max_separated_point_num = parallel_points.size() >= 4 ?
                parallel_points.size() - 4 : 0;
    size_t last_size = parallel_points.size();
    size_t update_size = last_size;
    for(int i = 0; i <= max_separated_point_num; i++ ){
        if(last_size == update_size){
            if(parallel_points.size() - 4 >= (size_t)(i)){
                adjustPoints(parallel_points, i);
                update_size = parallel_points.size();
            }
            else{
                break;
            }
        }
        else{
            i = 0;
            last_size = update_size;
            max_separated_point_num = parallel_points.size() >= 4 ?
                            parallel_points.size() - 4 : 0;
            if(parallel_points.size() - 4 >= (size_t)(i)){
                adjustPoints(parallel_points, i);
                update_size = parallel_points.size();
            }
            else{
                break;
            }
        }
    }
}

void
Utility::adjustPoints
(std::vector<Utility::Point>& parallel_points, int separated_point_num)
{
    if(separated_point_num < 0 || parallel_points.size() < 4) return;

    std::list<Utility::Point> points_list;
    for(auto& point : parallel_points){
        points_list.push_back(point);
    }

    auto it = points_list.begin();
    while(std::next(it, 3 + separated_point_num) != points_list.end()){
        auto p1 = *it;
        auto p2 = *std::next(it);
        auto p3 = *std::next(it, 2 + separated_point_num);
        auto p4 = *std::next(it, 3 + separated_point_num);
        if(segmentsIntersect(p1, p2, p3, p4)){
            for(int i = 0 ; i < separated_point_num + 2; i++){
                auto delete_it = std::next(it, 1);
                if(delete_it != points_list.end()){
                    points_list.erase(delete_it);
                }
                else {
                    break;
                }
            }
            auto intersection = computeIntersection(p1, p2, p3, p4);
            if(std::next(it, 1) != points_list.end())
                points_list.insert(std::next(it, 1), intersection);
        }
        else{
            ++it;
        }
    }

    parallel_points.clear();
    std::copy(points_list.begin(), points_list.end(),
                    std::back_inserter(parallel_points));
}
