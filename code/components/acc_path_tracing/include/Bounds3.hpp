#ifndef __BOUNDS3_HPP__
#define __BOUNDS3_HPP__

#include <vector>
#include <list>
#include <algorithm>
#include <array>
#include "scene/Scene.hpp"
#include "Ray.hpp"

using namespace NRenderer;
namespace AccPathTracer {
    class Bounds3 {
    public:
        Vec3 min;
        Vec3 max;
        enum class Type
        {
            SPHERE = 0x1,
            TRIANGLE = 0X2,
            PLANE = 0X3,
            MESH = 0X4
        };
        Bounds3() {
            double minNum = std::numeric_limits<double>::lowest();
            double maxNum = std::numeric_limits<double>::max();
            max = Vec3(minNum, minNum, minNum);
            min = Vec3(maxNum, maxNum, maxNum);
        }
        Bounds3(Sphere* sp) {
            type = Type::SPHERE;
            this->sp = sp;
            float r = sp->radius;
            min = Vec3{ sp->position.x - r, sp->position.y - r, sp->position.z - r };
            max = Vec3{ sp->position.x + r, sp->position.y + r, sp->position.z + r };
        }
        Bounds3(Triangle* tr) {
            type = Type::TRIANGLE;
            this->tr = tr;
            min = Vec3{ std::min(tr->v1.x, std::min(tr->v2.x, tr->v3.x)),
                std::min(tr->v1.y, std::min(tr->v2.y, tr->v3.y)),
                std::min(tr->v1.z, std::min(tr->v2.z, tr->v3.z)) };

            max = Vec3{ std::max(tr->v1.x, std::max(tr->v2.x, tr->v3.x)),
                std::max(tr->v1.y, std::max(tr->v2.y, tr->v3.y)),
                std::max(tr->v1.z, std::max(tr->v2.z, tr->v3.z)) };
        }
        Vec3 Centroid() {
            return 0.5f * min + 0.5f* max;
        }
        Vec3 Diagonal() const { return max - min; }
        inline bool IntersectP(const Ray& ray, const Vec3& invDir,
            const std::array<int, 3>& dirisNeg);
        int maxExtent() const
        {
            Vec3 d = Diagonal();
            if (d.x > d.y && d.x > d.z)
                return 0;
            else if (d.y > d.z)
                return 1;
            else
                return 2;
        }
        union
        {
            Sphere* sp;
            Triangle* tr;
        };
        Type type ;

    };
    inline bool Bounds3::IntersectP(const Ray& ray, const Vec3& invDir, const std::array<int, 3>& dirIsNeg) {
        float t1_min, t1_max, t2_min, t2_max, t3_min, t3_max;
        t1_min = (min[0] - ray.origin[0]) * invDir[0];
        t2_min = (min[1] - ray.origin[1]) * invDir[1];
        t3_min = (min[2] - ray.origin[2]) * invDir[2];
        t1_max = (max[0] - ray.origin[0]) * invDir[0];
        t2_max = (max[1] - ray.origin[1]) * invDir[1];
        t3_max = (max[2] - ray.origin[2]) * invDir[2];
        /*if (!dirIsNeg[0]) std::swap(t1_max, t1_min);
        if (!dirIsNeg[1]) std::swap(t2_max, t2_min);
        if (!dirIsNeg[2]) std::swap(t3_max, t3_min);*/
        float t_near = std::max(t1_min, std::max(t2_min, t3_min));
        float t_far = std::min(t1_max, std::min(t2_max, t3_max));
        if (t_near > 0 && t_far > 0 && t_near < t_far) {
            return true;
        }
        return false;
    }
}

#endif