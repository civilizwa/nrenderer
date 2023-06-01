#pragma once
#ifndef __CAMERA_HPP__
#define __CAMERA_HPP__

#include "scene/Camera.hpp"
#include "geometry/vec.hpp"

#include "samplers/SamplerInstance.hpp"

#include "Ray.hpp"

namespace Metropolis
{
    using namespace std;
    using namespace NRenderer;
    class Camera
    {
    private:


    public:
        const NRenderer::Camera& camera;
        float lenRadius; // 光圈半径
        float theta;
        Vec3 u, v, w; // w向量表示相机的观察方向，u和v向量与w垂直，用于表示相机的水平和垂直方向
        Vec3 vertical; // 相机的垂直方向
        Vec3 horizontal; // 相机的水平方向
        Vec3 lowerLeft; // 视平面左下角的位置
        Vec3 position; // 相机的位置
        float halfHeight;
        // float dist = 0.0f; // 相机到成像平面的距离

        Camera(const NRenderer::Camera& camera)
            : camera                (camera)
        {
            position = camera.position;
            lenRadius = camera.aperture / 2.f;
            auto vfov = camera.fov;
            vfov = clamp(vfov, 160.f, 20.f);
            theta = glm::radians(vfov);
            halfHeight = tan(theta/2.f);
            auto halfWidth = camera.aspect*halfHeight;
            Vec3 up = camera.up;
            w = glm::normalize(camera.position - camera.lookAt);
            u = glm::normalize(glm::cross(up, w));
            v = glm::cross(w, u);

            auto focusDis = camera.focusDistance; // 焦距

            lowerLeft = position - halfWidth*focusDis*u - halfHeight*focusDis*v - focusDis*w; // 相机位置减去成像平面左下角的位置
            horizontal = 2*halfWidth*focusDis*u;
            vertical = 2*halfHeight*focusDis*v;
            // dist =  / (2.0 * tan((camera.fov / 2.0) * (PI / 180.0)));
            
            
        }

        // 从摄像机中发射光线
        Ray shoot(float s, float t) const {
            auto r = defaultSamplerInstance<UniformInCircle>().sample2d();
            float rx = r.x * lenRadius;
            float ry = r.y * lenRadius;
            Vec3 offset = u*rx + v*ry;
            return Ray{
                position + offset,
                glm::normalize(
                    lowerLeft + s*horizontal + t*vertical - position - offset
                )
            };
        }
    };
}

#endif