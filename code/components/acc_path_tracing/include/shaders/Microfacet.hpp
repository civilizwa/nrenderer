#pragma once
#ifndef __MICROFACET_HPP__
#define __MICROFACET_HPP__

#include "Shader.hpp"

namespace AccPathTracer
{
    class Microfacet : public Shader
    {
    private:
        Vec3 albedo;
        RenderOption renderoption;
    public:
        Microfacet(Material& material, vector<Texture>& textures,RenderOption& renderopt);
        Scattered shade(const Ray& ray, const Vec3& hitPoint, const Vec3& normal) const;
    };
}

#endif