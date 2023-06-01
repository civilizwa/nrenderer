#pragma once
#ifndef __SHADER_CREATOR_HPP__
#define __SHADER_CREATOR_HPP__

#include "Shader.hpp"
#include "Lambertian.hpp"
#include "Microfacet.hpp"
#include "Glass.hpp"

namespace AccPathTracer
{
    class ShaderCreator
    {
    public:
        ShaderCreator() = default;
        SharedShader create(Material& material, vector<Texture>& t, RenderOption& renderoption) {
            SharedShader shader{nullptr};
            if (renderoption.shaderType == 0) {
                switch (material.type)
                {
                case 0:
                    shader = make_shared<Lambertian>(material, t);
                    break;
                case 2:
                    shader = make_shared<Glass>(material, t);
                    break;
                default:
                    shader = make_shared<Lambertian>(material, t);
                    break;
                }
            }
            else if (renderoption.shaderType == 1) {
                shader = make_shared<Microfacet>(material, t, renderoption);//微表面
            }
            
            return shader;
        }
    };
}

#endif