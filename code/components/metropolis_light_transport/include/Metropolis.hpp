#pragma once
#ifndef __METROPOLIS_HPP__
#define __METROPOLIS_HPP__

#include "scene/Scene.hpp"
#include "Ray.hpp"
#include "Camera.hpp"
#include "intersections/HitRecord.hpp"
#include "Timer.hpp"
#include "shaders/ShaderCreator.hpp"
#include "Bounds3.hpp"
#include "BVH.hpp"
#include "Onb.hpp"
#include "TKahanAdder.hpp"
#include "PathContribution.hpp"
#include "TMarkovChian.hpp"

#include <tuple>


namespace Metropolis
{


    using namespace NRenderer;
    using namespace std;

 
    class MetropolisRenderer
    {
    private:
        unsigned long sample_num = 0;
        unsigned long mutations = 10;
        int PathRndsOffset = 0;
        float prnds[NumStates];

        SharedScene spScene;
        Scene& scene;

        unsigned int width;
        unsigned int height;
        unsigned int depth;
        unsigned int samples;
        unsigned int acc_type;

        using SCam = Metropolis::Camera;
        SCam camera;

        vector<SharedShader> shaderPrograms;
    public:
        MetropolisRenderer(SharedScene spScene)
            : spScene(spScene)
            , scene(*spScene)
            , camera(spScene->camera)
        {
            width = scene.renderOption.width;
            height = scene.renderOption.height;
            depth = scene.renderOption.depth;
            samples = scene.renderOption.samplesPerPixel;
            acc_type = scene.renderOption.acc_type;

            // for BVH
            getBox();
            tree = new BVHTree(spScene);
            tree->root = tree->build(box);
        }
        ~MetropolisRenderer() = default;

        using RenderResult = tuple<RGBA*, unsigned int, unsigned int>;
        RenderResult render();
        void release(const RenderResult& r);

        // BVH
        vector<Bounds3> box;
        BVHTree* tree;
        void getBox();

        // 用两个随机数 in [0, 1)^2 生成一个终点在半球面上的向量
        Vec3 VecRandom(const float rnd1, const float rnd2);

        // 构造一个按余弦分布的随机方向，更倾向于与法线接近的方向
        Vec3 VecCosine(const Vec3 n, const float g, const float rnd1, const float rnd2);

        // 将float映射到[0, 255]，并进行了gamma校正
        int toInt(float x) { return int(pow(1 - exp(-x), 1 / 2.2) * 255 + .5); }

        Vec3 AccumulatePathContribution(const PathContribution pc, const float mScaling);

        // xorshift PRNG，生成[0, 1) float
        inline float rnd() {
            static unsigned int x = 123456789, y = 362436069, z = 521288629, w = 88675123;
            unsigned int t = x ^ (x << 11);
            x = y;
            y = z;
            z = w;
            return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8))) * (1.0 / 4294967296.0);
        }

        // 根据s1, s2扰动value
        inline float perturb(const float value, const float s1, const float s2)
        {
            float Result;
            float r = rnd();
            if (r < 0.5)
            {
                r = r * 2.0;
                Result = value + s2 * exp(-log(s2 / s1) * r);
                if (Result > 1.0)
                    Result -= 1.0;
            }
            else
            {
                r = (r - 0.5) * 2.0;
                Result = value - s2 * exp(-log(s2 / s1) * r);
                if (Result < 0.0)
                    Result += 1.0;
            }
            return Result;
        }

        // large_step就是生成一条与之前完全无关的光路
        TMarkovChain large_step(TMarkovChain MC)
        {
            TMarkovChain Result;
            Result.C = MC.C;
            for (int i = 0; i < NumStates; i++)
                Result.u[i] = rnd();
            return Result;
        }

        TMarkovChain mutate(TMarkovChain MC)
        {
            TMarkovChain Result;
            Result.C = MC.C;

            // pixel location
            Result.u[0] = perturb(MC.u[0], 2.0 / float(width + height), 0.1);
            Result.u[1] = perturb(MC.u[1], 2.0 / float(width + height), 0.1);

            // the rest
            for (int i = 2; i < NumStates; i++)
                Result.u[i] = perturb(MC.u[i], 1.0 / 1024.0, 1.0 / 64.0);
            return Result;
        }

        void InitRandomNumbersByChain(const TMarkovChain MC)
        {
            for (int i = 0; i < NumStates; i++)
                prnds[i] = MC.u[i];
        }
        void InitRandomNumbers()
        {
            for (int i = 0; i < NumStates; i++)
                prnds[i] = rnd();
        }

        // local sampling PDFs and standard terms
        inline float GeometryTerm(const Vert e0, const Vert e1)
        {
            const Vec3 dv = e1.p - e0.p;
            const float d2 = glm::dot(dv, dv);
            return fabs(glm::dot(e0.n, dv) * glm::dot(e0.n, dv)) / (d2 * d2);
        }
        inline float DirectionToArea(const Vert current, const Vert next)
        {
            const Vec3 dv = next.p - current.p;
            const float d2 = glm::dot(dv, dv);
            return fabs(glm::dot(next.n, dv)) / (d2 * sqrt(d2));
        }
        inline float GlossyBRDF(const Vec3 wi, const Vec3 n, const Vec3 wo)
        {
            const float won = glm::dot(wo, n);
            const float win = glm::dot(wi, n);
            const Vec3 r = Reflect(-wi, n);
            return (Glossiness + 2.0) / (2.0 * PI) * pow(MAX(glm::dot(r, wo), 0.0), Glossiness) / MAX(fabs(win), fabs(won));
        }
        inline float GlossyPDF(const Vec3 wi, const Vec3 n, const Vec3 wo)
        {
            const Vec3 r = Reflect(-wi, n);
            return (Glossiness + 1.0) / (2.0 * PI) * pow(MAX(glm::dot(r, wo), 0.0), Glossiness);
        }
        inline float LambertianBRDF(const Vec3 wi, const Vec3 n, const Vec3 wo)
        {
	        return 1.0 / PI;
        }
        inline float LambertianPDF(const Vec3 wi, const Vec3 n, const Vec3 wo)
        {
	        return fabs(glm::dot(wo, n)) / PI;
        }

        // measurement contribution function
        Vec3 PathThroughput(const Path Xb)
        {
            Vec3 f = Vec3(1.0, 1.0, 1.0);
            for (int i = 0; i < Xb.n; i++)
            {
                if (i == 0)
                {
                    float W = 1.0 / float(width * height);
                    Vec3 d0 = Xb.x[1].p - Xb.x[0].p;
                    const float dist2 = glm::dot(d0, d0);
                    d0 = d0 * (1.0f / sqrt(dist2));
                    const float c = glm::dot(d0, Camera.w);
                    const float ds2 = (Camera.dist / c) * (Camera.dist / c);
                    W = W / (c / ds2);
                    f = f * (W * fabs(glm::dot(d0, Xb.x[1].n) / dist2));
                }
                else if (i == (Xb.n - 1))
                {
                    if (sph[Xb.x[i].id].refl == LGHT)
                    {
                        const Vec3 d0 = glm::normalize((Xb.x[i - 1].p - Xb.x[i].p));
                        const float L = LambertianBRDF(d0, Xb.x[i].n, d0);
                        f = f. * (sph[Xb.x[i].id].c * L);
                    }
                    else
                    {
                        f = f * 0.0f;
                    }
                }
                else
                {
                    const Vec3 d0 = glm::normalize((Xb.x[i - 1].p - Xb.x[i].p));
                    const Vec3 d1 = glm::normalize((Xb.x[i + 1].p - Xb.x[i].p));
                    float BRDF = 0.0;
                    if (sph[Xb.x[i].id].refl == DIFF)
                    {
                        BRDF = LambertianBRDF(d0, Xb.x[i].n, d1);
                    }
                    else if (sph[Xb.x[i].id].refl == GLOS)
                    {
                        BRDF = GlossyBRDF(d0, Xb.x[i].n, d1);
                    }
                    f = f * (sph[Xb.x[i].id].c * BRDF * GeometryTerm(Xb.x[i], Xb.x[i + 1]));
                }
                if (MAX(MAX(f.x, f.y), f.z) == 0.0)
                    return f;
            }
            return f;
        }
    private:
        void renderTask(RGBA* pixels, int width, int height, int off, int step);

        RGB gamma(const RGB& rgb);
        RGB trace(const Ray& ray, int currDepth, int thread_id);
        HitRecord closestHitObject(const Ray& r);
        tuple<float, Vec3> closestHitLight(const Ray& r);

    };



}

#endif