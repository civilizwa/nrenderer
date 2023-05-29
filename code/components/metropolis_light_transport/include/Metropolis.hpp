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

        void AccumulatePathContribution(RGBA* pixels, const PathContribution pc, const float mScaling)
        {
            if (pc.sc == 0)
                return;
            for (int i = 0; i < pc.n; i++)
            {
                const int ix = int(pc.c[i].x), iy = int(pc.c[i].y); // 通过顶点找到其贡献的照片的像素位置
                Vec3 c = pc.c[i].c * mScaling;
                if ((ix < 0) || (ix >= width) || (iy < 0) || (iy >= height))
                    continue;
                pixels[ix + iy * width] += Vec4{ c, 1 }; // TODO 不确定应该怎么写
                // pixels[(height - ix - 1) * width + iy] += Vec4{ c, 1 };
            }
        }

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
            // TODO
            Vec3 f = Vec3(1.0, 1.0, 1.0);
            return f;
        }

        // check if the path can be connected or not (visibility term)
        bool isConnectable(const Path Xeye, const Path Xlight, float& px, float& py)
        {
            // TODO
            return true;
        }

        // path probability density
        // - take the sum of all possible probability densities if the numbers of subpath vertices are not specified
        float PathProbablityDensity(const Path SampledPath, const int PathLength, const int SpecifiedNumEyeVertices = -1, const int SpecifiedNumLightVertices = -1)
        {
            // TODO
            return 1.0;
        }

        Ray SampleLightSources(const float rnd1, const float rnd2)
        {
            // TODO
            const Vec3 d = VecRandom(rnd1, rnd2);
            return Ray{};
        }
        Path GenerateLightPath(const int MaxLightEvents)
        {
            // TODO

            // 初始化路径采样结果Result
            Path Result;
            return Result;
        }

        Ray SampleCamera(const float rnd1, const float rnd2)
        {
            // TODO
            return Ray{};
        }
        Path GenerateEyePath(const int MaxEyeEvents)
        {
            // TODO

            // 初始化路径采样结果Result
            Path Result;
            return Result;
        }

        // balance heuristic
        float MISWeight(const Path SampledPath, const int NumEyeVertices, const int NumLightVertices, const int PathLength)
        {
            const float p_i = PathProbablityDensity(SampledPath, PathLength, NumEyeVertices, NumLightVertices);
            const float p_all = PathProbablityDensity(SampledPath, PathLength);
            if ((p_i == 0.0) || (p_all == 0.0))
            {
                return 0.0;
            }
            else
            {
                return MAX(MIN(p_i / p_all, 1.0), 0.0);
            }
        }

        // BPT connections 双向路径追踪中的连接
        // - limit the connection to a specific technique if s and t are provided
        PathContribution CombinePaths(const Path EyePath, const Path LightPath, const int SpecifiedNumEyeVertices = -1, const int SpecifiedNumLightVertices = -1)
        {
            // TODO 是一个很重要但是很复杂的函数...
            PathContribution Result;
            return Result;
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