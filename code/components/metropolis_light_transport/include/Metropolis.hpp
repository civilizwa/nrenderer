#pragma once
#ifndef __SIMPLE_PATH_TRACER_HPP__
#define __SIMPLE_PATH_TRACER_HPP__

#include "scene/Scene.hpp"
#include "Ray.hpp"
#include "Camera.hpp"
#include "intersections/HitRecord.hpp"
#include "Timer.hpp"
#include "shaders/ShaderCreator.hpp"
#include "Bounds3.hpp"
#include "BVH.hpp"
#include <tuple>
#include <Bounds3.hpp>
#include <stack>

namespace Metropolis
{
    using namespace NRenderer;
    using namespace std;

    // 生成xor128随机数的类
    class xor128
    {
    public:
        unsigned int x, y, z, w;
        xor128()
        {
            static unsigned int x = 123456789;
            static unsigned int y = 362436069;
            static unsigned int z = 521288629;
            static unsigned int w = 88675123;
        }
        inline unsigned int step(void)
        {
            unsigned int t;
            t = x ^ (x << 11);
            x = y; y = z; z = w;
            return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
        }
        void setSeed(unsigned u)
        {
            x ^= u;
        }
        inline float rand01() { return static_cast<float>((double)step() / (UINT_MAX)); }
    };

    // 记录生成的随机数信息
    struct PrimarySample
    {
        int modifyTime;
        float value;
        PrimarySample(xor128& valxor)
        {
            modifyTime = 0;
            value = valxor.rand01();
        }
    };

    // 记录采样路径对图片的影响
    struct PathSample
    {
        int x;
        int y;
        glm::vec3 Color;
        float weight;

        PathSample(const int x = 0, const int y = 0, const glm::vec3& color = glm::vec3(0.0f), const float weight = 1.0f) :
            x(x), y(y), Color(color), weight(weight) {}
    };

    PathSample generate_new_path (Camera& camera, const int width, const int height, MLT& mlt) {
        int x, y;
        float weight = 4.0;
        weight *= width;

        x = static_cast<int>(mlt.NextSample() * width);

        if (x == width)
            x = 0;

        weight *= height;

        y = static_cast<int>(mlt.NextSample() * height);

        if (y == height)
            y = 0;

        int halfWidth = static_cast<int>(width * 0.5f);
        int halfHeight = static_cast<int>(height * 0.5f);

        Vec3 viewVec = Camera.v;
    }

    class MLT {

    public:
        xor128 xor128_;
        int global_time;
        int large_step;
        int large_step_time;
        int used_rand_coords;

        vector<PrimarySample> primary_samples;
        stack<PrimarySample> primary_samples_stack;

        MLT() {
            global_time = large_step = large_step_time = used_rand_coords = 0;
            for (int i = 0; i < 128; i++)
                primary_samples.push_back(PrimarySample(xor128_));
        }

        void ResetRandomCoords() { used_rand_coords = 0; }

        // 生成[0, 1]之间的float
        inline double rand01(void) { return xor128_.rand01(); }
        float NextSample();

    private:
        float Mutate(const float  x) {
            const float r = static_cast<float>(xor128_.rand01());
            const float s1 = 1.0f / 512.0f;
            const float	s2 = 1.0f / 16.0f;
            const float dx = s1 / (s1 / s2 + fabsf(2.0f * r - 1.0f)) - s1 / (s1 / s2 + 1.0f);
            if (r < 0.5f)
            {
                const float x1 = x + dx;
                return (x1 < 1.0f) ? x1 : x1 - 1.0f;
            }
            else {
                const float x1 = x - dx;
                return (x1 < 0.0f) ? x1 + 1.f : x1;
            }
        };
    };

    class MetropolisRenderer
    {
    public:
        unsigned long sample_num = 0;
        unsigned long mutations = 10;

    private:
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
        vector<Bounds3> box;
        BVHTree* tree;
        void getBox();

    private:
        void renderTask(RGBA* pixels, int width, int height, int off, int step);

        RGB gamma(const RGB& rgb);
        RGB trace(const Ray& ray, int currDepth, int thread_id);
        HitRecord closestHitObject(const Ray& r);
        tuple<float, Vec3> closestHitLight(const Ray& r);

    };
}

#endif