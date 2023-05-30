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
#include <algorithm>
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

        // W：我用了一种更简单的方法去判断某物体的类型是不是光源，你在此文件中搜索light_id就能找到了
        ////该函数由gpt生成
        //bool pointInAreaLight(const AreaLight& a, const Vec3& point, const float tMin, const float tMax) {
        //    Vec3 normal = glm::cross(a.u, a.v);
        //    Vec3 position = a.position;
        //    auto Np_dot_d = glm::dot(normal, point - position);
        //    if (Np_dot_d < 0.0000001f && Np_dot_d > -0.00000001f) {
        //        return false;
        //    }
        //    float dp = -glm::dot(position, normal);
        //    float t = (-dp - glm::dot(normal, point - position)) / Np_dot_d;
        //    if (t >= tMax || t < tMin) {
        //        return false;
        //    }
        //    // cross test    
        //    Vec3 hitPoint = point;
        //    Mat3x3 d{ a.u, a.v, glm::cross(a.u, a.v) };
        //    d = glm::inverse(d);
        //    auto res = d * (hitPoint - position);
        //    auto u = res.x, v = res.y;
        //    if ((u <= 1 && u >= 0) && (v <= 1 && v >= 0)) {
        //        return true;
        //    }
        //    return false;
        //}
        // 
  
        // measurement contribution function
        Vec3 PathThroughput(const Path Xb)
        {
            Vec3 f = Vec3(1.0, 1.0, 1.0);
            for (int i = 0; i < Xb.n; i++)
            {
                // 获取color的思路：用id找到它在nodes中的实体类型(plane/triangle/sphere)，然后可以通过这个实体的material.index()获取，然后通过if语句结合.scn文件中的数值得到颜色
                int color_index = -1;
                Vec3 color{ 0 };
                int id = Xb.x[i].id;
                auto node = scene.nodes[id];
                if (node.type == Node::Type::SPHERE) {
                    color_index = scene.sphereBuffer[node.entity].material.index();
                }
                else if (node.type == Node::Type::TRIANGLE) {
                    color_index = scene.triangleBuffer[node.entity].material.index();
                }
                else if (node.type == Node::Type::PLANE) {
                    color_index = scene.planeBuffer[node.entity].material.index();
                }
                else if (node.type == Node::Type::MESH) {
                    color_index = scene.meshBuffer[node.entity].material.index();
                }

                if (color_index == 0)
                    color = { 0.725, 0.71, 0.68 };
                else if (color_index == 1)
                    color = { 0.63, 0.065, 0.05 };
                else if (color_index == 2)
                    color = { 0.14, 0.45, 0.091 };

                if (i == 0)
                {
                    float W = 1.0 / float(width * height);
                    Vec3 d0 = Xb.x[1].p - Xb.x[0].p;
                    const float dist2 = glm::dot(d0, d0);
                    d0 = d0 * (1.0f / sqrt(dist2));
                    const float c = glm::dot(d0, camera.w);
                    float dist = width / (2.0 * camera.halfHeight);
                    const float ds2 = (dist / c) * (dist / c);
                    W = W / (c / ds2);
                    f = f * (W * fabs(glm::dot(d0, Xb.x[1].n) / dist2));
                }
                else if (i == (Xb.n - 1))
                {
                    if (Xb.x[i].id == light_id) // 这样写应该可以吧，判断一个点的reflect类型是不是光源的方法就是去看它的id
                    {
                        const Vec3 d0 = glm::normalize((Xb.x[i - 1].p - Xb.x[i].p));
                        const float L = LambertianBRDF(d0, Xb.x[i].n, d0);

                        f = f * (color * L);
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
                    BRDF = LambertianBRDF(d0, Xb.x[i].n, d1);
                    f = f * (color * BRDF * GeometryTerm(Xb.x[i], Xb.x[i + 1]));
                }
                if (Max(f) == 0.0)
                    return f;
            }
            return f;
        }

        // check if the path can be connected or not (visibility term)
        bool isConnectable(const Path Xeye, const Path Xlight, float& px, float& py)
        {
            // TODO
            Vec3 Direction;
            // 取出两条路线的最后一个点
            const Vert& Xeye_e = Xeye.x[Xeye.n - 1];
            const Vert& Xlight_e = Xlight.x[Xlight.n - 1];

            bool Result;
            if ((Xeye.n == 0) && (Xlight.n >= 2))
            {
                // no direct hit to the film (pinhole)
                Result = false;
            }
            else if ((Xeye.n >= 2) && (Xlight.n == 0))
            {
                // direct hit to the light source
                Result = (Xeye_e.id == light_id);
                Direction = glm::normalize((Xeye.x[1].p - Xeye.x[0].p));
            }
            else if ((Xeye.n == 1) && (Xlight.n >= 1))
            {
                // light tracing
                Ray r(Xeye.x[0].p, glm::normalize((Xlight_e.p - Xeye.x[0].p)));
                auto hitObject = closestHitObject(r);
                Result = (hitObject && (hitObject->id == Xlight_e.id));
                Direction = r.direction;
            }
            else
            {
                // shadow ray connection
                Ray r(Xeye_e.p, glm::normalize((Xlight_e.p - Xeye_e.p)));
                auto hitObject = closestHitObject(r);
                Result = (hitObject && (hitObject->id == Xlight_e.id));
                Direction = glm::normalize((Xeye.x[1].p - Xeye.x[0].p));
            }

            // get the pixel location
            float dist = width / (2.0 * camera.halfHeight);
            Vec3 ScreenCenter = camera.position + (camera.w * dist);
            Vec3 ScreenPosition = camera.position + (Direction * (dist / glm::dot(Direction, camera.w))) - ScreenCenter;
            px = glm::dot(camera.u, ScreenPosition) + (width * 0.5);
            py = glm::dot(-camera.v, ScreenPosition) + (height * 0.5);
            return Result && ((px >= 0) && (px < height) && (py >= 0) && (py < height));
        }

        // path probability density
        // - take the sum of all possible probability densities if the numbers of subpath vertices are not specified
        float PathProbablityDensity(const Path SampledPath, const int PathLength, const int SpecifiedNumEyeVertices = -1, const int SpecifiedNumLightVertices = -1)
        {
            TKahanAdder SumPDFs(0.0);
            bool Specified = (SpecifiedNumEyeVertices != -1) && (SpecifiedNumLightVertices != -1);

            float dist = width / (2.0 * camera.halfHeight);
            // number of eye subpath vertices
            for (int NumEyeVertices = 0; NumEyeVertices <= PathLength + 1; NumEyeVertices++)
            {
                // extended BPT
                float p = 1.0;

                // number of light subpath vertices
                int NumLightVertices = (PathLength + 1) - NumEyeVertices;

                // we have pinhole camera
                if (NumEyeVertices == 0)
                    continue;

                // add all?
                if (Specified && ((NumEyeVertices != SpecifiedNumEyeVertices) || (NumLightVertices != SpecifiedNumLightVertices)))
                    continue;

                // sampling from the eye
                for (int i = -1; i <= NumEyeVertices - 2; i++)
                {
                    if (i == -1)
                    {
                        // PDF of sampling the camera position (the same delta function with the scaling 1.0 for all the PDFs - they cancel out)
                        p = p * 1.0;
                    }
                    else if (i == 0)
                    {
                        p = p * 1.0 / float(width * height);
                        Vec3 Direction0 = glm::normalize(SampledPath.x[1].p - SampledPath.x[0].p);
                        float CosTheta = glm::dot(Direction0, camera.w);
                        float DistanceToScreen2 = dist / CosTheta;
                        DistanceToScreen2 = DistanceToScreen2 * DistanceToScreen2;
                        p = p / (CosTheta / DistanceToScreen2);

                        p = p * DirectionToArea(SampledPath.x[0], SampledPath.x[1]);
                    }
                    else
                    {
                        // PDF of sampling ith vertex
                        Vec3 Direction0 = glm::normalize(SampledPath.x[i - 1].p - SampledPath.x[i].p);
                        Vec3 Direction1 = glm::normalize(SampledPath.x[i + 1].p - SampledPath.x[i].p);
                        p = p * LambertianPDF(Direction0, SampledPath.x[i].n, Direction1);
                        p = p * DirectionToArea(SampledPath.x[i], SampledPath.x[i + 1]);
                    }
                }

                if (p != 0.0)
                {
                    // sampling from the light source
                    for (int i = -1; i <= NumLightVertices - 2; i++)
                    {
                        if (i == -1)
                        {
                            // PDF of sampling the light position (assume area-based sampling)
                            float lightArea = Max(scene.areaLightBuffer[0].u) * Max(scene.areaLightBuffer[0].v);
                            p = p * (1.0 / lightArea);
                        }
                        else if (i == 0)
                        {
                            Vec3 Direction0 = glm::normalize(SampledPath.x[PathLength - 1].p - SampledPath.x[PathLength].p);
                            p = p * LambertianPDF(SampledPath.x[PathLength].n, SampledPath.x[PathLength].n, Direction0);
                            p = p * DirectionToArea(SampledPath.x[PathLength], SampledPath.x[PathLength - 1]);
                        }
                        else
                        {
                            // PDF of sampling (PathLength - i)th vertex
                            Vec3 Direction0 = glm::normalize(SampledPath.x[PathLength - (i - 1)].p - SampledPath.x[PathLength - i].p);
                            Vec3 Direction1 = glm::normalize(SampledPath.x[PathLength - (i + 1)].p - SampledPath.x[PathLength - i].p);
                            p = p * LambertianPDF(Direction0, SampledPath.x[PathLength - i].n, Direction1);
                            p = p * DirectionToArea(SampledPath.x[PathLength - i], SampledPath.x[PathLength - (i + 1)]);
                        }
                    }
                }

                if (Specified && (NumEyeVertices == SpecifiedNumEyeVertices) && (NumLightVertices == SpecifiedNumLightVertices))
                    return p;

                // sum the probability density (use Kahan summation algorithm to reduce numerical issues)
                SumPDFs.add(p);
            }
            return SumPDFs.sum;
        }

        Ray SampleLightSources(const float rnd1, const float rnd2)
        {
            // base1, base2是面光源的对角线上的两个点
            auto& areaLight = scene.areaLightBuffer[0];
            Vec3 base1 = areaLight.position + areaLight.u;
            Vec3 base2 = areaLight.position + areaLight.v;
            float x_len = std::fabs(base1.x - base2.x), z_len = std::fabs(base1.z - base2.z); // 考虑到面光源与xz平面平行，就这么写了
            Vec3 p = base1 + rnd1 * x_len + rnd2 * z_len;

            // 用沿单位球面分布的向量对y加绝对值取负得到沿单位半球（xz平面以下的半球）分布的向量
            Vec3 d = VecRandom(rnd1, rnd2);
            if (d.y > 0)
                d.y *= -1.0f;
            return Ray{p, d};
        }
        Path GenerateLightPath(const int MaxLightEvents)
        {
            // 初始化路径采样结果Result
            Path Result;
            Result.n = 0;
            if (MaxLightEvents == 0)
                return Result;
            for (int i = 0; i < MaxEvents; i++)
                Result.x[i].id = -1;
            PathRndsOffset = NumStatesSubpath;

            // 从光源采样一条光线
            Ray r = SampleLightSources(prnds[PathRndsOffset + 0], prnds[PathRndsOffset + 1]);
            PathRndsOffset += NumRNGsPerEvent;

            // Result第一个点放光源的信息
            auto& areaLight = scene.areaLightBuffer[0];
            Vec3 n = glm::normalize(glm::cross(areaLight.u, areaLight.v));
            Result.x[0] = Vert(r.origin, n, light_id);
            Result.n++;

            // 对这条光线去采样生成完整路径
            trace(Result, r, 1);
            return Result;
        }

        Ray SampleCamera(const float rnd1, const float rnd2)
        {
            // 不用写，应该可以用camera.shoot替代
            return Ray{};
        }
        Path GenerateEyePath(const int MaxEyeEvents)
        {
            // 初始化路径采样结果Result
            Path Result;
            Result.n = 0;
            if (MaxEyeEvents == 0)
                return Result;
            // 初始化Result中的顶点
            for (int i = 0; i < MaxEvents; i++) {
                Result.x[i].id = -1; // -1 表示该顶点还没有东西，另外light_id = -3, camera_id = -2
            }
            PathRndsOffset = 0;

            // 从camera采样一束光线
            Ray r = camera.shoot(prnds[PathRndsOffset + 0], prnds[PathRndsOffset + 1]); // 这里应该可以用这个函数替代SampleCamera()吧
            PathRndsOffset += NumRNGsPerEvent;

            // Result第一个点放camera信息
            Result.x[0] = Vert(r.origin, camera.w, camera_id);
            Result.n++;

            // 对这条光线去采样生成完整路径
            trace(Result, r, 1);
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
            Result.n = 0;
            Result.sc = 0.0;

            for (int PathLength = MinPathLength; PathLength <= MaxPathLength; PathLength++) {
                for (int NumEyeVertices = 0; NumEyeVertices <= PathLength + 1; NumEyeVertices++)
                {
                    const int NumLightVertices = (PathLength + 1) - NumEyeVertices;

                    if (NumEyeVertices == 0)
                        continue; // no direct hit to the film (pinhole) 成像平面？
                    if (NumEyeVertices > EyePath.n)
                        continue;
                    if (NumLightVertices > LightPath.n)
                        continue;

                    // extract subpaths
                    Path Eyesubpath = EyePath;
                    Path Lightsubpath = LightPath;
                    Eyesubpath.n = NumEyeVertices;
                    Lightsubpath.n = NumLightVertices;

                    // check the path visibility
                    float px = -1.0, py = -1.0;
                    if (!isConnectable(Eyesubpath, Lightsubpath, px, py))
                        continue;

                    // construct a full path
                    Path SampledPath;
                    for (int i = 0; i < NumEyeVertices; i++)
                        SampledPath.x[i] = EyePath.x[i];
                    for (int i = 0; i < NumLightVertices; i++)
                        SampledPath.x[PathLength - i] = LightPath.x[i];
                    SampledPath.n = NumEyeVertices + NumLightVertices;

                    // evaluate the path
                    Vec3 f = PathThroughput(SampledPath);
                    float p = PathProbablityDensity(SampledPath, PathLength, NumEyeVertices, NumLightVertices);
                    float w = MISWeight(SampledPath, NumEyeVertices, NumLightVertices, PathLength);
                    if ((w <= 0.0) || (p <= 0.0))
                        continue;

                    Vec3 c = f * (w / p);
                    if (Max(c) <= 0.0)
                        continue;

                    // store the pixel contribution
                    Result.c[Result.n] = Contribution(px, py, c);
                    Result.n++;

                    // scalar contribution function
                    Result.sc = MAX(Max(c), Result.sc);
                }
            }

            return Result;
        }

        void AccumulatePathContribution(RGBA* pixels, const PathContribution pc, const float mScaling)
        {
            // TODO
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

    private:
        void renderTask(RGBA* pixels, int width, int height, int off, int step);

        RGB gamma(const RGB& rgb);
        void trace(Path &path, const Ray& ray, int currDepth);
        HitRecord closestHitObject(const Ray& r);
        tuple<float, Vec3> closestHitLight(const Ray& r);

    };



}

#endif