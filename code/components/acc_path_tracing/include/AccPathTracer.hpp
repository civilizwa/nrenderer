#pragma once
#ifndef __SIMPLE_PATH_TRACER_HPP__
#define __SIMPLE_PATH_TRACER_HPP__

#include "scene/Scene.hpp"
#include "Ray.hpp"
#include "Camera.hpp"
#include "intersections/HitRecord.hpp"

#include "shaders/ShaderCreator.hpp"
#include "Bounds3.hpp"
#include "BVH.hpp"
#include <tuple>
#include <Bounds3.hpp>
#include "KD.hpp"

namespace AccPathTracer
{
    using namespace NRenderer;
    using namespace std;

    class AccPathTracerRenderer
    {
    public:
    private:
        SharedScene spScene;
        Scene& scene;

        unsigned int width;
        unsigned int height;
        unsigned int depth;
        unsigned int samples;
        unsigned int acc_type;

        using SCam = AccPathTracer::Camera;
        SCam camera;

        vector<SharedShader> shaderPrograms;
    public:
        AccPathTracerRenderer(SharedScene spScene)
            : spScene(spScene)
            , scene(*spScene)
            , camera(spScene->camera)
        {
            width = scene.renderOption.width;
            height = scene.renderOption.height;
            depth = scene.renderOption.depth;
            samples = scene.renderOption.samplesPerPixel;
            acc_type = scene.renderOption.acc_type;
            if (acc_type == 1) {
                getBox();
                tree = new BVHTree(spScene);
                tree->root = tree->build(box);
            }
            if (acc_type == 2) {
                kd_tree.Insert(scene.sphereBuffer, scene.triangleBuffer);
            }
        }
        ~AccPathTracerRenderer() = default;

        using RenderResult = tuple<RGBA*, unsigned int, unsigned int>;
        RenderResult render();
        void release(const RenderResult& r);
        vector<Bounds3> box;
        BVHTree* tree;
        KDTree kd_tree;
        void getBox();

    private:
        void renderTask(RGBA* pixels, int width, int height, int off, int step);

        RGB gamma(const RGB& rgb);
        RGB trace(const Ray& ray, int currDepth);
        HitRecord closestHitObject(const Ray& r);
        tuple<float, Vec3> closestHitLight(const Ray& r);

    };
}

#endif