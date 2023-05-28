#include "Metropolis.hpp"
#include "server/Server.hpp"
#include "VertexTransformer.hpp"
#include "intersections/intersections.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "Raytracer.hpp"

const auto taskNums = 32;
Metropolis::Timer timers[taskNums]{};

namespace Metropolis {
    RGB MetropolisRenderer::gamma(const RGB& rgb) {
        return glm::sqrt(rgb);
    }

    void MetropolisRenderer::renderTask(RGBA* pixels, int width, int height, int off, int step) {
        // TODO 应该就改改这里
        std::vector<Vec3> tmp_image;
        tmp_image.resize(width * height, Vec3{ 0.0 });

        MLT mlt;
        int seed = rand();
        mlt.xor128_.setSeed(off + seed);

        float b;
        float p_large = 0.1;
        int acceptedPaths = 0, rejectPaths = 0;
        PathSample old_path;

        mlt.large_step = 1;

        // 这个循环找到一条有效的初始路径old_path
        while (1)
        {
            mlt.ResetRandomCoords();
            //Compute new path

            PathSample sample = generate_new_path(camera, info.m_width, info.m_height, mlt);
            mlt.global_time++;

            //Clear the stack
            while (!mlt.primary_samples_stack.empty())
                mlt.primary_samples_stack.pop();

            auto value = luminance(sample.Color);
            //Check if valid path
            if (value > 0.0f)
            {
                b = value;
                p_large = 0.5f;
                old_path = sample;
                break;
            }
        }
    }

    auto MetropolisRenderer::render() -> RenderResult {
        for (int i = 0; i < taskNums; i++) {
            timers[i].init();
        }
        // shaders
        shaderPrograms.clear();
        ShaderCreator shaderCreator{};
        for (auto& m : scene.materials) {
            shaderPrograms.push_back(shaderCreator.create(m, scene.textures));
        }

        RGBA* pixels = new RGBA[width * height]{};

        // 局部坐标转换成世界坐标
        VertexTransformer vertexTransformer{};
        vertexTransformer.exec(spScene);

        thread t[taskNums];
        for (int i = 0; i < taskNums; i++) {
            t[i] = thread(&MetropolisRenderer::renderTask,
                this, pixels, width, height, i, taskNums);
        }
        for (int i = 0; i < taskNums; i++) {
            t[i].join();
        }
        getServer().logger.log("Done...");

        float total_ms = 0.0;
        for (int i = 0; i < taskNums; i++) {
            total_ms += timers[i].getTime();
            // cout << "thread" << i << ": " << timers[i].getTime() << "ms." << endl;
        }
        cout << "threadNum = " << taskNums << ", closestHitObject time per thread with BVH: " << total_ms / taskNums / 1000.0 << "s." << endl;

        return { pixels, width, height };
    }

    void MetropolisRenderer::release(const RenderResult& r) {
        auto [p, w, h] = r;
        delete[] p;
    }

    HitRecord MetropolisRenderer::closestHitObject(const Ray& r) {
        HitRecord closestHit = nullopt;
        float closest = FLOAT_INF;

        // BVH
        auto hitRecord = tree->Intersect(r, closest);
        if (hitRecord && hitRecord->t < closest) {
            closest = hitRecord->t;
            closestHit = hitRecord;
        }

        return closestHit;
    }

    tuple<float, Vec3> MetropolisRenderer::closestHitLight(const Ray& r) {
        Vec3 v = {};
        HitRecord closest = getHitRecord(FLOAT_INF, {}, {}, {});
        for (auto& a : scene.areaLightBuffer) {
            auto hitRecord = Intersection::xAreaLight(r, a, 0.000001, closest->t);
            if (hitRecord && closest->t > hitRecord->t) {
                closest = hitRecord;
                v = a.radiance;
            }
        }
        return { closest->t, v };
    }
    void MetropolisRenderer::getBox() {
        BVHTree* tree = new BVHTree(this->spScene);
        BVHNode* temp = new BVHNode();
        vector<Node> objects = spScene->nodes;
        temp->buildBounds(spScene, objects);
        box = temp->bounds;
    }

    RGB MetropolisRenderer::trace(const Ray& r, int currDepth, int thread_id) {
        if (currDepth == depth) return scene.ambient.constant;

        timers[thread_id].start();
        auto hitObject = closestHitObject(r);
        timers[thread_id].stop();

        auto [t, emitted] = closestHitLight(r);

        if (hitObject && hitObject->t < t) {
            auto mtlHandle = hitObject->material;
            auto scattered = shaderPrograms[mtlHandle.index()]->shade(r, hitObject->hitPoint, hitObject->normal);
            auto scatteredRay = scattered.ray;
            auto attenuation = scattered.attenuation;
            auto emitted = scattered.emitted;
            auto next = trace(scatteredRay, currDepth + 1, thread_id);
            float n_dot_in = glm::dot(hitObject->normal, scatteredRay.direction);
            float pdf = scattered.pdf;
            /**
             * emitted      - Le(p, w_0)
             * next         - Li(p, w_i)
             * n_dot_in     - cos<n, w_i>
             * atteunation  - BRDF
             * pdf          - p(w)
             **/
            return emitted + attenuation * next * n_dot_in / pdf;
        }
        else if (t != FLOAT_INF) {
            return emitted;
        }
        else {
            return Vec3{ 0 };
        }
    }
}