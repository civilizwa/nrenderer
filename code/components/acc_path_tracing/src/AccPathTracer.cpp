#include "server/Server.hpp"

#include "AccPathTracer.hpp"

#include "VertexTransformer.hpp"
#include "intersections/intersections.hpp"

#include "glm/gtc/matrix_transform.hpp"

const auto taskNums = 32;
AccPathTracer::Timer timers[taskNums]{};

namespace AccPathTracer
{

    RGB AccPathTracerRenderer::gamma(const RGB& rgb) {
        return glm::sqrt(rgb);
    }

    void AccPathTracerRenderer::renderTask(RGBA* pixels, int width, int height, int off, int step) {
        for (int i = off; i < height; i += step) {
            for (int j = 0; j < width; j++) {
                Vec3 color{ 0, 0, 0 };
                for (int k = 0; k < samples; k++) {
                    auto r = defaultSamplerInstance<UniformInSquare>().sample2d(); // 随机生成采样点
                    float rx = r.x;
                    float ry = r.y;
                    float x = (float(j) + rx) / float(width);
                    float y = (float(i) + ry) / float(height);
                    auto ray = camera.shoot(x, y);
                    color += trace(ray, 0, off);
                }
                color /= samples;
                color = gamma(color);
                pixels[(height - i - 1) * width + j] = { color, 1 };
            }
        }
        

    }

    auto AccPathTracerRenderer::render() -> RenderResult {
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
        Model m = spScene->models[0];
        //在这里生成Bounds
        if (acc_type == 1) {
            getBox();
            tree = new BVHTree(spScene);
            tree->root = tree->build(box);
        }
        if (acc_type == 2) {
            // kd_tree.Insert(scene.sphereBuffer, scene.triangleBuffer);
        }
        //-------------
        thread t[taskNums];
        for (int i = 0; i < taskNums; i++) {
            t[i] = thread(&AccPathTracerRenderer::renderTask,
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

    void AccPathTracerRenderer::release(const RenderResult& r) {
        auto [p, w, h] = r;
        delete[] p;
    }

    HitRecord AccPathTracerRenderer::closestHitObject(const Ray& r) {
        HitRecord closestHit = nullopt;
        float closest = FLOAT_INF;

        // BVH
        auto hitRecord = tree->Intersect(r, closest);
        if (hitRecord && hitRecord->t < closest ) {
            closest = hitRecord->t;
            closestHit = hitRecord;
        }
        //for (auto& s : scene.sphereBuffer) {
        //    auto hitRecord = Intersection::xSphere(r, s, 0.000001, closest);
        //    if (hitRecord && hitRecord->t < closest) {
        //        closest = hitRecord->t;
        //        closestHit = hitRecord;
        //    }
        //}
        //for (auto& t : scene.triangleBuffer) {
        //    auto hitRecord = Intersection::xTriangle(r, t, 0.000001, closest);
        //    if (hitRecord && hitRecord->t < closest) {
        //        closest = hitRecord->t;
        //        closestHit = hitRecord;
        //    }
        //}
        //for (auto& p : scene.planeBuffer) {
        //    auto hitRecord = Intersection::xPlane(r, p, 0.000001, closest);
        //    if (hitRecord && hitRecord->t < closest) {
        //        closest = hitRecord->t;
        //        closestHit = hitRecord;
        //    }
        //}
           
        return closestHit;

        // KD-tree
        //else if(acc_type==2)
        //{
        //    auto result = kd_tree.find_Node(r);
        //    for (auto& n : result) {
        //        for (auto& b : n->Box_list) {
        //            HitRecord hitRecord;
        //            if (b->type == Bounds3::Type::SPHERE) hitRecord = Intersection::xSphere(r, *(b->sp), 0.01, closest);
        //            else if (b->type == Bounds3::Type::TRIANGLE) hitRecord = Intersection::xTriangle(r, *(b->tr), 0.01, closest);
        //            if (hitRecord && hitRecord->t < closest) {
        //                closest = hitRecord->t;
        //                closestHit = hitRecord;
        //            }
        //        }
        //    }
        //}
    }

    tuple<float, Vec3> AccPathTracerRenderer::closestHitLight(const Ray& r) {
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
    void AccPathTracerRenderer::getBox() {
        BVHTree* tree = new BVHTree(this->spScene);
        BVHNode* temp = new BVHNode();
        vector<Node> objects = spScene->nodes;
        temp->buildBounds(spScene, objects);
        box = temp->bounds;
    }

    RGB AccPathTracerRenderer::trace(const Ray& r, int currDepth, int thread_id) {
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