#include "Metropolis.hpp"
#include "server/Server.hpp"
#include "VertexTransformer.hpp"
#include "intersections/intersections.hpp"
#include "glm/gtc/matrix_transform.hpp"

const auto taskNums = 32;
Metropolis::Timer timers[taskNums]{};

namespace Metropolis {
    RGB MetropolisRenderer::gamma(const RGB& rgb) {
        return glm::sqrt(rgb);
    }

    void MetropolisRenderer::renderTask(RGBA* pixels, int width, int height, int off, int step) {
        // TODO 应该就改改这里
        float d = rnd();
 
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

    Vec3 MetropolisRenderer::VecRandom(const float rnd1, const float rnd2)
    {
        const float temp1 = 2.0 * PI * rnd1, temp2 = 2.0 * rnd2 - 1.0;				// temp1 in [0, 2pi), temp2 in [-1, 1)
        const float s = sin(temp1), c = cos(temp1), t = sqrt(1.0 - temp2 * temp2); // s, c确定两个方位角，t确定维度
        return Vec3(s * t, temp2, c * t);
    }

    Vec3 MetropolisRenderer::VecCosine(const Vec3 n, const float g, const float rnd1, const float rnd2)
    {
        // g->inf, temp2->1-, t->0+, 余弦分布就越均匀；反之g = 1，则余弦分布越容易得到法线方向
        const float temp1 = 2.0 * PI * rnd1, temp2 = pow(rnd2, 1.0 / (g + 1.0));
        const float s = sin(temp1), c = cos(temp1), t = sqrt(1.0 - temp2 * temp2);
        Onb base = Vec3(s * t, temp2, c * t);
        return base.local(n);
    }

    Vec3 MetropolisRenderer::AccumulatePathContribution(const PathContribution pc, const float mScaling)
    {
        if (pc.sc == 0)
            return Vec3{ 0 };
        for (int i = 0; i < pc.n; i++)
        {
            const int ix = int(pc.c[i].x), iy = int(pc.c[i].y); // 通过顶点找到其贡献的照片的像素位置
            Vec3 c = pc.c[i].c * mScaling;
            if ((ix < 0) || (ix >= width) || (iy < 0) || (iy >= height))
                continue;
            c += c;
        }
    }


}