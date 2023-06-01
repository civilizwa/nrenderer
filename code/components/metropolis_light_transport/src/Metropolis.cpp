#include "Metropolis.hpp"
#include "server/Server.hpp"
#include "VertexTransformer.hpp"
#include "intersections/intersections.hpp"
#include "glm/gtc/matrix_transform.hpp"

const auto taskNums = 1;
Metropolis::Timer timers[taskNums]{};

namespace Metropolis {
    RGB MetropolisRenderer::gamma(const RGB& rgb) {
        return glm::sqrt(rgb);
    }

    void MetropolisRenderer::renderTask(RGBA* pixels, int width, int height, int off, int step) {
        unsigned long samps = 0; // 采样数

        // 调试CombinePaths->PathProbablityDensity中

        // estimate normalization constant b应该是用双向路径追踪计算照片最后的平均亮度
        float b = 0.0;
        for (int i = 0; i < N_Init; i++) {
            if (i == 60) {
                cout << "here" << endl;
            }

            cout << i << endl;
            InitRandomNumbers();
            Path eyePath = GenerateEyePath(MaxEvents);
            Path lightPath = GenerateLightPath(MaxEvents);
            PathContribution pc = CombinePaths(eyePath, lightPath);
            b += pc.sc;
        }
        b /= float(N_Init);

        // initialize the Markov chain
        TMarkovChain current, proposal;
        InitRandomNumbersByChain(current);
        current.C = CombinePaths(GenerateEyePath(MaxEvents), GenerateLightPath(MaxEvents));

        // integration
        for (int i = 0; i < mutations; i++) {
            samps++;
            // sample the path
            float isLargeStepDone;
            if (rnd() <= LargeStepProb)
            {
                proposal = large_step(current);
                isLargeStepDone = 1.0;
            }
            else
            {
                proposal = mutate(current);
                isLargeStepDone = 0.0;
            }
            InitRandomNumbersByChain(proposal);
            proposal.C = CombinePaths(GenerateEyePath(MaxEvents), GenerateLightPath(MaxEvents));

            float a = 1.0; // 接受proposal的概率
            if (current.C.sc > 0.0)
                a = MAX(MIN(1.0, proposal.C.sc / current.C.sc), 0.0);

            // accumulate samples
            if (proposal.C.sc > 0.0)
                AccumulatePathContribution(pixels, proposal.C, (a + isLargeStepDone) / (proposal.C.sc / b + LargeStepProb));
            if (current.C.sc > 0.0)
                AccumulatePathContribution(pixels, current.C, (1.0 - a) / (current.C.sc / b + LargeStepProb));

            // update the chain
            if (rnd() <= a)
                current = proposal;
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
        HitRecord closest = getHitRecord(FLOAT_INF, {}, {}, {}, {});
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

    void MetropolisRenderer::trace(Path &path, const Ray& r, int currDepth) {
        if (currDepth == depth) return;

        auto hitObject = closestHitObject(r);
        auto [t, emitted] = closestHitLight(r);

        if (hitObject && hitObject->t < t) {
            //auto mtlHandle = hitObject->material;
            //auto scattered = shaderPrograms[mtlHandle.index()]->shade(r, hitObject->hitPoint, hitObject->normal);
            //auto scatteredRay = scattered.ray;
            //auto attenuation = scattered.attenuation;
            //auto emitted = scattered.emitted;

            // 把相交点信息添加进path
            Vec3 p = hitObject->hitPoint, n = hitObject->normal;
            n = glm::dot(n, r.direction) < 0 ? n : n * -1.0f; // 让n与r方向相反
            path.x[path.n] = Vert(p, n, hitObject->id);
            path.n++;
            const float rnd0 = prnds[(currDepth - 1) * NumRNGsPerEvent + 0 + PathRndsOffset];
            const float rnd1 = prnds[(currDepth - 1) * NumRNGsPerEvent + 1 + PathRndsOffset];

            // 对从交点发出的散射光线进行渲染
            Ray nr;
            nr.origin = p;
            nr.direction = VecCosine(n, 1.0, rnd0, rnd1); // 漫反射材质的反射光方向

            trace(path, nr, currDepth + 1);

            //auto next = trace(path, nr, currDepth + 1);
            //float n_dot_in = glm::dot(hitObject->normal, scatteredRay.direction);
            //float pdf = scattered.pdf;

        }
        else if (t != FLOAT_INF) {
            cout << "int trace(), hit the light source." << endl;
            path.x[path.n] = Vert(hitObject->hitPoint, hitObject->normal, light_id);
            path.n++;
        }
    }

    // 沿单位球面分布的向量
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



}