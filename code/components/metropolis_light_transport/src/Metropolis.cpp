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

        int invalid = 0;
        // estimate normalization constant b应该是用双向路径追踪计算照片最后的平均亮度
        double b = 0.0;
        for (int i = 0; i < N_Init; i++) {

            InitRandomNumbers();
            double sc = CombinePaths(GenerateEyePath(MaxEvents), GenerateLightPath(MaxEvents)).sc; // 我需要解决在这里面出现nan的问题
            //cout << "i = " << i << ", sc = " << sc << endl;
            b += sc;

            //Path eyePath = GenerateEyePath(MaxEvents);
            //Path lightPath = GenerateLightPath(MaxEvents);
            //PathContribution pc = CombinePaths(eyePath, lightPath);
            //b += pc.sc;
            //if (eyePath.n + lightPath.n <= 10) {
            //    invalid++;
            //}
            // cout << "i = " << i << ", sumN = " << eyePath.n + lightPath.n << ", sc = " << pc.sc << endl;            
            //if (abs(pc.sc) == DOUBLE_INF) {
            //    cout << "error" << endl;
            //    return;
            //}
        }

        b /= double(N_Init);
        cout << "b = " << b << ", invalid = " << invalid << endl; // b = 0.136614, invalid = 47

        // initialize the Markov chain
        TMarkovChain current, proposal;
        InitRandomNumbersByChain(current);
        current.C = CombinePaths(GenerateEyePath(MaxEvents), GenerateLightPath(MaxEvents));

        // integration
        for (int i = 0; i < mutations; i++) {
            samps++;
            // sample the path
            double isLargeStepDone;
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

            //for (int i = 0; i < current.C.n; i++) {
            //    cout << "i = " << i << ", current.C.c[i].color = " << current.C.c[i].c << endl;
            //}
            //cout << endl << endl;
            //for (int i = 0; i < proposal.C.n; i++) {
            //    cout << "i = " << i << ", proposal.C.c[i].color = " << proposal.C.c[i].c << endl;
            //}
            //cout << endl << endl;

            double a = 1.0; // 接受proposal的概率
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

        renderTask(pixels, width, height, 0, 1);

        //thread t[taskNums];
        //for (int i = 0; i < taskNums; i++) {
        //    t[i] = thread(&MetropolisRenderer::renderTask,
        //        this, pixels, width, height, i, taskNums);
        //}
        //for (int i = 0; i < taskNums; i++) {
        //    t[i].join();
        //}
        getServer().logger.log("Done...");

        for (int ix = 0; ix < width; ix = ix + 10) {
            for (int iy = 0; iy < height; iy = iy + 10) {
                cout << "ix = " << ix << ", iy = " << iy << ", color = " << pixels[(height - iy - 1) * width + ix] << endl;
            }
        }

        return { pixels, width, height };
    }

    void MetropolisRenderer::release(const RenderResult& r) {
        auto [p, w, h] = r;
        delete[] p;
    }

    HitRecord MetropolisRenderer::closestHitObject(const Ray& r) {
        HitRecord closestHit = nullopt;
        double closest = DOUBLE_INF;

        // BVH
        auto hitRecord = tree->Intersect(r, closest);
        if (hitRecord && hitRecord->t < closest) {
            closest = hitRecord->t;
            closestHit = hitRecord;
        }

        return closestHit;
    }

    tuple<double, Vec3d, HitRecord> MetropolisRenderer::closestHitLight(const Ray& r) {
        Vec3d v = {};
        HitRecord closest = getHitRecord(DOUBLE_INF, {}, {}, {}, {});
        for (auto& a : scene.areaLightBuffer) {
            auto hitRecord = Intersection::xAreaLight(r, a, 0.000001, closest->t);
            if (hitRecord && closest->t > hitRecord->t) {
                closest = hitRecord;
                v = a.radiance;
            }
        }
        return { closest->t, v, closest};
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
        auto [t, emitted, hitLight] = closestHitLight(r);

        if (hitObject && hitObject->t < t) {
            //auto mtlHandle = hitObject->material;
            //auto scattered = shaderPrograms[mtlHandle.index()]->shade(r, hitObject->hitPoint, hitObject->normal);
            //auto scatteredRay = scattered.ray;
            //auto attenuation = scattered.attenuation;
            //auto emitted = scattered.emitted;

            // 把相交点信息添加进path
            Vec3d p = hitObject->hitPoint, n = hitObject->normal;
            Vec3d dir = { r.direction.x, r.direction.y, r.direction.z };
            n = glm::dot(n, dir) < 0 ? n : n * -1.0; // 让n与r方向相反
            path.x[path.n] = Vert(p, n, hitObject->id);
            path.n++;
            const double rnd0 = prnds[(currDepth - 1) * NumRNGsPerEvent + 0 + PathRndsOffset];
            const double rnd1 = prnds[(currDepth - 1) * NumRNGsPerEvent + 1 + PathRndsOffset];

            // 对从交点发出的散射光线进行渲染
            Ray nr;
            nr.origin = p;
            nr.direction = VecCosine(n, 1.0, rnd0, rnd1); // 漫反射材质的反射光方向

            trace(path, nr, currDepth + 1);

            //auto next = trace(path, nr, currDepth + 1);
            //double n_dot_in = glm::dot(hitObject->normal, scatteredRay.direction);
            //double pdf = scattered.pdf;

        }
        else if (t != DOUBLE_INF) {
            // cout << "int trace(), hit the light source." << endl;
            path.x[path.n] = Vert(hitLight->hitPoint, hitLight->normal, light_id);
            path.n++;
        }
    }

    // 沿单位球面分布的向量
    Vec3d MetropolisRenderer::VecRandom(const double rnd1, const double rnd2)
    {
        const double temp1 = 2.0 * PI * rnd1, temp2 = 2.0 * rnd2 - 1.0;				// temp1 in [0, 2pi), temp2 in [-1, 1)
        const double s = sin(temp1), c = cos(temp1), t = sqrt(1.0 - temp2 * temp2); // s, c确定两个方位角，t确定维度
        return Vec3d(s * t, temp2, c * t);
    }

    Vec3d MetropolisRenderer::VecCosine(const Vec3d n, const double g, const double rnd1, const double rnd2)
    {
        // g->inf, temp2->1-, t->0+, 余弦分布就越均匀；反之g = 1，则余弦分布越容易得到法线方向
        const double temp1 = 2.0 * PI * rnd1, temp2 = pow(rnd2, 1.0 / (g + 1.0));
        const double s = sin(temp1), c = cos(temp1), t = sqrt(1.0 - temp2 * temp2);
        Vec3d base = Vec3d(s * t, temp2, c * t);
        return onb(base, n);
    }



}