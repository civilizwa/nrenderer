#include "Metropolis.hpp"
#include "server/Server.hpp"
#include "VertexTransformer.hpp"
#include "intersections/intersections.hpp"
#include "glm/gtc/matrix_transform.hpp"

const auto taskNums = 4;
std::mutex mtx;

namespace Metropolis {
    RGB MetropolisRenderer::gamma(const RGB& rgb) {
        return glm::sqrt(rgb);
    }
    
    RGBA MetropolisRenderer::gamma(const RGBA& rgba) {
        return glm::sqrt(rgba);
    }

    RGBA MetropolisRenderer::gamma(const RGBA& rgba, unsigned long samples) {
        return Vec4{ pow(1 - exp(-rgba[0] * samples), 1 / 2.2), pow(1 - exp(-rgba[1] * samples), 1 / 2.2) ,
            pow(1 - exp(-rgba[2] * samples), 1 / 2.2), 0};
    }


    void MetropolisRenderer::renderTask(RGBA* pixels, int width, int height, int off, int step) {
        // initialize the Markov chain
        TMarkovChain current, proposal;
        InitRandomNumbersByChain(current); // 将current的已经生成好的随机数数组u[]赋值给prnd[]
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

        //renderTask(pixels, width, height, 0, 1);

        // estimate normalization constant b是用双向路径追踪计算照片最后的平均亮度
        for (int i = 0; i < N_Init; i++) {
            InitRandomNumbers();
            double sc = CombinePaths(GenerateEyePath(MaxEvents), GenerateLightPath(MaxEvents)).sc;
            b += sc;
            std::lock_guard<std::mutex> lock(mtx); // 加锁

            //cout << "i = " << i << ", sc = " << sc << endl;

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
        cout << "b = " << b << endl; // b = 0.136614

        thread t[taskNums];
        for (int i = 0; i < taskNums; i++) {
            t[i] = thread(&MetropolisRenderer::renderTask,
                this, pixels, width, height, i, taskNums);
        }
        for (int i = 0; i < taskNums; i++) {
            t[i].join();
        }




        getServer().logger.log("Done...");

        // 设置透明度以及gamma校正
        cout << "samps = " << samps << endl;
        double s = double(width * height) / double(samps);

        for (int ix = 0; ix < spScene->renderOption.width; ix++) {
            for (int iy = 0; iy < spScene->renderOption.height; iy++) {
                // 框架中的校正
                // pixels[ix + iy * width] = gamma(pixels[ix + iy * width]);

                // MLT中的校正
                double gamma = 2.2;
                pixels[ix + iy * width][0] = pow(1 - exp(-pixels[ix + iy * width][0] * s), 1 / gamma);
                pixels[ix + iy * width][1] = pow(1 - exp(-pixels[ix + iy * width][1] * s), 1 / gamma);
                pixels[ix + iy * width][2] = pow(1 - exp(-pixels[ix + iy * width][2] * s), 1 / gamma);
                pixels[ix + iy * width][3] = 1.0;
            }
        }
            

        return { pixels, width, height };
    }

    void MetropolisRenderer::release(const RenderResult& r) {
        auto [p, w, h] = r;
        delete[] p;
    }

    HitRecord MetropolisRenderer::closestHitObject(const Ray& r) {
        //HitRecord closestHit = nullopt;
        //float closest = FLOAT_INF;

        // BVH
        //auto hitRecord = tree->Intersect(r, closest);
        //if (hitRecord && hitRecord->t < closest) {
        //    closest = hitRecord->t;
        //    closestHit = hitRecord;
        //}

        //return closestHit;

        int id = 0;

        HitRecord closestHit = nullopt;
        float closest = FLOAT_INF;
        for (auto& p : scene.planeBuffer) {
            auto hitRecord = Intersection::xPlane(r, p, 0.000001, closest, id);
            id++;
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                // cout << "in PLANE: " << closest << endl;
                closestHit = hitRecord;
            }
        }
        for (auto& s : scene.sphereBuffer) {
            auto hitRecord = Intersection::xSphere(r, s, 0.000001, closest, id);
            id++;
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                // cout << "in SPHERE: " << closest << endl;
                closestHit = hitRecord;
            }
        }
        for (auto& t : scene.triangleBuffer) {
            auto hitRecord = Intersection::xTriangle(r, t, 0.000001, closest, id);
            id++;
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                // cout << "in TRIANGLE: " << closest << endl;
                closestHit = hitRecord;
            }
        }

        return closestHit;
    }

    tuple<double, Vec3d, HitRecord> MetropolisRenderer::closestHitLight(const Ray& r) {
        Vec3d v = {};
        HitRecord closest = getHitRecord(FLOAT_INF, {}, {}, {}, {});
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
            Vec3d p = hitObject->hitPoint, n = glm::normalize(hitObject->normal);
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
                

            // 这里和源程序相比比较奇怪的一点是，trace只是往path中添加了交点信息，而对于相交的颜色的结果是选择丢弃的，到了PathThroughput函数中才去计算某一条Path的color contribution
            // 这里是不是一个可以改进的地方呢？

            // auto next = trace(path, scatteredRay, currDepth + 1);
            //float n_dot_in = glm::dot(hitObject->normal, scatteredRay.direction);//cos(N法向量,L光源)
            //float pdf = scattered.pdf;
            ///**
            //    * emitted      - Le(p, w_0)
            //    * next         - Li(p, w_i)
            //    * n_dot_in     - cos<n, w_i>
            //    * atteunation  - BRDF
            //    * pdf          - p(w)
            //    **/
            //    //emitted:直接光照发出的强度
            //    //后部分是漫反射
            //return emitted + attenuation * next * n_dot_in / pdf;

        }
        else if (t != FLOAT_INF) {
            path.x[path.n] = Vert(hitLight->hitPoint, hitLight->normal, light_id);
            path.n++;
        }
    }

    // 沿单位球面分布的向量
    //Vec3d MetropolisRenderer::VecRandom(const double rnd1, const double rnd2)
    //{
    //    const double temp1 = 2.0 * PI * rnd1, temp2 = 2.0 * rnd2 - 1.0;				// temp1 in [0, 2pi), temp2 in [-1, 1)
    //    const double s = sin(temp1), c = cos(temp1), t = sqrt(1.0 - temp2 * temp2); // s, c确定两个方位角，t确定维度
    //    return Vec3d(s * t, temp2, c * t);
    //}

    // 返回沿下半球面一定角度范围内分布的单位向量
    Vec3d MetropolisRenderer::VecRandom(const double rnd1, const double rnd2)
    {
        const double temp1 = 2.0 * PI * rnd1; // [0, 2pi]
        const double temp2 = -0.6 + rnd2 * -0.2;

        const double s = sin(temp1);
        const double c = cos(temp1);
        const double t = sqrt(1.0 - temp2 * temp2);

        return Vec3d(s * t, temp2, c * t);
    }

    // 返回余弦分布的向量
    Vec3d MetropolisRenderer::VecCosine(const Vec3d n, const double g, const double rnd1, const double rnd2)
    {
        // g趋向于无穷，则函数返回结果就是n；g趋向于-1，函数返回（应该？）在球面上均匀分布的向量
        const double temp1 = 2.0 * PI * rnd1, temp2 = pow(rnd2, 1.0 / (g + 1.0));
        const double s = sin(temp1), c = cos(temp1), t = sqrt(1.0 - temp2 * temp2);
        Vec3d base = Vec3d(s * t, temp2, c * t);
        return onb(base, n);
    }
}