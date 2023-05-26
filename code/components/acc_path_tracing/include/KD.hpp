#ifndef __KDTREE_HPP__
#define __KDTREE_HPP__

#define MIN_BOX_SIZE 10
#define MAX_OBJ_BOX 100

#include <vector>
#include <list>
#include <algorithm>

#include "scene/Scene.hpp"
#include "Ray.hpp"
#include "Bounds3.hpp"
using namespace NRenderer;
namespace AccPathTracer {
    struct KDNode {
        KDNode() {
            double minNum = std::numeric_limits<double>::lowest();
            double maxNum = std::numeric_limits<double>::max();
            min = Vec3{ minNum ,minNum ,minNum };
            max = Vec3{ maxNum ,maxNum ,maxNum };
            left = right = nullptr;
        }
        KDNode(Bounds3* box) {
            double minNum = std::numeric_limits<double>::lowest();
            double maxNum = std::numeric_limits<double>::max();
            min = Vec3{ minNum ,minNum ,minNum };
            max = Vec3{ maxNum ,maxNum ,maxNum };
            Update(box);
            left = right = nullptr;
            Box_list.push_back(box);
        }

        /*KDNode(KDNode* n, int index, bool flag, int v) {
            double minNum = std::numeric_limits<double>::lowest();
            double maxNum = std::numeric_limits<double>::max();
            min = n->min;
            max = n->max;
            if (flag) {
                max[index] = v;
            }
            else {
                min[index] = v;
            }
            left = right = nullptr;
        }*/

        void Update(Bounds3* box) {
            min.x = std::min(min.x, box->min.x);
            min.y = std::min(min.y, box->min.y);
            min.z = std::min(min.z, box->min.z);
            max.x = std::max(max.x, box->max.x);
            max.y = std::max(max.y, box->max.y);
            max.z = std::max(max.z, box->max.z);

        }

        bool IsInNode(const Bounds3* box) {
            if (box->min.x >= min.x && box->min.y >= -min.y && box->min.z >= min.z 
                && box->max.x <= max.x && box->max.y < max.y && box->max.z <= max.z) {
                return true;
            }
            else {
                return false;
            }
        }

        void Insert(Bounds3* box) {
            Update(box);
            Box_list.push_back(box);
        }

        Vec3 min;
        Vec3 max;
        std::list<Bounds3*> Box_list;
        KDNode* left, * right;
        bool is_leaf = true;
        float r;
    };

    struct KDTree {
        KDTree() { root = new KDNode; }
        void Insert(std::vector<Bounds3*>& boxs, int s, int e, KDNode* n) {
            for (int i = s; i < e; i++) {
                n->Update(boxs[i]);
            }
            if (s == e) return;
            else if (e - s <= MAX_SIZE) {
                for (int i = s; i < e; i++) {
                    n->Insert(boxs[i]);
                }
            }
            else {
                for (int i = s; i < e; i++) {
                    n->Update(boxs[i]);
                }
                n->is_leaf = false;
                n->left = new KDNode;
                n->right = new KDNode;
                std::sort(boxs.begin() + s, boxs.begin() + e, [&](const Bounds3* A, const Bounds3* B) {
                    if (index == 1) {
                        return A->min.x < B->min.x;
                    }
                    else if (index == 2) {
                        return A->min.y < B->min.y;
                    }
                    else {
                        return A->min.z < B->min.z;
                    }
                    });
                index = (index + 1) % 3;
                int mid = (s + e) / 2;
                Insert(boxs, s, mid, n->left);
                Insert(boxs, mid, e, n->right);
            }
        }
        //void Insert(std::vector<NRenderer::Sphere>& sps, std::vector<NRenderer::Triangle>& trs) {
        //    std::vector<Bounds3*> boxs(sps.size() + trs.size());
        //    int index = 0;
        //    for (auto& sp : sps) {
        //        Mat4x4 t{ 1 };
        //        boxs[index++] = new Bounds3(&sp);
        //    }
        //    for (auto& tr : trs) {
        //        boxs[index++] = new Bounds3(&tr);
        //    }
        //    Insert(boxs, 0, boxs.size(), root);
        //}

        bool IsHit(const AccPathTracer::Ray& r, KDNode* n) {
            auto inv_dir = 1.0f / r.direction;
            Vec3 cube_min(n->min[0], n->min[1], n->min[2]), cube_max(n->max[0], n->max[1], n->max[2]);

            auto tMin = (cube_min - r.origin) * inv_dir;
            auto tMax = (cube_max - r.origin) * inv_dir;
            auto t1 = min(tMin, tMax);
            auto t2 = max(tMin, tMax);
            float tNear = std::max(std::max(t1.x, t1.y), t1.z);
            float tFar = std::min(std::min(t2.x, t2.y), t2.z);

            return  tNear < tFar;
        }

        std::list<KDNode*> find_Node(const AccPathTracer::Ray& r) {
            std::list<KDNode*> result;
            find_Node(r, root, result);
            return result;
        }

        void find_Node(const AccPathTracer::Ray& r, KDNode* n, std::list<KDNode*>& result) {
            if (IsHit(r, n)) {
                if (n->is_leaf) {
                    result.push_back(n);
                }
                else {
                    find_Node(r, n->left, result);
                    find_Node(r, n->right, result);
                }
            }
        }

        KDNode* root = nullptr;
        int MAX_SIZE = 3;
        int index = 0;
    };
}
#endif