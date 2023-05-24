#ifndef __BVHTREE_HPP__
#define __BVHREE_HPP__

#include <optional>
#include <vector>
#include <list>
#include <algorithm>
#include "scene/Scene.hpp"
#include "Bounds3.hpp"
#include "Ray.hpp"
#include "intersections/HitRecord.hpp"
#include "intersections/intersections.hpp"

using namespace std;

namespace AccPathTracer {
	class BVHNode {
	public:
		Bounds3 bound;
		vector<Bounds3> bounds;//非叶子节点后面所有的盒
		BVHNode* left;
		BVHNode* right;
		BVHNode() {
			bound = Bounds3();
			left = nullptr;
			right = nullptr;
		}
		BVHNode(Bounds3 b, BVHNode* l, BVHNode* r) {
			bound = b;
			left = l;
			right = r;
		}
		void buildBounds(SharedScene spscene, vector<Node> objects) {
			for (int i = 0; i < objects.size(); i++) {
				if (objects[i].type == Node::Type::SPHERE) {
					bounds.push_back(Bounds3(&spscene->sphereBuffer[objects[i].entity]));
				}
				if (objects[i].type == Node::Type::TRIANGLE) {
					bounds.push_back(Bounds3(&spscene->triangleBuffer[objects[i].entity]));
				}
			}
		}
	};

	class BVHTree {
	public:
		BVHNode* root;
		SharedScene spscene;
		inline BVHNode* build(vector<Bounds3> bound);
		vector<Bounds3> bounds;
		BVHTree() {};
		BVHTree(SharedScene spScene) {
			spscene = spScene;
			root = nullptr;
		}
		inline HitRecord Intersect(const Ray& ray);
		inline HitRecord getIntersect(const Ray& ray, BVHNode* root);
	};

	inline Bounds3 Union(const Bounds3& b1, const Bounds3& b2) {
		Bounds3 ret;
		ret.min = Min(b1.min, b2.min);
		ret.max = Max(b1.max, b2.max);
		return ret;
	}

	inline Bounds3 Union(const Bounds3& b, const Vec3& p)
	{
		Bounds3 ret;
		ret.min = Min(b.min, p);
		ret.max = Max(b.max, p);
		return ret;
	}

	inline HitRecord BVHTree::getIntersect(const Ray& ray,BVHNode* node) {
		HitRecord inter;
		std::array<int, 3> dirIsNeg;
		dirIsNeg[0] = int(ray.direction[0] > 0);
		dirIsNeg[1] = int(ray.direction[1] > 0);
		dirIsNeg[2] = int(ray.direction[2] > 0);
		auto direction_inv = Vec3{ 1. / ray.direction.x, 1. / ray.direction.y, 1. / ray.direction.z };
		if (!node->bound.IntersectP(ray, direction_inv, dirIsNeg)) {
			inter = getMissRecord();
			return inter;
		}
		if (node->left == nullptr && node->right == nullptr) {
			//叶子节点，求改节点对应实体与光线的交点
			if (node->bound.type == Bounds3::Type::SPHERE) {
				inter = Intersection::xSphere(ray, *node->bound.sp, 0.00001, FLOAT_INF);
			}
			if (node->bound.type == Bounds3::Type::TRIANGLE) {
				inter = Intersection::xTriangle(ray, *node->bound.tr, 0.00001, FLOAT_INF);
			}
			return inter;
		}
		HitRecord hit1 = this->getIntersect(ray,node->left);
		HitRecord hit2 = this->getIntersect(ray,node->right);
		if (hit1 == nullopt && hit2 == nullopt) {
			return nullopt;
		}
		else if (hit1!= nullopt) {
			return hit1;
		}
		else if (hit2!= nullopt) {
			return hit2;
		}
		else {
			hit1->t < hit2->t ? hit1 : hit2;
		}
	}

	inline HitRecord BVHTree::Intersect(const Ray& ray) {
		if (!root) {
			return getMissRecord();
		}
		root->bounds = bounds;
		HitRecord hit = BVHTree::getIntersect(ray, root);
		return hit;
	}

	inline BVHNode* BVHTree::build(vector<Bounds3> bounds) {
		BVHNode* node = new BVHNode;
		node->bounds = bounds;
		if (node->bounds.size() == 1) {
			node->bound = node->bounds[0];
			node->left = nullptr;
			node->right = nullptr;
			//只有在叶子节点才为object赋值
			return node;
		}
		else if (node->bounds.size() == 2) {
			node->left = build(vector{ node->bounds[0] });
			node->right = build(vector{ node->bounds[1] });
			node->bound = Union(node->left->bound, node->right->bound);
			return node;
		}
		else {
			Bounds3 centroidBounds;
			for (int i = 0; i < node->bounds.size(); ++i) {
				centroidBounds = Union(centroidBounds, node->bounds[i].Centroid());
			}
			int dim = centroidBounds.maxExtent();
			switch (dim) {
			case 0:
				std::sort(node->bounds.begin(), node->bounds.end(), [](auto f1, auto f2) {
					return f1.Centroid().x < f2.Centroid().x;
					});
				break;
			case 1:
				std::sort(node->bounds.begin(), node->bounds.end(), [](auto f1, auto f2) {
					return f1.Centroid().y < f2.Centroid().y;
					});
				break;
			case 2:
				std::sort(node->bounds.begin(), node->bounds.end(), [](auto f1, auto f2) {
					return f1.Centroid().z < f2.Centroid().z;
					});
				break;
			}
			auto beginning = node->bounds.begin();
			auto middling = node->bounds.begin() + (node->bounds.size() / 2);
			auto ending = node->bounds.end();
			auto leftshapes = vector<Bounds3>(beginning, middling);
			auto rightshapes = vector<Bounds3>(middling, ending);

			assert(node->bounds.size() == (leftshapes.size() + rightshapes.size()));

			node->left = build(leftshapes);
			node->right = build(rightshapes);

			node->bound = Union(node->left->bound, node->right->bound);
		}
		return node;
	}
}
#endif