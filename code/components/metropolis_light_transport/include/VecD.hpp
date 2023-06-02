#pragma once
#ifndef __VEC_D__
#define __VEC_D__

#include "Metropolis.hpp"

namespace Metropolis
{
	struct VecD
	{
		double x, y, z;
		VecD(double x_ = 0, double y_ = 0, double z_ = 0)
		{
			x = x_;
			y = y_;
			z = z_;
		}
		inline VecD operator-() const { return VecD(-x, -y, -z); }
		inline VecD operator+(const VecD& b) const { return VecD(x + b.x, y + b.y, z + b.z); }
		inline VecD operator-(const VecD& b) const { return VecD(x - b.x, y - b.y, z - b.z); }
		inline VecD operator+(double b) const { return VecD(x + b, y + b, z + b); }
		inline VecD operator-(double b) const { return VecD(x - b, y - b, z - b); }
		inline VecD operator*(double b) const { return VecD(x * b, y * b, z * b); }
		inline VecD mul(const VecD& b) const { return VecD(x * b.x, y * b.y, z * b.z); }							   // 对应元素乘
		inline VecD norm() { return (*this) * (1.0 / sqrt(x * x + y * y + z * z)); }								   // 归一化
		inline double dot(const VecD& b) const { return x * b.x + y * b.y + z * b.z; }							   // 点乘
		VecD operator%(const VecD& b) const { return VecD(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); } // 叉乘
		inline VecD reflect(const VecD& n) const { return (*this) - n * 2.0 * n.dot((*this)); }					   // 给定法向量n计算反射向量
		inline double Max() const { return MAX(MAX(x, y), z); }
	};

	std::ostream& operator<<(std::ostream& out, const VecD& v)  // 友元函数定义
	{
		return out << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
	}
}


#endif