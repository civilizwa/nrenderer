#pragma once
#ifndef __PATH_CONTRIBUTION__
#define __PATH_CONTRIBUTION__

#include "Metropolis.hpp"

#define PI ((double)3.1415926535898)
#define MAX(x, y) ((x > y) ? x : y)
#define MIN(x, y) ((x < y) ? x : y)

// parameters
#define MinPathLength 3 // avoid sampling direct illumination. light subpath + eye subpath + 连接，一条完整的BPT path至少3段
#define MaxPathLength 20
#define Glossiness 25.0
#define N_Init 10000
#define LargeStepProb 0.3

// scene independent constants
#define NumRNGsPerEvent 2									 // 每个事件使用的随机数数量
#define MaxEvents (MaxPathLength + 1)						 // 路径中允许的最大事件（发射、散射、接收光线）数量
#define NumStatesSubpath ((MaxEvents + 2) * NumRNGsPerEvent) // 子路径的状态信息（例如光的颜色、方向）数量
#define NumStates (NumStatesSubpath * 2)					 // 完整路径的状态数量，一个完整路径有两个子路径

#define light_id -3
#define light_area (4.0 * PI * sph[light_id].rad * sph[light_id].rad) // 光源表面积
#define camera_id -2


namespace Metropolis
{
	// path data
	struct Vert
	{
		Vec3 p, n;
		int id;
		Vert() {};
		Vert(Vec3 p_, Vec3 n_, int id_) : p(p_), n(n_), id(id_) {}
	};
	struct Path
	{
		Vert x[MaxEvents];
		int n; // 一条路径中的vertex数量
		Path() { n = 0; }
	};
	struct Contribution
	{
		double x, y;
		Vec3 c; // 颜色
		Contribution() {};
		Contribution(double x_, double y_, Vec3 c_) : x(x_), y(y_), c(c_) {}
	};
	struct PathContribution
	{
		Contribution c[MaxEvents * MaxEvents];
		int n;
		double sc; // scalar contribution function
		PathContribution()
		{
			n = 0;
			sc = 0.0;
		}
	};
}

#endif