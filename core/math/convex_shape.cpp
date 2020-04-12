/*************************************************************************/
/*  collision_shape_3d.cpp                                               */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2020 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2020 Godot Engine contributors (cf. AUTHORS.md).   */
/*                                                                       */
/* Permission is hereby granted, free of charge, to any person obtaining */
/* a copy of this software and associated documentation files (the       */
/* "Software"), to deal in the Software without restriction, including   */
/* without limitation the rights to use, copy, modify, merge, publish,   */
/* distribute, sublicense, and/or sell copies of the Software, and to    */
/* permit persons to whom the Software is furnished to do so, subject to */
/* the following conditions:                                             */
/*                                                                       */
/* The above copyright notice and this permission notice shall be        */
/* included in all copies or substantial portions of the Software.       */
/*                                                                       */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*/
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE     */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                */
/*************************************************************************/

#include "convex_shape.h"

#include <string.h>

ConvexShape::ConvexShape(const Plane *p_planes, int p_plane_count) {
	set_planes(p_planes, p_plane_count);
}

ConvexShape::ConvexShape(const Plane *p_planes, int p_plane_count, const Vector3 *p_points, int p_point_count) {
	set_planes_and_points(p_planes, p_plane_count, p_points, p_point_count);
}

ConvexShape::ConvexShape() {
}

void ConvexShape::set_planes(const Plane *p_planes, int p_plane_count) {
	planes.resize(p_plane_count);
	memcpy(planes.ptrw(), p_planes, p_plane_count * sizeof(Plane));

	_compute_points_from_planes();
}

void ConvexShape::set_planes_and_points(const Plane *p_planes, int p_plane_count, const Vector3 *p_points, int p_point_count) {
	planes.resize(p_plane_count);
	memcpy(planes.ptrw(), p_planes, p_plane_count * sizeof(Plane));

	points.resize(p_point_count);
	memcpy(points.ptrw(), p_points, p_point_count * sizeof(Vector3));
}

void ConvexShape::_compute_points_from_planes() {

	WARN_PRINT("Doing slow point calculation for ConvexShape.");

	points.clear();

	// Do initial intersection tests.
	for (int i = planes.size() - 1; i >= 0; i--) {
		for (int j = i - 1; j >= 0; j--) {
			for (int k = j - 1; k >= 0; k--) {

				//std::cout << "  " << plane_names[i] << ", " << plane_names[j] << ", " << plane_names[k] << std::endl;
				Vector3 convex_shape_point;
				if (planes[i].intersect_3(planes[j], planes[k], &convex_shape_point)) {

					bool excluded = false;

					/*std::cout << "    GOOD: "
						<< convex_shape_point.x << ", "
						<< convex_shape_point.y << ", "
						<< convex_shape_point.z << std::endl;*/

					// See if any other plane excludes this point.
					for (int n = 0; n < planes.size(); n++) {
						if (n != i && n != j && n != k) {
							real_t dp = planes[n].normal.dot(convex_shape_point);
							if (dp - planes[n].d > CMP_EPSILON) {
								/*std::cout << "      BUT... Excluded by " << plane_names[n] << std::endl; */
								excluded = true;
								break;
							}
						}
					}

					if (!excluded) {
						points.push_back(convex_shape_point);
					}

				} else {
					//std::cout << "    BAD" << std::endl;
				}
			}
		}
	}
}
