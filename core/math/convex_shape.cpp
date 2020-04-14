/*************************************************************************/
/*  convex_shape.cpp                                                     */
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
	for (int i = 0; i < p_plane_count; i++) {
		planes.write[i] = p_planes[i];
	}

	_compute_points_from_planes();
}

void ConvexShape::set_planes_and_points(const Plane *p_planes, int p_plane_count, const Vector3 *p_points, int p_point_count) {
	planes.resize(p_plane_count);
	for (int i = 0; i < p_plane_count; i++) {
		planes.write[i] = p_planes[i];
	}

	points.resize(p_point_count);
	for (int i = 0; i < p_point_count; i++) {
		points.write[i] = p_points[i];
	}
}

void ConvexShape::_compute_points_from_planes() {

	points.clear();

	// Iterate through every unique combination of any three planes.
	for (int i = planes.size() - 1; i >= 0; i--) {
		for (int j = i - 1; j >= 0; j--) {
			for (int k = j - 1; k >= 0; k--) {

				// Find the point where these planes all cross over (if they
				// do at all).
				Vector3 convex_shape_point;
				if (planes[i].intersect_3(planes[j], planes[k], &convex_shape_point)) {

					// See if any *other* plane excludes this point because it's
					// on the wrong side.
					bool excluded = false;
					for (int n = 0; n < planes.size(); n++) {
						if (n != i && n != j && n != k) {
							real_t dp = planes[n].normal.dot(convex_shape_point);
							if (dp - planes[n].d > CMP_EPSILON) {
								excluded = true;
								break;
							}
						}
					}

					// Only add the point if it passed all tests.
					if (!excluded) {
						points.push_back(convex_shape_point);
					}
				}
			}
		}
	}
}
