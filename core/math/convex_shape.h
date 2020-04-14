/*************************************************************************/
/*  convex_shape.h                                                       */
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

#ifndef CONVEX_SHAPE_H
#define CONVEX_SHAPE_H

#include "core/math/plane.h"
#include "core/math/vector3.h"

class ConvexShape {
public:
	ConvexShape();
	ConvexShape(const Plane *p_planes, int p_plane_count);
	ConvexShape(const Plane *p_planes, int p_plane_count, const Vector3 *p_points, int p_point_count);

	// Set the shape based on an array of planes. Points are calculated based on planes.
	void set_planes(const Plane *p_planes, int p_plane_count);

	// Set the shape, and use precomputed points. In case the creator already knows where the corners are.
	void set_planes_and_points(const Plane *p_planes, int p_plane_count, const Vector3 *p_points, int p_point_count);

	Vector<Vector3> points;
	Vector<Plane> planes;

private:
	void _compute_points_from_planes();
};

#endif // CONVEX_SHAPE_H
