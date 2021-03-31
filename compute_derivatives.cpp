/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2010-2011, Willow Garage, Inc.
 *  Copyright (c) 2012-, Open Perception, Inc.
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the copyright holder(s) nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * $Id$
 *
 */
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <vector>

#include <radius_search.h>

/** \brief Precompute angular components of derivatives.
 * \note Equation 6.19 and 6.21 [Magnusson 2009].
 * \param[in] p the current transform vector
 * \param[in] j_ang_a precomputed angular gradient
 * \param[in] j_ang_b precomputed angular gradient
 * \param[in] j_ang_c precomputed angular gradient
 * \param[in] j_ang_d precomputed angular gradient
 * \param[in] j_ang_e precomputed angular gradient
 * \param[in] j_ang_f precomputed angular gradient
 * \param[in] j_ang_g precomputed angular gradient
 * \param[in] j_ang_h precomputed angular gradient
 * \param[in] j_ang precomputed angular gradient
 * \param[in] h_ang_a2 precomputed angular hessian
 * \param[in] h_ang_a3 precomputed angular hessian
 * \param[in] h_ang_b2 precomputed angular hessian
 * \param[in] h_ang_b3 precomputed angular hessian
 * \param[in] h_ang_c2 precomputed angular hessian
 * \param[in] h_ang_c3 precomputed angular hessian
 * \param[in] h_ang_d1 precomputed angular hessian
 * \param[in] h_ang_d2 precomputed angular hessian
 * \param[in] h_ang_d3 precomputed angular hessian
 * \param[in] h_ang_e1 precomputed angular hessian
 * \param[in] h_ang_e2 precomputed angular hessian
 * \param[in] h_ang_e3 precomputed angular hessian
 * \param[in] h_ang_f1 precomputed angular hessian
 * \param[in] h_ang_f2 precomputed angular hessian
 * \param[in] h_ang_f3 precomputed angular hessian
 * \param[in] h_ang precomputed angular hessian
 */
void computeAngleDerivatives(
  const float p[6], float j_ang_a[3], float j_ang_b[3], float j_ang_c[3], float j_ang_d[3],
  float j_ang_e[3], float j_ang_f[3], float j_ang_g[3], float j_ang_h[3], float j_ang[8][4],
  float h_ang_a2[3], float h_ang_a3[3], float h_ang_b2[3], float h_ang_b3[3], float h_ang_c2[3],
  float h_ang_c3[3], float h_ang_d1[3], float h_ang_d2[3], float h_ang_d3[3], float h_ang_e1[3],
  float h_ang_e2[3], float h_ang_e3[3], float h_ang_f1[3], float h_ang_f2[3], float h_ang_f3[3],
  float h_ang[16][4])
{
  // Simplified math for near 0 angles
  float cx, cy, cz, sx, sy, sz;
  if (fabs(p[3]) < 10e-5) {
    // p(3) = 0;
    cx = 1.0;
    sx = 0.0;
  } else {
    cx = cos(p[3]);
    sx = sin(p[3]);
  }
  if (fabs(p[4]) < 10e-5) {
    // p(4) = 0;
    cy = 1.0;
    sy = 0.0;
  } else {
    cy = cos(p[4]);
    sy = sin(p[4]);
  }

  if (fabs(p[5]) < 10e-5) {
    // p(5) = 0;
    cz = 1.0;
    sz = 0.0;
  } else {
    cz = cos(p[5]);
    sz = sin(p[5]);
  }

  // Precomputed angular gradient components. Letters correspond to Equation 6.19 [Magnusson 2009]
  j_ang_a[0] = (-sx * sz + cx * sy * cz);
  j_ang_a[1] = (-sx * cz - cx * sy * sz);
  j_ang_a[2] = (-cx * cy);
  j_ang_b[0] = (cx * sz + sx * sy * cz);
  j_ang_b[1] = (cx * cz - sx * sy * sz);
  j_ang_b[2] = (-sx * cy);
  j_ang_c[0] = (-sy * cz);
  j_ang_c[1] = sy * sz;
  j_ang_c[0] = cy;
  j_ang_d[0] = sx * cy * cz;
  j_ang_d[1] = (-sx * cy * sz);
  j_ang_d[2] = sx * sy;
  j_ang_e[0] = (-cx * cy * cz);
  j_ang_e[1] = cx * cy * sz;
  j_ang_e[2] = (-cx * sy);
  j_ang_f[0] = (-cy * sz);
  j_ang_f[1] = (-cy * cz);
  j_ang_f[2] = 0;
  j_ang_g[0] = (cx * cz - sx * sy * sz);
  j_ang_g[1] = (-cx * sz - sx * sy * cz);
  j_ang_g[2] = 0;
  j_ang_h[0] = (sx * cz + cx * sy * sz);
  j_ang_h[1] = (cx * sy * cz - sx * sz);
  j_ang_h[2] = 0;

  j_ang[0][0] = (-sx * sz + cx * sy * cz);
  j_ang[0][1] = (-sx * cz - cx * sy * sz);
  j_ang[0][2] = (-cx * cy);
  j_ang[0][3] = 0.0f;
  j_ang[1][0] = (cx * sz + sx * sy * cz);
  j_ang[1][1] = (cx * cz - sx * sy * sz);
  j_ang[1][2] = (-sx * cy);
  j_ang[1][3] = 0.0f;
  j_ang[2][0] = (-sy * cz);
  j_ang[2][1] = sy * sz;
  j_ang[2][2] = cy;
  j_ang[2][3] = 0.0f;
  j_ang[3][0] = sx * cy * cz;
  j_ang[3][1] = (-sx * cy * sz);
  j_ang[3][2] = sx * sy, j_ang[3][3] = 0.0f;

  j_ang[4][0] = (-cx * cy * cz);
  j_ang[4][1] = cx * cy * sz;
  j_ang[4][2] = (-cx * sy);
  j_ang[4][3] = 0.0f;

  j_ang[5][0] = (-cy * sz);
  j_ang[5][1] = (-cy * cz);
  j_ang[5][2] = 0;
  j_ang[5][3] = 0.0f;
  j_ang[6][0] = (cx * cz - sx * sy * sz);
  j_ang[6][1] = (-cx * sz - sx * sy * cz);
  j_ang[6][2] = 0;
  j_ang[6][3] = 0.0f;
  j_ang[7][0] = (sx * cz + cx * sy * sz);
  j_ang[7][1] = (cx * sy * cz - sx * sz);
  j_ang[7][2] = 0;
  j_ang[7][3] = 0.0f;

  // Precomputed angular hessian components. Letters correspond to Equation 6.21 and numbers correspond to row index
  // [Magnusson 2009]
  h_ang_a2[0] = (-cx * sz - sx * sy * cz);
  h_ang_a2[1] = (-cx * cz + sx * sy * sz);
  h_ang_a2[2] = sx * cy;
  h_ang_a3[0] = (-sx * sz + cx * sy * cz);
  h_ang_a3[1] = (-cx * sy * sz - sx * cz);
  h_ang_a3[2] = (-cx * cy);

  h_ang_b2[0] = (cx * cy * cz);
  h_ang_b2[1] = (-cx * cy * sz);
  h_ang_b2[2] = (cx * sy);
  h_ang_b3[0] = (sx * cy * cz);
  h_ang_b3[1] = (-sx * cy * sz);
  h_ang_b3[2] = (sx * sy);

  h_ang_c2[0] = (-sx * cz - cx * sy * sz);
  h_ang_c2[1] = (sx * sz - cx * sy * cz);
  h_ang_c2[2] = 0;
  h_ang_c3[0] = (cx * cz - sx * sy * sz);
  h_ang_c3[1] = (-sx * sy * cz - cx * sz);
  h_ang_c3[2] = 0;

  h_ang_d1[0] = (-cy * cz);
  h_ang_d1[1] = (cy * sz);
  h_ang_d1[2] = (sy);
  h_ang_d2[0] = (-sx * sy * cz);
  h_ang_d2[1] = (sx * sy * sz);
  h_ang_d2[2] = (sx * cy);
  h_ang_d3[0] = (cx * sy * cz);
  h_ang_d3[1] = (-cx * sy * sz);
  h_ang_d3[2] = (-cx * cy);

  h_ang_e1[0] = (sy * sz);
  h_ang_e1[1] = (sy * cz);
  h_ang_e1[2] = 0;
  h_ang_e2[0] = (-sx * cy * sz);
  h_ang_e1[1] = (-sx * cy * cz);
  h_ang_e1[2] = 0;
  h_ang_e3[0] = (cx * cy * sz);
  h_ang_e3[1] = (cx * cy * cz);
  h_ang_e3[2] = 0;

  h_ang_f1[0] = (-cy * cz);
  h_ang_f1[1] = (cy * sz);
  h_ang_f1[2] = 0;
  h_ang_f2[0] = (-cx * sz - sx * sy * cz);
  h_ang_f2[1] = (-cx * cz + sx * sy * sz);
  h_ang_f2[2] = 0;
  h_ang_f3[0] = (-sx * sz + cx * sy * cz);
  h_ang_f3[1] = (-cx * sy * sz - sx * cz);
  h_ang_f3[2] = 0;

  // a2
  h_ang[0][0] = (-cx * sz - sx * sy * cz);
  h_ang[0][1] = (-cx * cz + sx * sy * sz);
  h_ang[0][2] = sx * cy;
  h_ang[0][3] = 0.0f;
  // a3
  h_ang[1][0] = (-sx * sz + cx * sy * cz);
  h_ang[1][1] = (-cx * sy * sz - sx * cz);
  h_ang[1][2] = (-cx * cy);
  h_ang[1][3] = 0.0f;
  // b2
  h_ang[2][0] = (cx * cy * cz);
  h_ang[2][1] = (-cx * cy * sz);
  h_ang[2][2] = (cx * sy);
  h_ang[2][3] = 0.0f;
  // b3
  h_ang[3][0] = (sx * cy * cz);
  h_ang[3][1] = (-sx * cy * sz);
  h_ang[3][2] = (sx * sy);
  h_ang[3][3] = 0.0f;
  // c2
  h_ang[4][0] = (-sx * cz - cx * sy * sz);
  h_ang[4][1] = (sx * sz - cx * sy * cz);
  h_ang[4][2] = 0;
  h_ang[4][3] = 0.0f;
  // c3
  h_ang[5][0] = (cx * cz - sx * sy * sz);
  h_ang[5][1] = (-sx * sy * cz - cx * sz);
  h_ang[5][2] = 0;
  h_ang[5][3] = 0.0f;

  // d1
  h_ang[6][0] = (-cy * cz);
  h_ang[6][1] = (cy * sz);
  h_ang[6][2] = (sy);
  h_ang[6][3] = 0.0f;
  // d2
  h_ang[7][0] = (-sx * sy * cz);
  h_ang[7][1] = (sx * sy * sz);
  h_ang[7][2] = (sx * cy);
  h_ang[7][3] = 0.0f;

  // d3
  h_ang[8][0] = (cx * sy * cz);
  h_ang[8][1] = (-cx * sy * sz);
  h_ang[8][2] = (-cx * cy);
  h_ang[8][3] = 0.0f;
  // e1
  h_ang[9][0] = (sy * sz);
  h_ang[9][1] = (sy * cz);
  h_ang[9][2] = 0;
  h_ang[9][3] = 0.0f;
  // e2
  h_ang[10][0] = (-sx * cy * sz);
  h_ang[10][1] = (-sx * cy * cz);
  h_ang[10][2] = 0;
  h_ang[10][3] = 0.0f;
  // e3
  h_ang[11][0] = (cx * cy * sz);
  h_ang[11][1] = (cx * cy * cz);
  h_ang[11][2] = 0;
  h_ang[11][3] = 0.0f;

  // f1
  h_ang[12][0] = (-cy * cz);
  h_ang[12][1] = (cy * sz);
  h_ang[12][2] = 0;
  h_ang[12][3] = 0.0f;
  // f2
  h_ang[13][0] = (-cx * sz - sx * sy * cz);
  h_ang[13][1] = (-cx * cz + sx * sy * sz);
  h_ang[13][2] = 0;
  h_ang[13][3] = 0.0f;
  // f2
  h_ang[14][0] = (-sx * sz + cx * sy * cz);
  h_ang[14][1] = (-cx * sy * sz - sx * cz);
  h_ang[14][2] = 0;
  h_ang[14][3] = 0.0f;
  // f3
  h_ang[15][0] = 0.0f;
  h_ang[15][1] = 0.0f;
  h_ang[15][2] = 0.0f;
  h_ang[15][3] = 0.0f;
}

/**
 * @brief return the indexes of 7 neighbor voxels of the reference point
 * @param lidar_points_x       (input) 1-D array of x value of input points (n elements)
 * @param lidar_points_y       (input) 1-D array of y value of input points (n elements)
 * @param lidar_points_z       (input) 1-D array of z value of input points (n elements)
 * @param map_points_x         (input) 1-D array of x value of map points (? elements)
 * @param map_points_y         (input) 1-D array of y value of map points (? elements)
 * @param map_points_z         (input) 1-D array of z value of map points (? elements)
 * @param node_indexes         (input) 1-D array of indexes of map points (? elements)
 * @param root_node            (input) kdtree root node
 * @param n                    (input) Number of input points
 * @param limit                (input) limit of neighbors
 * @param radius               (input) search radius
 * @param neighbor_candidate_indexes (output) 1-D array of indexes of neighbor candidates (? elements)
 * @param neighbor_candidate_dists   (output) 1-D array of distances of neighbor candidates (? elements)
 */

void radiusSearchC(
  const float * lidar_points_x, const float * lidar_points_y, const float * lidar_points_z,
  const float * map_points_x, const float * map_points_y, const float * map_points_z,
  const int * node_indexes, const kdtree_node * root_node, const int n, const int limit,
  const float radius, int * neighbor_candidate_indexes, float * neighbor_candidate_dists)
{
  // 1d range kernel for the point cloud
  for (int item_index = 0; item_index < n; item_index++) {
    float reference_point[4] = {
      lidar_points_x[item_index], lidar_points_y[item_index], lidar_points_z[item_index], 0.0f};

    const kdtree_node * current_node = root_node;
    int neighbors_count = 0;

    /*
    if (item_index == 0) {
      printf("item:%d, root_node (at :%p)(%f, %f, %f), child1:%p, child2:%p\n", item_index, current_node, current_node->location.x, current_node->location.y, current_node->location.z, current_node->child1, current_node->child2);
      printf("item:%d, map_points (%f, %f, %f), \n", item_index, map_points_x[item_index], map_points_y[item_index], map_points_z[item_index]);
    }
    */

    // radius search
    float dist_square[4];
    float dist;
    while (true) {
      if ((current_node->child1 == NULL) && (current_node->child2 == NULL)) {
        //printf("point %d (%d to %d)\n", item_index, current_node->left_index, current_node->right_index);
        for (int i = current_node->left_index; i < current_node->right_index; ++i) {
          int index = node_indexes[i];
          float map_point[4] = {
            map_points_x[index], map_points_y[index], map_points_z[index], 0.0f};
          //float dist_square [4];
          for (int d = 0; d < 4; d++) {
            dist_square[d] =
              (map_point[d] - reference_point[d]) * (map_point[d] - reference_point[d]);
          }
          dist = sqrtf(dist_square[0] + dist_square[1] + dist_square[2]);
          if (dist < radius) {
            //printf("leaf node, get radius, %d\n", limit*item_index + neighbors_count);
            neighbor_candidate_indexes[limit * item_index + neighbors_count] = index;
            neighbor_candidate_dists[limit * item_index + neighbors_count] = dist;
            neighbors_count++;
          }
          if (neighbors_count >= limit) {
            break;
          }
        }
        break;
      } else {
        int index = node_indexes
          [current_node->left_index +
           (current_node->right_index - current_node->left_index - 1) / 2];
        float map_point[4] = {map_points_x[index], map_points_y[index], map_points_z[index], 0.0f};
        for (int d = 0; d < 4; d++) {
          dist_square[d] =
            (map_point[d] - reference_point[d]) * (map_point[d] - reference_point[d]);
        }
        float dist = sqrtf(dist_square[0] + dist_square[1] + dist_square[2]);
        if (dist < radius) {
          //printf("non-leaf node, get radius, %d\n", limit*item_index + neighbors_count);
          neighbor_candidate_indexes[limit * item_index + neighbors_count] = index;
          neighbor_candidate_dists[limit * item_index + neighbors_count] = dist;
          neighbors_count++;
        }
        if (neighbors_count >= limit) {
          break;
        }
        float val = reference_point[current_node->axis];
        float diff = val - current_node->axis_val;

        const kdtree_node * next;
        if (diff < 0 && current_node->child1) {
          next = current_node->child1;
          // printf("%d in (%d), select child1\n", item_index, n);
        } else {
          next = current_node->child2;
          // printf("%d in (%d), select child2\n", item_index, n);
        }
        current_node = next;
      }
    }
  }
}

/** \brief Compute point derivatives.
 * \note Equation 6.18-21 [Magnusson 2009].
 * \param[in] x point from the input cloud
 * \param[in] compute_hessian flag to calculate hessian, unnessissary for step calculation.
 * \param[in] j_ang precomputed angular gradient
 * \param[in] h_ang precomputed angular hessian
 */
void computePointDerivatives(const float x[4], float point_gradient_[4][6], float point_hessian_[24][6], const float j_ang[8][4], const float h_ang[16][4])
{
  float x4[4];
  x4[0] = x[0];
  x4[1] = x[1];
  x4[2] = x[2];
  x4[3] = 0.0f;

  // Calculate first derivative of Transformation Equation 6.17 w.r.t. transform vector p.
  // Derivative w.r.t. ith element of transform vector corresponds to column i, Equation 6.18 and 6.19 [Magnusson 2009]
  float x_j_ang[8];
  {
    x_j_ang[0] =
      j_ang[0][0] * x4[0] + j_ang[0][1] * x4[1] + j_ang[0][2] * x4[2] + j_ang[0][3] * x4[3];
    x_j_ang[1] =
      j_ang[1][0] * x4[0] + j_ang[1][1] * x4[1] + j_ang[1][2] * x4[2] + j_ang[1][3] * x4[3];
    x_j_ang[2] =
      j_ang[2][0] * x4[0] + j_ang[2][1] * x4[1] + j_ang[2][2] * x4[2] + j_ang[2][3] * x4[3];
    x_j_ang[3] =
      j_ang[3][0] * x4[0] + j_ang[3][1] * x4[1] + j_ang[3][2] * x4[2] + j_ang[3][3] * x4[3];
    x_j_ang[4] =
      j_ang[4][0] * x4[0] + j_ang[4][1] * x4[1] + j_ang[4][2] * x4[2] + j_ang[4][3] * x4[3];
    x_j_ang[5] =
      j_ang[5][0] * x4[0] + j_ang[5][1] * x4[1] + j_ang[5][2] * x4[2] + j_ang[5][3] * x4[3];
    x_j_ang[6] =
      j_ang[6][0] * x4[0] + j_ang[6][1] * x4[1] + j_ang[6][2] * x4[2] + j_ang[6][3] * x4[3];
    x_j_ang[7] =
      j_ang[7][0] * x4[0] + j_ang[7][1] * x4[1] + j_ang[7][2] * x4[2] + j_ang[7][3] * x4[3];
  }

  point_gradient_[1][3] = x_j_ang[0];
  point_gradient_[2][3] = x_j_ang[1];
  point_gradient_[0][4] = x_j_ang[2];
  point_gradient_[1][4] = x_j_ang[3];
  point_gradient_[2][4] = x_j_ang[4];
  point_gradient_[0][5] = x_j_ang[5];
  point_gradient_[1][5] = x_j_ang[6];
  point_gradient_[2][5] = x_j_ang[7];

  float x_h_ang[16];
  {
    x_h_ang[0] =
      h_ang[0][0] * x4[0] + h_ang[0][1] * x4[1] + h_ang[0][2] * x4[2] + h_ang[0][3] * x4[3];
    x_h_ang[1] =
      h_ang[1][0] * x4[0] + h_ang[1][1] * x4[1] + h_ang[1][2] * x4[2] + h_ang[1][3] * x4[3];
    x_h_ang[2] =
      h_ang[2][0] * x4[0] + h_ang[2][1] * x4[1] + h_ang[2][2] * x4[2] + h_ang[2][3] * x4[3];
    x_h_ang[3] =
      h_ang[3][0] * x4[0] + h_ang[3][1] * x4[1] + h_ang[3][2] * x4[2] + h_ang[3][3] * x4[3];
    x_h_ang[4] =
      h_ang[4][0] * x4[0] + h_ang[4][1] * x4[1] + h_ang[4][2] * x4[2] + h_ang[4][3] * x4[3];
    x_h_ang[5] =
      h_ang[5][0] * x4[0] + h_ang[5][1] * x4[1] + h_ang[5][2] * x4[2] + h_ang[5][3] * x4[3];
    x_h_ang[6] =
      h_ang[6][0] * x4[0] + h_ang[6][1] * x4[1] + h_ang[6][2] * x4[2] + h_ang[6][3] * x4[3];
    x_h_ang[7] =
      h_ang[7][0] * x4[0] + h_ang[7][1] * x4[1] + h_ang[7][2] * x4[2] + h_ang[7][3] * x4[3];
    x_h_ang[8] =
      h_ang[8][0] * x4[0] + h_ang[8][1] * x4[1] + h_ang[8][2] * x4[2] + h_ang[8][3] * x4[3];
    x_h_ang[9] =
      h_ang[9][0] * x4[0] + h_ang[9][1] * x4[1] + h_ang[9][2] * x4[2] + h_ang[9][3] * x4[3];
    x_h_ang[10] =
      h_ang[10][0] * x4[0] + h_ang[10][1] * x4[1] + h_ang[10][2] * x4[2] + h_ang[10][3] * x4[3];
    x_h_ang[11] =
      h_ang[11][0] * x4[0] + h_ang[11][1] * x4[1] + h_ang[11][2] * x4[2] + h_ang[11][3] * x4[3];
    x_h_ang[12] =
      h_ang[12][0] * x4[0] + h_ang[12][1] * x4[1] + h_ang[12][2] * x4[2] + h_ang[12][3] * x4[3];
    x_h_ang[13] =
      h_ang[13][0] * x4[0] + h_ang[13][1] * x4[1] + h_ang[13][2] * x4[2] + h_ang[13][3] * x4[3];
    x_h_ang[14] =
      h_ang[14][0] * x4[0] + h_ang[14][1] * x4[1] + h_ang[14][2] * x4[2] + h_ang[14][3] * x4[3];
    x_h_ang[15] =
      h_ang[15][0] * x4[0] + h_ang[15][1] * x4[1] + h_ang[15][2] * x4[2] + h_ang[15][3] * x4[3];
  }

  // Vectors from Equation 6.21 [Magnusson 2009]
  float a[4];
  float b[4];
  float c[4];
  float d[4];
  float e[4];
  float f[4];
  a[0] = 0.0f;
  a[1] = x_h_ang[0];
  a[2] = x_h_ang[1];
  a[3] = 0.0f;
  b[0] = 0.0f;
  b[1] = x_h_ang[2];
  b[2] = x_h_ang[3];
  b[3] = 0.0f;
  c[0] = 0.0f;
  c[1] = x_h_ang[4];
  c[2] = x_h_ang[5];
  c[3] = 0.0f;
  d[0] = x_h_ang[6];
  d[1] = x_h_ang[7];
  d[2] = x_h_ang[8];
  d[3] = 0.0f;
  e[0] = x_h_ang[9];
  e[1] = x_h_ang[10];
  e[2] = x_h_ang[11];
  e[3] = 0.0f;
  f[0] = x_h_ang[12];
  f[1] = x_h_ang[13];
  f[2] = x_h_ang[14];
  f[3] = x_h_ang[15];

  // Calculate second derivative of Transformation Equation 6.17 w.r.t. transform vector p.
  // Derivative w.r.t. ith and jth elements of transform vector corresponds to the 3x1 block matrix starting at
  // (3i,j), Equation 6.20 and 6.21 [Magnusson 2009]
  point_hessian_[12][3] = a[0];
  point_hessian_[13][3] = a[1];
  point_hessian_[14][3] = a[2];
  point_hessian_[15][3] = a[3];

  point_hessian_[16][3] = b[0];
  point_hessian_[17][3] = b[1];
  point_hessian_[18][3] = b[2];
  point_hessian_[19][3] = b[3];

  point_hessian_[20][3] = c[0];
  point_hessian_[21][3] = c[1];
  point_hessian_[22][3] = c[2];
  point_hessian_[23][3] = c[3];

  point_hessian_[12][4] = b[0];
  point_hessian_[13][4] = b[1];
  point_hessian_[14][4] = b[2];
  point_hessian_[15][4] = b[3];

  point_hessian_[16][4] = d[0];
  point_hessian_[17][4] = d[1];
  point_hessian_[18][4] = d[2];
  point_hessian_[19][4] = d[3];

  point_hessian_[20][4] = e[0];
  point_hessian_[21][4] = e[1];
  point_hessian_[22][4] = e[2];
  point_hessian_[23][4] = e[3];

  point_hessian_[12][5] = c[0];
  point_hessian_[13][5] = c[1];
  point_hessian_[14][5] = c[2];
  point_hessian_[15][5] = c[3];

  point_hessian_[16][5] = e[0];
  point_hessian_[17][5] = e[1];
  point_hessian_[18][5] = e[2];
  point_hessian_[19][5] = e[3];

  point_hessian_[20][5] = f[0];
  point_hessian_[21][5] = f[1];
  point_hessian_[22][5] = f[2];
  point_hessian_[23][5] = f[3];
}

float updateDerivatives(
  float score_gradient[6], float hessian[6][6], const float point_gradient4[4][6],
  const float point_hessian_[24][6], const float x_trans[3], const float c_inv[3][3], const float gauss_d1_,
  const float gauss_d2_, const float gauss_d3_)
{
  float x_trans4[4];
  x_trans4[0] = x_trans[0];
  x_trans4[1] = x_trans[1];
  x_trans4[2] = x_trans[2];
  x_trans4[3] = 0.0f;

  float c_inv4[4][4];
  c_inv4[0][0] = c_inv[0][0];
  c_inv4[0][1] = c_inv[0][1];
  c_inv4[0][2] = c_inv[0][2];
  c_inv4[0][3] = 0.0f;
  c_inv4[1][0] = c_inv[1][0];
  c_inv4[1][1] = c_inv[1][1];
  c_inv4[1][2] = c_inv[1][2];
  c_inv4[1][3] = 0.0f;
  c_inv4[2][0] = c_inv[2][0];
  c_inv4[2][1] = c_inv[2][1];
  c_inv4[2][2] = c_inv[2][2];
  c_inv4[2][3] = 0.0f;
  c_inv4[3][0] = 0.0f;
  c_inv4[3][1] = 0.0f;
  c_inv4[3][2] = 0.0f;
  c_inv4[3][3] = 0.0f;

  float x_trans4_x_c_inv4[4];
  x_trans4_x_c_inv4[0] = x_trans4[0] * c_inv4[0][0] + x_trans4[1] * c_inv4[0][1] +
                         x_trans4[2] * c_inv4[0][2] + x_trans4[3] * c_inv4[0][3];

  x_trans4_x_c_inv4[1] = x_trans4[0] * c_inv4[1][0] + x_trans4[1] * c_inv4[1][1] +
                         x_trans4[2] * c_inv4[1][2] + x_trans4[3] * c_inv4[1][3];

  x_trans4_x_c_inv4[2] = x_trans4[0] * c_inv4[2][0] + x_trans4[1] * c_inv4[2][1] +
                         x_trans4[2] * c_inv4[2][2] + x_trans4[3] * c_inv4[2][3];

  x_trans4_x_c_inv4[3] = x_trans4[0] * c_inv4[3][0] + x_trans4[1] * c_inv4[3][1] +
                         x_trans4[2] * c_inv4[3][2] + x_trans4[3] * c_inv4[3][3];

  float x_trans4_dot_x_trans4_x_c_inv4 =
    x_trans4[0] * x_trans4_x_c_inv4[0] + x_trans4[1] * x_trans4_x_c_inv4[1] +
    x_trans4[2] * x_trans4_x_c_inv4[2] + x_trans4[3] * x_trans4_x_c_inv4[3];
  // e^(-d_2/2 * (x_k - mu_k)^T Sigma_k^-1 (x_k - mu_k)) Equation 6.9 [Magnusson 2009]
  // float e_x_cov_x = exp(-gauss_d2_ * x_trans4_dot_x_trans4_x_c_inv4 * 0.5f);
  float e_x_cov_x = 2.7182282f * (-gauss_d2_ * x_trans4_dot_x_trans4_x_c_inv4 * 0.5f);
  // Calculate probability of transtormed points existance, Equation 6.9 [Magnusson 2009]
  float score_inc = -gauss_d1_ * e_x_cov_x;

  e_x_cov_x = gauss_d2_ * e_x_cov_x;

  // Error checking for invalid values.
  if (e_x_cov_x > 1.0 || e_x_cov_x < 0.0 || e_x_cov_x != e_x_cov_x) return 0;

  // Reusable portion of Equation 6.12 and 6.13 [Magnusson 2009]
  e_x_cov_x *= gauss_d1_;

  float c_inv4_x_point_gradient4[4][6];
  c_inv4_x_point_gradient4[0][0] =
    c_inv4[0][0] * point_gradient4[0][0] + c_inv4[0][1] * point_gradient4[1][0] +
    c_inv4[0][2] * point_gradient4[2][0] + c_inv4[0][3] * point_gradient4[3][0];

  c_inv4_x_point_gradient4[0][1] =
    c_inv4[0][0] * point_gradient4[0][1] + c_inv4[0][1] * point_gradient4[1][1] +
    c_inv4[0][2] * point_gradient4[2][1] + c_inv4[0][3] * point_gradient4[3][1];

  c_inv4_x_point_gradient4[0][2] =
    c_inv4[0][0] * point_gradient4[0][2] + c_inv4[0][1] * point_gradient4[1][2] +
    c_inv4[0][2] * point_gradient4[2][2] + c_inv4[0][3] * point_gradient4[3][2];

  c_inv4_x_point_gradient4[0][3] =
    c_inv4[0][0] * point_gradient4[0][3] + c_inv4[0][1] * point_gradient4[1][3] +
    c_inv4[0][2] * point_gradient4[2][3] + c_inv4[0][3] * point_gradient4[3][3];

  c_inv4_x_point_gradient4[0][4] =
    c_inv4[0][0] * point_gradient4[0][4] + c_inv4[0][1] * point_gradient4[1][4] +
    c_inv4[0][2] * point_gradient4[2][4] + c_inv4[0][3] * point_gradient4[3][4];

  c_inv4_x_point_gradient4[0][5] =
    c_inv4[0][0] * point_gradient4[0][5] + c_inv4[0][1] * point_gradient4[1][5] +
    c_inv4[0][2] * point_gradient4[2][5] + c_inv4[0][3] * point_gradient4[3][5];

  ////

  c_inv4_x_point_gradient4[1][0] =
    c_inv4[1][0] * point_gradient4[0][0] + c_inv4[1][1] * point_gradient4[1][0] +
    c_inv4[1][2] * point_gradient4[2][0] + c_inv4[1][3] * point_gradient4[3][0];

  c_inv4_x_point_gradient4[1][1] =
    c_inv4[1][0] * point_gradient4[0][1] + c_inv4[1][1] * point_gradient4[1][1] +
    c_inv4[1][2] * point_gradient4[2][1] + c_inv4[1][3] * point_gradient4[3][1];

  c_inv4_x_point_gradient4[1][2] =
    c_inv4[1][0] * point_gradient4[0][2] + c_inv4[1][1] * point_gradient4[1][2] +
    c_inv4[1][2] * point_gradient4[2][2] + c_inv4[1][3] * point_gradient4[3][2];

  c_inv4_x_point_gradient4[1][3] =
    c_inv4[1][0] * point_gradient4[0][3] + c_inv4[1][1] * point_gradient4[1][3] +
    c_inv4[1][2] * point_gradient4[2][3] + c_inv4[1][3] * point_gradient4[3][3];

  c_inv4_x_point_gradient4[1][4] =
    c_inv4[1][0] * point_gradient4[0][4] + c_inv4[1][1] * point_gradient4[1][4] +
    c_inv4[1][2] * point_gradient4[2][4] + c_inv4[1][3] * point_gradient4[3][4];

  c_inv4_x_point_gradient4[1][5] =
    c_inv4[1][0] * point_gradient4[0][5] + c_inv4[1][1] * point_gradient4[1][5] +
    c_inv4[1][2] * point_gradient4[2][5] + c_inv4[1][3] * point_gradient4[3][5];

  ////

  c_inv4_x_point_gradient4[2][0] =
    c_inv4[2][0] * point_gradient4[0][0] + c_inv4[2][1] * point_gradient4[1][0] +
    c_inv4[2][2] * point_gradient4[2][0] + c_inv4[2][3] * point_gradient4[3][0];

  c_inv4_x_point_gradient4[2][1] =
    c_inv4[2][0] * point_gradient4[0][1] + c_inv4[2][1] * point_gradient4[1][1] +
    c_inv4[2][2] * point_gradient4[2][1] + c_inv4[2][3] * point_gradient4[3][1];

  c_inv4_x_point_gradient4[2][2] =
    c_inv4[2][0] * point_gradient4[0][2] + c_inv4[2][1] * point_gradient4[1][2] +
    c_inv4[2][2] * point_gradient4[2][2] + c_inv4[2][3] * point_gradient4[3][2];

  c_inv4_x_point_gradient4[2][3] =
    c_inv4[2][0] * point_gradient4[0][3] + c_inv4[2][1] * point_gradient4[1][3] +
    c_inv4[2][2] * point_gradient4[2][3] + c_inv4[2][3] * point_gradient4[3][3];

  c_inv4_x_point_gradient4[2][4] =
    c_inv4[2][0] * point_gradient4[0][4] + c_inv4[2][1] * point_gradient4[1][4] +
    c_inv4[2][2] * point_gradient4[2][4] + c_inv4[2][3] * point_gradient4[3][4];

  c_inv4_x_point_gradient4[2][5] =
    c_inv4[2][0] * point_gradient4[0][5] + c_inv4[2][1] * point_gradient4[1][5] +
    c_inv4[2][2] * point_gradient4[2][5] + c_inv4[2][3] * point_gradient4[3][5];

  ////

  c_inv4_x_point_gradient4[3][0] =
    c_inv4[3][0] * point_gradient4[0][0] + c_inv4[3][1] * point_gradient4[1][0] +
    c_inv4[3][2] * point_gradient4[2][0] + c_inv4[3][3] * point_gradient4[3][0];

  c_inv4_x_point_gradient4[3][1] =
    c_inv4[3][0] * point_gradient4[0][1] + c_inv4[3][1] * point_gradient4[1][1] +
    c_inv4[3][2] * point_gradient4[2][1] + c_inv4[3][3] * point_gradient4[3][1];

  c_inv4_x_point_gradient4[3][2] =
    c_inv4[3][0] * point_gradient4[0][2] + c_inv4[3][1] * point_gradient4[1][2] +
    c_inv4[3][2] * point_gradient4[2][2] + c_inv4[3][3] * point_gradient4[3][2];

  c_inv4_x_point_gradient4[3][3] =
    c_inv4[3][0] * point_gradient4[0][3] + c_inv4[3][1] * point_gradient4[1][3] +
    c_inv4[3][2] * point_gradient4[2][3] + c_inv4[3][3] * point_gradient4[3][3];

  c_inv4_x_point_gradient4[3][4] =
    c_inv4[3][0] * point_gradient4[0][4] + c_inv4[3][1] * point_gradient4[1][4] +
    c_inv4[3][2] * point_gradient4[2][4] + c_inv4[3][3] * point_gradient4[3][4];

  c_inv4_x_point_gradient4[3][5] =
    c_inv4[3][0] * point_gradient4[0][5] + c_inv4[3][1] * point_gradient4[1][5] +
    c_inv4[3][2] * point_gradient4[2][5] + c_inv4[3][3] * point_gradient4[3][5];

  ////

  float x_trans4_dot_c_inv4_x_point_gradient4[6];
  x_trans4_dot_c_inv4_x_point_gradient4[0] =
    x_trans4[0] * c_inv4_x_point_gradient4[0][0] + x_trans4[1] * c_inv4_x_point_gradient4[1][0] +
    x_trans4[2] * c_inv4_x_point_gradient4[2][0] + x_trans4[3] * c_inv4_x_point_gradient4[3][0];

  x_trans4_dot_c_inv4_x_point_gradient4[1] =
    x_trans4[0] * c_inv4_x_point_gradient4[0][1] + x_trans4[1] * c_inv4_x_point_gradient4[1][1] +
    x_trans4[2] * c_inv4_x_point_gradient4[2][1] + x_trans4[3] * c_inv4_x_point_gradient4[3][1];

  x_trans4_dot_c_inv4_x_point_gradient4[2] =
    x_trans4[0] * c_inv4_x_point_gradient4[0][2] + x_trans4[1] * c_inv4_x_point_gradient4[1][2] +
    x_trans4[2] * c_inv4_x_point_gradient4[2][2] + x_trans4[3] * c_inv4_x_point_gradient4[3][2];

  x_trans4_dot_c_inv4_x_point_gradient4[3] =
    x_trans4[0] * c_inv4_x_point_gradient4[0][3] + x_trans4[1] * c_inv4_x_point_gradient4[1][3] +
    x_trans4[2] * c_inv4_x_point_gradient4[2][3] + x_trans4[3] * c_inv4_x_point_gradient4[3][3];

  x_trans4_dot_c_inv4_x_point_gradient4[4] =
    x_trans4[0] * c_inv4_x_point_gradient4[0][4] + x_trans4[1] * c_inv4_x_point_gradient4[1][4] +
    x_trans4[2] * c_inv4_x_point_gradient4[2][4] + x_trans4[3] * c_inv4_x_point_gradient4[3][4];

  x_trans4_dot_c_inv4_x_point_gradient4[5] =
    x_trans4[0] * c_inv4_x_point_gradient4[0][5] + x_trans4[1] * c_inv4_x_point_gradient4[1][5] +
    x_trans4[2] * c_inv4_x_point_gradient4[2][5] + x_trans4[3] * c_inv4_x_point_gradient4[3][5];

  score_gradient[0] += e_x_cov_x * x_trans4_dot_c_inv4_x_point_gradient4[0];
  score_gradient[1] += e_x_cov_x * x_trans4_dot_c_inv4_x_point_gradient4[1];
  score_gradient[2] += e_x_cov_x * x_trans4_dot_c_inv4_x_point_gradient4[2];
  score_gradient[3] += e_x_cov_x * x_trans4_dot_c_inv4_x_point_gradient4[3];
  score_gradient[4] += e_x_cov_x * x_trans4_dot_c_inv4_x_point_gradient4[4];
  score_gradient[5] += e_x_cov_x * x_trans4_dot_c_inv4_x_point_gradient4[5];

  ///
  float point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[6][6];
  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[0][0] =
    point_gradient4[0][0] * c_inv4_x_point_gradient4[0][0] +
    point_gradient4[1][0] * c_inv4_x_point_gradient4[1][0] +
    point_gradient4[2][0] * c_inv4_x_point_gradient4[2][0] +
    point_gradient4[3][0] * c_inv4_x_point_gradient4[3][0];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[0][1] =
    point_gradient4[0][0] * c_inv4_x_point_gradient4[0][1] +
    point_gradient4[1][0] * c_inv4_x_point_gradient4[1][1] +
    point_gradient4[2][0] * c_inv4_x_point_gradient4[2][1] +
    point_gradient4[3][0] * c_inv4_x_point_gradient4[3][1];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[0][2] =
    point_gradient4[0][0] * c_inv4_x_point_gradient4[0][2] +
    point_gradient4[1][0] * c_inv4_x_point_gradient4[1][2] +
    point_gradient4[2][0] * c_inv4_x_point_gradient4[2][2] +
    point_gradient4[3][0] * c_inv4_x_point_gradient4[3][2];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[0][3] =
    point_gradient4[0][0] * c_inv4_x_point_gradient4[0][3] +
    point_gradient4[1][0] * c_inv4_x_point_gradient4[1][3] +
    point_gradient4[2][0] * c_inv4_x_point_gradient4[2][3] +
    point_gradient4[3][0] * c_inv4_x_point_gradient4[3][3];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[0][4] =
    point_gradient4[0][0] * c_inv4_x_point_gradient4[0][4] +
    point_gradient4[1][0] * c_inv4_x_point_gradient4[1][4] +
    point_gradient4[2][0] * c_inv4_x_point_gradient4[2][4] +
    point_gradient4[3][0] * c_inv4_x_point_gradient4[3][4];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[0][5] =
    point_gradient4[0][0] * c_inv4_x_point_gradient4[0][5] +
    point_gradient4[1][0] * c_inv4_x_point_gradient4[1][5] +
    point_gradient4[2][0] * c_inv4_x_point_gradient4[2][5] +
    point_gradient4[3][0] * c_inv4_x_point_gradient4[3][5];

  ////

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[1][0] =
    point_gradient4[0][1] * c_inv4_x_point_gradient4[0][0] +
    point_gradient4[1][1] * c_inv4_x_point_gradient4[1][0] +
    point_gradient4[2][1] * c_inv4_x_point_gradient4[2][0] +
    point_gradient4[3][1] * c_inv4_x_point_gradient4[3][0];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[1][1] =
    point_gradient4[0][1] * c_inv4_x_point_gradient4[0][1] +
    point_gradient4[1][1] * c_inv4_x_point_gradient4[1][1] +
    point_gradient4[2][1] * c_inv4_x_point_gradient4[2][1] +
    point_gradient4[3][1] * c_inv4_x_point_gradient4[3][1];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[1][2] =
    point_gradient4[0][1] * c_inv4_x_point_gradient4[0][2] +
    point_gradient4[1][1] * c_inv4_x_point_gradient4[1][2] +
    point_gradient4[2][1] * c_inv4_x_point_gradient4[2][2] +
    point_gradient4[3][1] * c_inv4_x_point_gradient4[3][2];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[1][3] =
    point_gradient4[0][1] * c_inv4_x_point_gradient4[0][3] +
    point_gradient4[1][1] * c_inv4_x_point_gradient4[1][3] +
    point_gradient4[2][1] * c_inv4_x_point_gradient4[2][3] +
    point_gradient4[3][1] * c_inv4_x_point_gradient4[3][3];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[1][4] =
    point_gradient4[0][1] * c_inv4_x_point_gradient4[0][4] +
    point_gradient4[1][1] * c_inv4_x_point_gradient4[1][4] +
    point_gradient4[2][1] * c_inv4_x_point_gradient4[2][4] +
    point_gradient4[3][1] * c_inv4_x_point_gradient4[3][4];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[1][5] =
    point_gradient4[0][1] * c_inv4_x_point_gradient4[0][5] +
    point_gradient4[1][1] * c_inv4_x_point_gradient4[1][5] +
    point_gradient4[2][1] * c_inv4_x_point_gradient4[2][5] +
    point_gradient4[3][1] * c_inv4_x_point_gradient4[3][5];

  ////

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[2][0] =
    point_gradient4[0][2] * c_inv4_x_point_gradient4[0][0] +
    point_gradient4[1][2] * c_inv4_x_point_gradient4[1][0] +
    point_gradient4[2][2] * c_inv4_x_point_gradient4[2][0] +
    point_gradient4[3][2] * c_inv4_x_point_gradient4[3][0];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[2][1] =
    point_gradient4[0][2] * c_inv4_x_point_gradient4[0][1] +
    point_gradient4[1][2] * c_inv4_x_point_gradient4[1][1] +
    point_gradient4[2][2] * c_inv4_x_point_gradient4[2][1] +
    point_gradient4[3][2] * c_inv4_x_point_gradient4[3][1];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[2][2] =
    point_gradient4[0][2] * c_inv4_x_point_gradient4[0][2] +
    point_gradient4[1][2] * c_inv4_x_point_gradient4[1][2] +
    point_gradient4[2][2] * c_inv4_x_point_gradient4[2][2] +
    point_gradient4[3][2] * c_inv4_x_point_gradient4[3][2];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[2][3] =
    point_gradient4[0][2] * c_inv4_x_point_gradient4[0][3] +
    point_gradient4[1][2] * c_inv4_x_point_gradient4[1][3] +
    point_gradient4[2][2] * c_inv4_x_point_gradient4[2][3] +
    point_gradient4[3][2] * c_inv4_x_point_gradient4[3][3];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[2][4] =
    point_gradient4[0][2] * c_inv4_x_point_gradient4[0][4] +
    point_gradient4[1][2] * c_inv4_x_point_gradient4[1][4] +
    point_gradient4[2][2] * c_inv4_x_point_gradient4[2][4] +
    point_gradient4[3][2] * c_inv4_x_point_gradient4[3][4];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[2][5] =
    point_gradient4[0][2] * c_inv4_x_point_gradient4[0][5] +
    point_gradient4[1][2] * c_inv4_x_point_gradient4[1][5] +
    point_gradient4[2][2] * c_inv4_x_point_gradient4[2][5] +
    point_gradient4[3][2] * c_inv4_x_point_gradient4[3][5];

  ////

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[3][0] =
    point_gradient4[0][3] * c_inv4_x_point_gradient4[0][0] +
    point_gradient4[1][3] * c_inv4_x_point_gradient4[1][0] +
    point_gradient4[2][3] * c_inv4_x_point_gradient4[2][0] +
    point_gradient4[3][3] * c_inv4_x_point_gradient4[3][0];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[3][1] =
    point_gradient4[0][3] * c_inv4_x_point_gradient4[0][1] +
    point_gradient4[1][3] * c_inv4_x_point_gradient4[1][1] +
    point_gradient4[2][3] * c_inv4_x_point_gradient4[2][1] +
    point_gradient4[3][3] * c_inv4_x_point_gradient4[3][1];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[3][2] =
    point_gradient4[0][3] * c_inv4_x_point_gradient4[0][2] +
    point_gradient4[1][3] * c_inv4_x_point_gradient4[1][2] +
    point_gradient4[2][3] * c_inv4_x_point_gradient4[2][2] +
    point_gradient4[3][3] * c_inv4_x_point_gradient4[3][2];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[3][3] =
    point_gradient4[0][3] * c_inv4_x_point_gradient4[0][3] +
    point_gradient4[1][3] * c_inv4_x_point_gradient4[1][3] +
    point_gradient4[2][3] * c_inv4_x_point_gradient4[2][3] +
    point_gradient4[3][3] * c_inv4_x_point_gradient4[3][3];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[3][4] =
    point_gradient4[0][3] * c_inv4_x_point_gradient4[0][4] +
    point_gradient4[1][3] * c_inv4_x_point_gradient4[1][4] +
    point_gradient4[2][3] * c_inv4_x_point_gradient4[2][4] +
    point_gradient4[3][3] * c_inv4_x_point_gradient4[3][4];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[3][5] =
    point_gradient4[0][3] * c_inv4_x_point_gradient4[0][5] +
    point_gradient4[1][3] * c_inv4_x_point_gradient4[1][5] +
    point_gradient4[2][3] * c_inv4_x_point_gradient4[2][5] +
    point_gradient4[3][3] * c_inv4_x_point_gradient4[3][5];

  ////

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[4][0] =
    point_gradient4[0][4] * c_inv4_x_point_gradient4[0][0] +
    point_gradient4[1][4] * c_inv4_x_point_gradient4[1][0] +
    point_gradient4[2][4] * c_inv4_x_point_gradient4[2][0] +
    point_gradient4[3][4] * c_inv4_x_point_gradient4[3][0];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[4][1] =
    point_gradient4[0][4] * c_inv4_x_point_gradient4[0][1] +
    point_gradient4[1][4] * c_inv4_x_point_gradient4[1][1] +
    point_gradient4[2][4] * c_inv4_x_point_gradient4[2][1] +
    point_gradient4[3][4] * c_inv4_x_point_gradient4[3][1];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[4][2] =
    point_gradient4[0][4] * c_inv4_x_point_gradient4[0][2] +
    point_gradient4[1][4] * c_inv4_x_point_gradient4[1][2] +
    point_gradient4[2][4] * c_inv4_x_point_gradient4[2][2] +
    point_gradient4[3][4] * c_inv4_x_point_gradient4[3][2];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[4][3] =
    point_gradient4[0][4] * c_inv4_x_point_gradient4[0][3] +
    point_gradient4[1][4] * c_inv4_x_point_gradient4[1][3] +
    point_gradient4[2][4] * c_inv4_x_point_gradient4[2][3] +
    point_gradient4[3][4] * c_inv4_x_point_gradient4[3][3];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[4][4] =
    point_gradient4[0][4] * c_inv4_x_point_gradient4[0][4] +
    point_gradient4[1][4] * c_inv4_x_point_gradient4[1][4] +
    point_gradient4[2][4] * c_inv4_x_point_gradient4[2][4] +
    point_gradient4[3][4] * c_inv4_x_point_gradient4[3][4];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[4][5] =
    point_gradient4[0][4] * c_inv4_x_point_gradient4[0][5] +
    point_gradient4[1][4] * c_inv4_x_point_gradient4[1][5] +
    point_gradient4[2][4] * c_inv4_x_point_gradient4[2][5] +
    point_gradient4[3][4] * c_inv4_x_point_gradient4[3][5];

  ////

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[5][0] =
    point_gradient4[0][5] * c_inv4_x_point_gradient4[0][0] +
    point_gradient4[1][5] * c_inv4_x_point_gradient4[1][0] +
    point_gradient4[2][5] * c_inv4_x_point_gradient4[2][0] +
    point_gradient4[3][5] * c_inv4_x_point_gradient4[3][0];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[5][1] =
    point_gradient4[0][5] * c_inv4_x_point_gradient4[0][1] +
    point_gradient4[1][5] * c_inv4_x_point_gradient4[1][1] +
    point_gradient4[2][5] * c_inv4_x_point_gradient4[2][1] +
    point_gradient4[3][5] * c_inv4_x_point_gradient4[3][1];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[5][2] =
    point_gradient4[0][5] * c_inv4_x_point_gradient4[0][2] +
    point_gradient4[1][5] * c_inv4_x_point_gradient4[1][2] +
    point_gradient4[2][5] * c_inv4_x_point_gradient4[2][2] +
    point_gradient4[3][5] * c_inv4_x_point_gradient4[3][2];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[5][3] =
    point_gradient4[0][5] * c_inv4_x_point_gradient4[0][3] +
    point_gradient4[1][5] * c_inv4_x_point_gradient4[1][3] +
    point_gradient4[2][5] * c_inv4_x_point_gradient4[2][3] +
    point_gradient4[3][5] * c_inv4_x_point_gradient4[3][3];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[5][4] =
    point_gradient4[0][5] * c_inv4_x_point_gradient4[0][4] +
    point_gradient4[1][5] * c_inv4_x_point_gradient4[1][4] +
    point_gradient4[2][5] * c_inv4_x_point_gradient4[2][4] +
    point_gradient4[3][5] * c_inv4_x_point_gradient4[3][4];

  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[5][5] =
    point_gradient4[0][5] * c_inv4_x_point_gradient4[0][5] +
    point_gradient4[1][5] * c_inv4_x_point_gradient4[1][5] +
    point_gradient4[2][5] * c_inv4_x_point_gradient4[2][5] +
    point_gradient4[3][5] * c_inv4_x_point_gradient4[3][5];

  float x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[6];

  {  // i = 0;
    int i = 0;
    ///
    int j = 0;
    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[0] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][0] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][0] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][0] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][0];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[1] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][1] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][1] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][1] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][1];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[2] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][2] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][2] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][2] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][2];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[3] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][3] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][3] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][3] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][3];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[4] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][4] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][4] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][4] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][4];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[5] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][5] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][5] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][5] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][5];

    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    ////
    j = 1;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 2;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 3;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 4;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 5;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
  }

  {  // i = 1;
    int i = 1;
    ///
    int j = 0;
    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[0] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][0] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][0] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][0] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][0];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[1] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][1] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][1] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][1] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][1];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[2] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][2] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][2] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][2] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][2];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[3] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][3] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][3] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][3] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][3];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[4] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][4] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][4] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][4] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][4];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[5] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][5] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][5] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][5] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][5];

    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    ////
    j = 2;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 3;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 4;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 5;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 6;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
  }

  {  // i = 2;
    int i = 2;
    ///
    int j = 0;
    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[0] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][0] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][0] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][0] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][0];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[1] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][1] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][1] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][1] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][1];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[2] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][2] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][2] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][2] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][2];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[3] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][3] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][3] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][3] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][3];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[4] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][4] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][4] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][4] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][4];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[5] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][5] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][5] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][5] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][5];

    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    ////
    j = 2;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 3;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 4;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 5;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 6;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
  }

  {  // i = 3;
    int i = 3;
    ///
    int j = 0;
    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[0] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][0] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][0] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][0] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][0];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[1] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][1] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][1] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][1] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][1];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[2] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][2] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][2] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][2] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][2];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[3] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][3] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][3] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][3] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][3];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[4] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][4] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][4] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][4] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][4];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[5] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][5] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][5] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][5] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][5];

    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    ////
    j = 1;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 2;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 3;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 4;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 5;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
  }

  {  // i = 4;
    int i = 4;
    ///
    int j = 0;
    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[0] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][0] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][0] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][0] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][0];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[1] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][1] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][1] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][1] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][1];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[2] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][2] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][2] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][2] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][2];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[3] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][3] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][3] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][3] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][3];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[4] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][4] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][4] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][4] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][4];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[5] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][5] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][5] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][5] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][5];

    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    ////
    j = 2;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 3;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 4;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 5;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 6;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
  }

  {  // i = 5;
    int i = 5;
    ///
    int j = 0;
    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[0] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][0] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][0] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][0] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][0];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[1] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][1] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][1] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][1] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][1];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[2] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][2] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][2] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][2] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][2];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[3] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][3] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][3] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][3] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][3];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[4] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][4] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][4] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][4] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][4];

    x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[5] =
      x_trans4_x_c_inv4[0] * point_hessian_[i * 4 + 0][5] +
      x_trans4_x_c_inv4[1] * point_hessian_[i * 4 + 1][5] +
      x_trans4_x_c_inv4[2] * point_hessian_[i * 4 + 2][5] +
      x_trans4_x_c_inv4[3] * point_hessian_[i * 4 + 3][5];

    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    ////
    j = 2;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 3;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 4;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 5;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
    j = 6;
    hessian[i][j] += e_x_cov_x * (-gauss_d2_ * x_trans4_dot_c_inv4_x_point_gradient4[i] *
                                    x_trans4_dot_c_inv4_x_point_gradient4[j] +
                                  x_trans4_dot_c_inv4_x_ext_point_hessian_4ij[j] +
                                  point_gradient4_colj_dot_c_inv4_x_point_gradient4_col_i[j][i]);
  }

  return score_inc;
}

typedef struct leaf
{
  int nr_points;
  float mean[3];
  float * centroid;
  int centroid_size;
  float cov[3][3];
  float icov[3][3];
  float evecs[3][3];
  float evals[3];
} leaf;

/** \brief Compute derivatives of probability function w.r.t. the transformation vector.
   * \note Equation 6.10, 6.12 and 6.13 [Magnusson 2009].
   * \param[out] score_gradient the gradient vector of the probability function w.r.t. the transformation vector
   * \param[out] hessian the hessian matrix of the probability function w.r.t. the transformation vector
   * \param[in] trans_cloud_x transformed point cloud x
   * \param[in] trans_cloud_y transformed point cloud y
   * \param[in] trans_cloud_z transformed point cloud z
   * \param[in] p the current transform vector
   * \param[in] input_points_x the input point cloud dataset x
   * \param[in] input_points_y the input point cloud dataset y
   * \param[in] input_points_z the input point cloud dataset z
   * \param[in] input_points_size the input point cloud dataset size
   */
float computeDerivatives(
  float score_gradient[6], float hessian[6][6], float * trans_cloud_x, float * trans_cloud_y,
  float * trans_cloud_z, float p[6], float * input_points_x, float * input_points_y,
  float * input_points_z, int input_points_size)
{
  memset(score_gradient, 0, sizeof(float) * 6);
  memset(hessian, 0, sizeof(float) * 6 * 6);
  float score = 0;

  float score_gradients[6] = {};
  float hessians[6][6] = {};

  // Precompute Angular Derivatives (eq. 6.19 and 6.21)[Magnusson 2009]
  float j_ang_a[3];
  float j_ang_b[3];
  float j_ang_c[3];
  float j_ang_d[3];
  float j_ang_e[3];
  float j_ang_f[3];
  float j_ang_g[3];
  float j_ang_h[3];
  float j_ang[8][4];
  float h_ang_a2[3];
  float h_ang_a3[3];
  float h_ang_b2[3];
  float h_ang_b3[3];
  float h_ang_c2[3];
  float h_ang_c3[3];
  float h_ang_d1[3];
  float h_ang_d2[3];
  float h_ang_d3[3];
  float h_ang_e1[3];
  float h_ang_e2[3];
  float h_ang_e3[3];
  float h_ang_f1[3];
  float h_ang_f2[3];
  float h_ang_f3[3];
  float h_ang[16][4];
  // TODO: Supposed to call kernel function
  computeAngleDerivatives(
    p, j_ang_a, j_ang_b, j_ang_c, j_ang_d, j_ang_e, j_ang_f, j_ang_g, j_ang_h, j_ang, h_ang_a2,
    h_ang_a3, h_ang_b2, h_ang_b3, h_ang_c2, h_ang_c3, h_ang_d1, h_ang_d2, h_ang_d3, h_ang_e1,
    h_ang_e2, h_ang_e3, h_ang_f1, h_ang_f2, h_ang_f3, h_ang);

  std::vector<leaf> neighborhoods;
  std::vector<float> distancess;

  // Update gradient and hessian for each point, line 17 in Algorithm 2 [Magnusson 2009]
  for (int idx = 0; idx < input_points_size; idx++) {
    // Original Point and Transformed Point
    float x_pt[3];
    float x_trans_pt[3];
    // Original Point and Transformed Point (for math)
    float x[3];
    float x_trans[3];
    // Occupied Voxel
    leaf cell;
    // Inverse Covariance of Occupied Voxel
    float c_inv[3][3];

    // Initialize Point Gradient and Hessian
    float point_gradient_[4][6] = {
      {1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
      {0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f},
      {0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f},
      {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}};
    float point_hessian_[24][6] = {};

    x_trans_pt[0] = trans_cloud_x[idx];
    x_trans_pt[1] = trans_cloud_y[idx];
    x_trans_pt[2] = trans_cloud_z[idx];

    int leaf_size = 0;
    //target_cells_.radiusSearch(x_trans_pt, resolution_, neighborhoods[idx], distances[idx]);
    // TODO: Supposed to call kernel function
    //radiusSearchC()

    float score_pt = 0;
    float score_gradient_pt[6] = {};
    float hessian_pt[6][6] = {};

    for (int i = 0; i < leaf_size; ++i) {
      x_pt[0] = input_points_x[idx];
      x_pt[1] = input_points_y[idx];
      x_pt[2] = input_points_z[idx];
      x[0] = x_pt[0];
      x[1] = x_pt[1];
      x[2] = x_pt[2];
      x_trans[0] = x_trans_pt[0];
      x_trans[1] = x_trans_pt[1];
      x_trans[2] = x_trans_pt[2];

      // Denorm point, x_k' in Equations 6.12 and 6.13 [Magnusson 2009]
      x_trans[0] -= neighborhoods[i].mean[0];
      x_trans[1] -= neighborhoods[i].mean[1];
      x_trans[2] -= neighborhoods[i].mean[2];
      // Uses precomputed covariance for speed.
      memcpy(c_inv, neighborhoods[i].icov, sizeof(c_inv));

      // Compute derivative of transform function w.r.t. transform vector, J_E
      // and H_E in Equations 6.18 and 6.20 [Magnusson 2009]
      computePointDerivatives(x, point_gradient_, point_hessian_, j_ang, h_ang);
      // Update score, gradient and hessian, lines 19-21 in Algorithm 2,
      // according to Equations 6.10, 6.12 and 6.13, respectively [Magnusson
      // 2009]
      float gauss_d1_;
      float gauss_d2_;
      float gauss_d3_;
      // TODO: Supposed to call kernel function
      score_pt += updateDerivatives(
        score_gradient_pt, hessian_pt, point_gradient_, point_hessian_,
        x_trans, c_inv, gauss_d1_, gauss_d2_, gauss_d3_);
    }
    score += score_pt;

    int j = 0;
    for (j = 0; j < 6; ++j) score_gradients[j] += score_gradient_pt[j];
    for (j = 0; j < 6; ++j) {
      for (int k = 0; k < 6; ++k) hessians[j][k] += hessian_pt[j][k];
    }
  }

  return (score);
}
