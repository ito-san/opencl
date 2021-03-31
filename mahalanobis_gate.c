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
/**
 * @brief check whether a measurement value falls within the mahalanobis distance threshold
 * @param dist_max mahalanobis distance threshold
 * @param estimated current estimated state
 * @param measured measured state
 * @param estimated_cov current estimation covariance
 * @return whether it falls within the mahalanobis distance threshold
 */
int mahalanobisGate(
  const float dist_max, const float x[2], const float obj_x[2], const float cov[2][2])
{
  float a[2] = {x[0] - obj_x[0], x[1] - obj_x[1]};
  float det = cov[0][0] * cov[1][1] - cov[0][1] * cov[1][0];
  float inv[2][2];
  inv[0][0] = cov[1][1] / det;
  inv[0][1] = cov[0][1] * (-1.0f) / det;
  inv[1][0] = cov[1][0] * (-1.0f) / det;
  inv[1][1] = cov[0][0] / det;
  float b[2] = {a[0] * inv[0][0] + a[1] * inv[0][1], a[0] * inv[1][0] + a[1] * inv[1][1]};
  float mahalanobis_squared = b[0] * a[0] + b[1] * a[1];

  //printf(
  //  "!measurement update: mahalanobis = %f, gate limit = %lf", sqrt(mahalanobis_squared), dist_max);
  if (mahalanobis_squared > dist_max * dist_max) {
    return 0;
  }

  return 1;
}
