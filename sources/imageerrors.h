#ifndef IMAGEERRORS_H
#define IMAGEERRORS_H

#include <math.h>
#include <stdio.h>

#include "format/format.h"
#include "config.h"

enum EErrorMetric {
  EMSE = 0,
  ERMSE = 1,
  EMSE_LOG = 2,
  ERMSE_LOG = 3,
  ETVI = 4,
  ERelative = 5,
  ERelMSE = 6,
  ESMAPE = 7,
  ENUMBER_METRIC = 8,
};

double luminance(double R, double G, double B) {
  return 0.212671f * R + 0.715160f * G + 0.072169f * B;
}

void falseColor(float v, float vmin, float vmax, unsigned char& red,
                unsigned char& green, unsigned char& blue) {
  float Fred = 1.f, Fgreen = 1.f, Fblue = 1.f;

  // Compute vMin, vMax range
  if (v < vmin)
    v = vmin;
  if (v > vmax)
    v = vmax;
  float dv = vmax - vmin;

  // Color ramp manipulation
  if (v < (vmin + 0.25 * dv)) {
    Fred = 0;
    Fgreen = 4 * (v - vmin) / dv;
  } else if (v < (vmin + 0.5 * dv)) {
    Fred = 0;
    Fblue = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
  } else if (v < (vmin + 0.75 * dv)) {
    Fred = 4 * (v - vmin - 0.5 * dv) / dv;
    Fblue = 0;
  } else {
    Fgreen = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
    Fblue = 0;
  }

  // Scale into 255
  red = Fred * 255;
  green = Fgreen * 255;
  blue = Fblue * 255;
}

float tvi(const float Y) {
  float logY = log10f(std::max(0.0f, Y));
  float logResult = 0;

  if (logY < -3.94f) {
    logResult = -2.86f;
  }  // end if
  else if (logY < -1.44f) {
    logResult = powf(0.405f * logY + 1.6f, 2.18f) - 2.86f;
  }  // end else if
  else if (logY < -0.0184f) {
    logResult = logY - 0.395f;
  }  // end else
  else if (logY < 1.9f) {
    logResult = powf(0.249f * logY + 0.65f, 2.7f) - 0.72f;
  }  // end else if
  else {
    logResult = logY - 1.255f;
  }  // end else

  return powf(10.0f, logResult);
}

static float metricPix(RGB& pix1, RGB& pix2,
                       EErrorMetric metric) {
    float Rdiff = 0.f, Gdiff = 0.f, Bdiff = 0.f;

    float r1 = std::get<0>(pix1), r2 = std::get<0>(pix2);
    float g1 = std::get<1>(pix1), g2 = std::get<1>(pix2);
    float b1 = std::get<2>(pix1), b2 = std::get<2>(pix2);


  // Compute the pixel differences
  if (metric == EMSE || metric == ERMSE ||
    metric == ERelative || metric == ETVI ||
    metric == ERelMSE || metric == ESMAPE) {
    // Difference between the pixels
    Rdiff = (r1 - r2);
    Gdiff = (g1 - g2);
    Bdiff = (b1 - b2);
  } else if (metric == EMSE_LOG || metric == ERMSE_LOG) {
    // Logarithm difference
    Rdiff = log(r1 + 0.00001f) - log(r2 + 0.00001f);
    Gdiff = log(g1 + 0.00001f) - log(g2 + 0.00001f);
    Bdiff = log(b1 + 0.00001f) - log(b2 + 0.00001f);
  } else {
    printf("ERROR, No valid metric is found\n");
    return 0.f;
  }

  // The different way to compute the difference
  float diff = 0.f;
  if (metric == ETVI) {
    diff = luminance(fabs(Rdiff),fabs(Gdiff),fabs(Bdiff));  //< Luminance computation
    diff /= tvi(luminance(r2, g2, b2));
  } else if(metric == ERelMSE) {
    float refGray = (r2*r2+g2*g2+b2*b2);
    diff = (Rdiff*Rdiff+Gdiff*Gdiff+Bdiff*Bdiff);
    diff *= diff;
    diff /= (refGray + 0.0001);
  } else if (metric == ERelative) {
    diff = luminance(fabs(Rdiff),fabs(Gdiff),fabs(Bdiff));
    diff /= (luminance(r2, g2, b2) + 0.0001);
  } else if(metric == EMSE || metric == ERMSE ||
    metric == EMSE_LOG || metric == ERMSE_LOG){
    diff = Rdiff * Rdiff + Gdiff * Gdiff + Bdiff * Bdiff;
  } else if (metric == ESMAPE) {
    diff = fabs(Rdiff) + fabs(Gdiff) + fabs(Bdiff);
    diff /= (r1 + r2 + g1 + g2 + b1 + b2 + 0.0001);
    diff *= 2.0;
  } else {
    printf("ERROR, No valid metric is found (Diff computation)\n");
    return 0.f;
  }
  return diff;
}

static float errorNorm(float error,
                       int n, // Number of pixels treated
                       EErrorMetric metric) {
  if (n != 0) {
    error /= (float) n;
  }
  if (metric == ERMSE || metric == ERMSE_LOG) {
    error = sqrt(error);
  }
  return error;
}

float metric(imgRGB& img1, imgRGB& img2,
			 std::vector<unsigned char>& imgdiff,
             const std::vector<bool>& mask = std::vector<bool>(),
             EErrorMetric metric = EMSE) {
#if VERBOSE
  printf("Computing mse with a mult of %f\n",mult);
#endif
  int n = 0.f;
  float error = 0.f;
  float divFactor = 0.f;

  for (std::size_t i = 0; i < img1.size(); ++i) {
    // If we have a mask and the current mask pixel is non-zero
    if (mask.size() != 0 && mask[i]) {
      // We need to skip that pixel
      if (imgdiff.size() != 0) {
        imgdiff[i * 3] = 0;
        imgdiff[i * 3 + 1] = 0;
        imgdiff[i * 3 + 2] = 0;
      }
      continue;
    }

    float diff = metricPix(img1[i], img2[i], metric);
    if (diff < 0) {
      printf("Error, diff negative !!!!\n ");
    }

    // Compute the total error
    error += diff;
    ++n;

    // We compute the fake color
    if (imgdiff.size() != 0) {
      float distance = diff;
      if (metric == ERMSE || metric == ERMSE_LOG) {
        distance = sqrt(distance);
      }
      falseColor(distance, 0.f, 1.f, imgdiff[i * 3], imgdiff[i * 3 + 1],
                 imgdiff[i * 3 + 2]);
    }

  }

  error = errorNorm(error, n, metric);

  return error;
}

#endif // IMAGEERRORS_H
