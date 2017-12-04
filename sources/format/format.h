#pragma once

#include <vector>
#include <tuple>

#include <pybind11/stl.h>

#include <Python.h>

#define EXT_UNKNOW  -1
#define EXT_RGBE  0
#define EXT_PFM   1

void upper_string(char *string) {
  while (*string) {
    if (*string >= 'a' && *string <= 'z') {
      *string = *string - 32;
    }
    string++;
  }
}

int helper_get_ext(const char* path) {
  char * pathTmp;
  char * posPoint;

  pathTmp = strdup(path);
  upper_string(pathTmp);
  posPoint = strrchr(pathTmp, '.');

  if (posPoint == NULL) {
    return EXT_UNKNOW;
  }

  if (strcmp(posPoint, ".PFM") == 0) {
    return EXT_PFM;
  } else if (strcmp(posPoint, ".HDR") == 0) {
    return EXT_RGBE;
  } else if (strcmp(posPoint, ".RGBE") == 0) {
    return EXT_RGBE;
  }

  free(pathTmp);

  return EXT_UNKNOW;
}

typedef std::tuple<float, float, float> RGB;
typedef std::vector<RGB> imgRGB;
std::vector<unsigned char> imgMask;
typedef std::tuple<std::tuple<int, int>, imgRGB> readReturn;

class Format {
protected:
  int width;
  int height;
  float* data;

public:
  Format(int w, int h, const std::vector<std::vector<float>>& d):
    width(w), height(h)
  {
      data = (float *) malloc(sizeof(float) * 3 *width*height);
      for (unsigned int i = 0; i < width*height; i++) {
          for (unsigned int j = 0; j < 3; j++) {
              // XXX Protect the reading
              data[i * 3 + j] = d[i][j];
          }
      }
  }
  Format(const std::string& path):
    width(0), height(0), data(NULL)
  {
  }

  virtual ~Format() {
    if(data) {
      free(data);
    }
  }

  virtual void write(const std::string& path) = 0;
  virtual void read(const std::string& path) = 0;

  int getWidth() const {
    return width;
  }
  int getHeight() const {
    return height;
  }
  float* getData() const {
    return data;
  }

  // returns ((width, height), listImage)
  readReturn pack() const {
    std::vector<std::tuple<float, float, float>> pixels(width * height);
    for (int i = 0; i < width * height; i++) {
      pixels[i] = std::make_tuple( data[i * 3],  data[i * 3 + 1],  data[i * 3 + 2]);
    }
    return std::make_tuple(std::make_tuple(width, height),
                           pixels);
  }

  double* toDouble(float mult=1.f) const {
    double * v = new double[width*height*3];
    for(int i = 0; i < width*height*3; i++) {
      v[i] = data[i]*mult;
    }
    return v;
  }
};

#include "rgbe.h"
#include "pfm.h"

Format* loadImage(const std::string& path) {
  int ext = helper_get_ext(path.c_str());
  if (ext == EXT_UNKNOW) {
    printf("ERROR: Unkown extension, QUIT: %s\n", path.c_str());
    return NULL;
  }

  Format* f = NULL;
  // Header reading
  if (ext == EXT_RGBE) {
    f = new FormatRGBE(path);
  } else if (ext == EXT_PFM) {
    f = new FormatPFM(path);
  } else {
    printf("ERROR: No Reader, QUIT: %s\n", path.c_str());
  }
  return f;
}


