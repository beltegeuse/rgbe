#pragma once

#include <vector>
#include <tuple>

#include <pybind11/stl.h>

#include <Python.h>

enum class Extension: char {
  Unknown = 0,
  RGBE = 1,
  PFM = 2,
  OpenEXR = 3,
};

void upper_string(char *string) {
  while (*string) {
    if (*string >= 'a' && *string <= 'z') {
      *string = *string - 32;
    }
    string++;
  }
}

Extension helper_get_ext(const char* path) {
  char * pathTmp;
  char * posPoint;

  pathTmp = strdup(path);
  upper_string(pathTmp);
  posPoint = strrchr(pathTmp, '.');

  if (posPoint == NULL) {
    return Extension::Unknown;
  }

  if (strcmp(posPoint, ".PFM") == 0) {
    return Extension::PFM;
  } else if (strcmp(posPoint, ".HDR") == 0) {
    return Extension::RGBE;
  } else if (strcmp(posPoint, ".RGBE") == 0) {
    return Extension::RGBE;
  } else if (strcmp(posPoint, ".EXR") == 0) {
    return Extension::OpenEXR;
  }

  free(pathTmp);

  return Extension::Unknown;
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
#include "openexr.h"

Format* loadImage(const std::string& path) {
  auto ext = helper_get_ext(path.c_str());
  if (ext == Extension::Unknown) {
    printf("ERROR: Unkown extension, QUIT: %s\n", path.c_str());
    return NULL;
  }

  Format* f = NULL;
  // Header reading
  if (ext == Extension::RGBE) {
    f = new FormatRGBE(path);
  } else if (ext == Extension::PFM) {
    f = new FormatPFM(path);
  } else if (ext == Extension::OpenEXR) {
    f = new FormatOpenEXR(path);
  } else {
    printf("ERROR: No Reader, QUIT: %s\n", path.c_str());
  }
  return f;
}


