#pragma once

#include "format.h"

#include <ImfInputFile.h>
#include <ImfChannelList.h>
#include <ImfFrameBuffer.h>
#include <ImfRgbaFile.h>
#include <half.h>

using namespace Imf;
using namespace Imath;

class FormatOpenEXR : public Format {
public:
  FormatOpenEXR(int w, int h, const std::vector<std::vector<float>> &d) :
      Format(w, h, d) {
  }
  FormatOpenEXR(const std::string &path) :
      Format(path) {
    read(path);
  }

  virtual void write(const std::string &path) {
    Rgba *hrgba = new Rgba[width * height];
    for (int i = 0; i < width * height; ++i)
      hrgba[i] = Rgba(data[3*i], data[3*i+1], data[3*i+2], 1.0);

    Box2i displayWindow(V2i(0,0), V2i(width-1, height-1));
    Box2i dataWindow = displayWindow;

    RgbaOutputFile file(path.c_str(), displayWindow, dataWindow, WRITE_RGB);
    file.setFrameBuffer(hrgba, 1, width);
    try {
      file.writePixels(height);
    }
    catch (const std::exception &e) {
      fprintf(stderr, "Unable to write image file \"%s\": %s", path,
              e.what());
    }

    delete[] hrgba;
  }

  virtual void read(const std::string &path) {
    try {
      InputFile file(path.c_str());
      Box2i dw = file.header().dataWindow();
      width = dw.max.x - dw.min.x + 1;
      height = dw.max.y - dw.min.y + 1;

      half *hrgba = new half[3 * height * width];
      int nChannels = 3;

      hrgba -= 3 * (dw.min.x + dw.min.y * width);
      FrameBuffer frameBuffer;
      frameBuffer.insert("R", Slice(HALF, (char *) hrgba,
                                    3 * sizeof(half), width * 3 * sizeof(half), 1, 1, 0.0));
      frameBuffer.insert("G", Slice(HALF, (char *) hrgba + sizeof(half),
                                    3 * sizeof(half), width * 3 * sizeof(half), 1, 1, 0.0));
      frameBuffer.insert("B", Slice(HALF, (char *) hrgba + 2 * sizeof(half),
                                    3 * sizeof(half), width * 3 * sizeof(half), 1, 1, 0.0));

      file.setFrameBuffer(frameBuffer);
      file.readPixels(dw.min.y, dw.max.y);

      hrgba += 3 * (dw.min.x + dw.min.y * width);
      data = (float *) malloc(sizeof(float) * 3 * width * height);
      for (int i = 0; i < nChannels * height * width; ++i)
        data[i] = hrgba[i];
      delete[] hrgba;
    } catch (const std::exception &e) {
      fprintf(stderr, "Unable to read image file \"%s\": %s", path, e.what());
    }
  }
};