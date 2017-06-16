#include <string>

#include "config.h"
#include "imageerrors.h"
#include "format/format.h"

// Python warpper
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Other include
#include <Python.h>
#include <omp.h>

// STL include
#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>

void applyScale(imgRGB& img, float mult = 1.f) {
    for(std::size_t i = 0; i < img.size(); i++) {
        img[i] = std::make_tuple(std::get<0>(img[i])*mult,
                                 std::get<1>(img[i])*mult,
                                 std::get<2>(img[i])*mult);
    }
}

int rgbe_rmse(int width, int height,
              imgRGB& imgHDR1,
              imgRGB& imgHDR2,
              float mult = 1.f,
              int typeMetric = 0,
              PyObject *imgMask = NULL) {
    std::cerr << "not implemented...\n";
    return 0;

//  // Check parameters
//  if (typeMetric < 0 || typeMetric > 4) {
//    printf("Error when reading the arguments, aborting.\n");
//    return 0;
//  }
//#if VERBOSE
//  printf("Reading of the parameters successful, mult is %f\n", mult);
//#endif
//
//  // Read all images
//    applyScale(imgHDR1, mult);
//    applyScale(imgHDR2, mult);
//    unsigned char * imgMaskData = 0;
//    if (imgMask != 0 && imgMask != Py_None) {
//        printf("No support for masking for now");
//        return 0;
//     }
//
//  // Compute the difference and fill diff image
//  unsigned char * imgDiffData = new unsigned char[width * height * 3];
//  float error = metric(imgHDR1Data, imgHDR2Data, imgMaskData, imgDiffData,
//      width, height, (EErrorMetric) typeMetric);
//
//  // We write back the difference image to a python object
//  // And pack the results
//  PyObject * imgDiff = PyList_New(width * height);
//  for (int i = 0; i < width * height; ++i) {
//    PyList_SetItem(imgDiff, i,
//        Py_BuildValue("(i,i,i)", imgDiffData[i * 3], imgDiffData[i * 3 + 1],
//            imgDiffData[i * 3 + 2]));
//  }
//  PyObject * resultsPy = Py_BuildValue("fO", error, imgDiff);
//  Py_DECREF(imgDiff);
//
//  // Free memory
//  delete[] imgDiffData;
//  if (imgMaskData) {
//    delete[] imgMaskData;
//  }
//
//  // We return the rmse and the difference image
//  return 0;
}

std::vector<float> rgbe_rmse_all(int width, int height, imgRGB& imgHDR1,
                                 imgRGB& imgHDR2, float mult = 1.f,
                                 py::object imgMask = py::none()) {
#if VERBOSE
    printf("Reading of the parameters successful, mult is %f\n", mult);
#endif
    std::vector<float> errors(6, 0.f);

    // Read all images
    applyScale(imgHDR1, mult);
    applyScale(imgHDR2, mult);
    unsigned char * imgMaskData = 0;
    if (!imgMask.is_none()) {
        std::cerr << "No support for masking for now\n";
        return errors;
    }

    for (int idError = 0; idError < 6; idError++) {
        // We compute the RMSE
        float error = metric(imgHDR1, imgHDR2, imgMaskData, NULL, (EErrorMetric) idError);
        errors[idError] = error;
    }

    // We return the rmse and the difference image
    return errors;
}

std::vector<std::vector<float>> rgbe_rmse_all_images(int width, int height,
                                                     const std::vector<std::string>& images,
                                                     imgRGB& imgHDRRef, float mult = 1.f,
                                                     py::object imgMask = py::none()) {
    std::vector<std::vector<float> > metrics(images.size());

    // Read the reference image and mask
    applyScale(imgHDRRef, mult);
    unsigned char * imgMaskData = 0;
    if (!imgMask.is_none()) {
//        imgMaskData = convertListIntoArrayMask(imgMask, width * height);
        std::cerr << *imgMask;
        std::cerr << "No support for masking for now\n";
        return metrics;
    }

    const int NBMETRIC = 7;

#pragma omp parallel for
    for (int i = 0; i < (int) images.size(); i++) {
        // Open the image
        Format* f = loadImage(images[i]);
        if (f != NULL && f->getWidth() == width && f->getHeight() == height) {
            // FIXME: avoid to do this conversion
            imgRGB imageHDR = std::get<1>(f->pack());
            applyScale(imageHDR, mult);
            // Else, compute the metric
            for (int idError = 0; idError < NBMETRIC; idError++) {
                // We compute the RMSE
                float error = metric(imageHDR, imgHDRRef, imgMaskData, NULL, (EErrorMetric) idError);
                metrics[i].push_back(error);
            }
        } else {
            // If there is an error on the image reading
            for (int i = 0; i < NBMETRIC; i++) {
                metrics[i].push_back(-1.f); // Invalid number
            }
        }

        delete f;
        std::cout << images[i] << ": ( " << metrics[i][0] << "; " << metrics[i][1]
                  << "; " << metrics[i][2] << "; " << metrics[i][3] << "; "
                  << metrics[i][4] << "; " << metrics[i][5] << "; " << metrics[i][6] << ")\n";
    }

    // Free memory
    if (imgMaskData)
        delete[] imgMaskData;

    return metrics;
}

std::vector<std::vector<float>> rgbe_rmse_all_images_percentage(int width, int height, float percentage,
                                                                const std::vector<std::string>& images,
                                                                imgRGB& imgHDR2, float mult = 1.f) {
    std::vector<std::vector<float> > metrics(images.size());
    std::cerr << "not implemented\n";
    return metrics;

//    if(percentage > 1.f || percentage <= 0.f) {
//        std::cerr << "Invalid percentage: " << percentage << "\n";
//        return 0;
//    }
//
//    // Read the reference image and mask
//    double * imgHDRRef = convertListIntoArray(imgHDR2, width * height, mult);
//    unsigned char * imgMaskData = 0; // No mask for now
//
//    const int NBMETRIC = 7;
//
//    // Launch the computation
//    std::vector<std::vector<float> > metrics(images.size());
//
//#pragma omp parallel for
//    for (int i = 0; i < (int) images.size(); i++) {
//        // Open the image
//        Format* f = loadImage(images[i]);
//        if (f != NULL && f->getWidth() == width && f->getHeight() == height) {
//            double * imageHDRDouble = f->toDouble(mult);
//            // Else, compute the metric
//            for (int idError = 0; idError < NBMETRIC; idError++) {
//                // Compute error pixel wise
//                std::vector<float> pixelErrors(f->getWidth()*f->getHeight(), 0.0f);
//                for (int i = 0; i < width * height; ++i) {
//                    pixelErrors[i] = metricPix(imageHDRDouble, imgHDRRef, i, (EErrorMetric) idError);
//                }
//
//                // Sort the metric vector
//                // and compute the metric for the percentage of selected pixels
//                std::sort(pixelErrors.begin(), pixelErrors.end());
//                float error = 0.f;
//                for (int i = 0; i < width * height * percentage; ++i) {
//                    error += pixelErrors[i];
//                }
//                error = errorNorm(error, width * height * percentage, (EErrorMetric) idError);
//
//                // We compute the RMSE
//                metrics[i].push_back(error);
//            }
//            delete[] imageHDRDouble;
//        } else {
//            // If there is an error on the image reading
//            for (int i = 0; i < NBMETRIC; i++) {
//                metrics[i].push_back(-1.f); // Invalid number
//            }
//        }
//
//        delete f;
//        std::cout << images[i] << ": ( " << metrics[i][0] << "; " << metrics[i][1]
//                  << "; " << metrics[i][2] << "; " << metrics[i][3] << "; "
//                  << metrics[i][4] << "; " << metrics[i][5] << "; " << metrics[i][6] << ")\n";
//    }
//
//    // Free memory
//    delete[] imgHDRRef;
//    if (imgMaskData)
//        delete[] imgMaskData;
//
//    // Construct the python results
//    PyObject * resPy = PyList_New(images.size());
//    for(int i = 0; i < images.size(); i++) {
//        PyList_SetItem(resPy, i, Py_BuildValue("fffffff", metrics[i][0], metrics[i][1],
//                                               metrics[i][2], metrics[i][3], metrics[i][4], metrics[i][5],metrics[i][6]));
//    }
//
//    return resPy;
}


float* convertImage(int width, int height, PyObject *imagePy) {
    float* image = (float *) malloc(sizeof(float) * 3 *width*height);
    for (unsigned int i = 0; i < width*height; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            // XXX Protect the reading
            image[i * 3 + j] = (float) PyFloat_AsDouble(
              PyTuple_GetItem(PyList_GetItem(imagePy, i), j));
        }
    }
    return image;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Write Techniques
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/// write Radiance files
bool rgbe_write_hdr(const std::string& path, int width, int height, PyObject *imagePy) {
    float* image = convertImage(width, height, imagePy);
    FormatRGBE f(width, height, image);
    f.write(path);

    return true;
}

/// write PFM files
bool rgbe_write_pfm(const std::string& path, int width, int height, PyObject *imagePy) {
    float* image = convertImage(width, height, imagePy);
    FormatPFM f(width, height, image);
    f.write(path);

    return true;
}

/// write HDR files
// Automatically choose the format
bool rgbe_write(const std::string& path, int width, int height, PyObject *imagePy) {
    // Read extension
    int ext = helper_get_ext(path.c_str());
    if (ext == EXT_UNKNOW) {
        printf("ERROR: Unkown extension, QUIT: %s\n", path.c_str());
        return false;
    }

    // Write the format by calling the dedicated function
    if (ext == EXT_RGBE) {
        return rgbe_write_hdr(path, width, height, imagePy);
    } else if (ext == EXT_PFM) {
        return rgbe_write_pfm(path, width, height, imagePy);
    } else {
        printf("ERROR: No writter, QUIT: %s\n", path.c_str());
        return false;
    }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Read Techniques
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/// read HDR files
readReturn rgbe_read_hdr(const std::string& path) {
    FormatRGBE f(path);
    return f.pack();
}

/// read HDR files
readReturn rgbe_read_pfm(const std::string& path) {
    FormatPFM f(path);
    return f.pack();
}

/// read HDR files
readReturn rgbe_read(const std::string& path) {
    int ext = helper_get_ext(path.c_str());
    if (ext == EXT_UNKNOW) {
        printf("ERROR: Unkown extension, QUIT: %s\n", path.c_str());
        std::vector<std::tuple<float,float,float>> pixels(0);
        return std::make_tuple(std::make_tuple(0,0), pixels);
    }

    // === Read info
    if (ext == EXT_RGBE) {
        return rgbe_read_hdr(path);
    } else if (ext == EXT_PFM) {
        return rgbe_read_pfm(path);
    } else {
        printf("ERROR: No Reader, QUIT: %s\n", path.c_str());
        std::vector<std::tuple<float,float,float>> pixels(0);
        return std::make_tuple(std::make_tuple(0,0), pixels);
    }
}

PyObject *rgbe_read_tonemap(const std::string& path, float gamma = 2.2f, float exposure = 0.0f) {
    Format* f = loadImage(path);
    // Header reading
    if (f == NULL) {
        return Py_BuildValue("(i,i)", 0, 0);
    }

    // List creation
    PyObject* pixels = PyList_New(f->getWidth() * f->getHeight());
    exposure = powf(2, exposure);
    float* image = f->getData();
    for (int i = 0; i < f->getWidth()*f->getHeight(); i++) {
        PyList_SetItem(pixels, i,
                       Py_BuildValue("(i,i,i)",
                                     (int) (powf(image[i * 3] * exposure, 1.f / gamma) * 255),
                                     (int) (powf(image[i * 3 + 1] * exposure, 1.f / gamma) * 255),
                                     (int) (powf(image[i * 3 + 2] * exposure, 1.f / gamma) * 255)));
    }

    PyObject* resultsPy = Py_BuildValue("((i,i),O)", f->getWidth(), f->getHeight(), pixels);
    Py_DECREF(pixels);

    return resultsPy;
}

void rgbe_merge(const std::string& path, int w, int h, int blockSize,
                const std::string& outpath, const std::string& ext = "hdr") {
    // read local ...
}

void rgbe_applyExposureGamma(pybind11::list img, double exposure = 1.0, double gamma = 2.2) {
    double exp = std::pow(2.0,exposure);
    double invG = 1.0/gamma;
    for(std::size_t i = 0; i < img.size(); i++) {
        pybind11::tuple pix = pybind11::tuple(img[i]);
        img[i] = std::make_tuple(std::pow(pybind11::float_(pix[0]).cast<float>()*exp, invG),
                                 std::pow(pybind11::float_(pix[1]).cast<float>()*exp, invG),
                                 std::pow(pybind11::float_(pix[2]).cast<float>()*exp, invG));
    }
}

PYBIND11_MODULE(rgbe, m) {

    using namespace pybind11::literals;

    pybind11::module mFast = m.def_submodule("fast", "A submodule of 'example'");
    mFast.doc() = "fast operation on hdr images (i.e metric computation)";

    // Read methods
    mFast.def("rmse", &rgbe_rmse, "compute a designed metric");
    mFast.def("rmse_all", &rgbe_rmse_all, "compute all metric provided by this module");
    mFast.def("rmse_all_images", &rgbe_rmse_all_images, "compute all metric for a list of images");
    mFast.def("rmse_all_images_percentage", &rgbe_rmse_all_images_percentage, "compute all metric with metric clamping");

    pybind11::module mIO = m.def_submodule("io", "A submodule of 'example'");
    mIO.doc() = "read and write HDR files (.pfm, .hdr)";
    mIO.def("read", &rgbe_read, "read hdr file (.pfm, .hdr)");
    mIO.def("read_pfm", &rgbe_read_pfm, "read pfm hdr file (.pfm)");
    mIO.def("read_hdr", &rgbe_read_hdr, "read rgbe hdr file (.hdr)");

    // Special read methods
    mIO.def("read_tonemap", &rgbe_read_tonemap, "read and tonemap an hdr image");

    // Write methods
    mIO.def("write", &rgbe_write, "write hdr file (.pfm, .hdr)");
    mIO.def("write_pfm", &rgbe_write_pfm, "write pfm hdr file (.pfm)");
    mIO.def("write_hdr", &rgbe_write_hdr, "write rgbe hdr file (.pfm)");

    pybind11::module mMerge = m.def_submodule("merge", "A submodule of 'example'");
    mMerge.doc() = "merge multiple hdr files into a big hdr file";
    mMerge.def("merge", &rgbe_merge, "merge several small hdr image into a big one");

    pybind11::module mUtils = m.def_submodule("utils", "A submodule of 'example'");
    mUtils.doc() = "utility package for hdr common manipulation (i.e. tonemap)";
    mUtils.def("applyExposureGamma", &rgbe_applyExposureGamma, "apply exposure and gamma correction",
               "pixels"_a, "exposure"_a = 1.0, "gamma"_a = 2.2);
}
