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

std::vector<unsigned char> rgbe_rmse(int width, int height,
              imgRGB& imgHDR1,
              imgRGB& imgHDR2,
              int typeMetric = 0,
              const std::vector<bool>& mask = std::vector<bool>()) {
	std::vector<unsigned char> imgDiff(imgHDR1.size() * 3);

	if (typeMetric < 0 || typeMetric > ENUMBER_METRIC) {
      printf("Error when reading the arguments, aborting.\n");
      return imgDiff;
    }

    float error = metric(imgHDR1, imgHDR2, imgDiff, mask, (EErrorMetric) typeMetric);
	std::cout << "error: " << error << "\n";

    // We return the rmse and the difference image
    return imgDiff;
}

std::vector<float> computeErrors(imgRGB& imageHDR, 
				 imgRGB& imgHDRRef, 
				 const std::vector<bool>& mask = std::vector<bool>()) {
    std::vector<float> errors(ENUMBER_METRIC, 0.f);
    std::vector<unsigned char> imgDiff;
	for (int idError = 0; idError < ENUMBER_METRIC; idError++) {
        // We compute the RMSE
        errors[idError] = metric(imageHDR, imgHDRRef, imgDiff, mask, (EErrorMetric) idError);
    }
    return errors;
}

std::vector<float> rgbe_rmse_all(int width, int height, imgRGB& imgHDR1,
                                 imgRGB& imgHDR2, float mult = 1.f,
                                 const std::vector<bool>& mask = std::vector<bool>()) {
#if VERBOSE
    printf("Reading of the parameters successful, mult is %f\n", mult);
#endif
    std::vector<float> errors(ENUMBER_METRIC, -1.f);

    // Read all images
    applyScale(imgHDR1, mult);
    applyScale(imgHDR2, mult);

    errors = computeErrors(imgHDR1, imgHDR2, mask);

    // We return the rmse and the difference image
    return errors;
}

std::vector<std::vector<float>> rgbe_rmse_all_images(int width, int height,
                                                     const std::vector<std::string>& images,
                                                     imgRGB& imgHDRRef, float mult = 1.f,
                                                     const std::vector<bool>& mask = std::vector<bool>()) {
    std::vector<std::vector<float> > metrics(images.size());

    // Read the reference image and mask
    applyScale(imgHDRRef, mult);

#pragma omp parallel for
    for (int i = 0; i < (int) images.size(); i++) {
        // Open the image
        Format* f = loadImage(images[i]);
        std::vector<float> errors(ENUMBER_METRIC, -1.f);

        if (f != NULL && f->getWidth() == width && f->getHeight() == height) {
            // FIXME: avoid to do this conversion
            imgRGB imageHDR = std::get<1>(f->pack());
            applyScale(imageHDR, mult);
            // Else, compute the metric
            errors = computeErrors(imageHDR, imgHDRRef, mask);
        }

        metrics[i] = errors;

        if(f != NULL) delete f;

        // Display computation results
        std::cout << images[i] << ": ( ";
        for(auto m : errors) {
            std::cout << m << "; ";
        }
        std::cout << ")\n";
    }

    return metrics;
}

std::vector<std::vector<float>> rgbe_rmse_all_images_percentage(int width, int height, float percentage,
                                                                const std::vector<std::string>& images,
                                                                imgRGB& imgHDRRef, float mult = 1.f) {
    std::vector<std::vector<float> > metrics(images.size());

    if(percentage > 1.f || percentage <= 0.f) {
        std::cerr << "Invalid percentage: " << percentage << "\n";
        return metrics; // Empty all
        // FIXME: Add the -1 for all
    }

    // Read the reference image and mask
    applyScale(imgHDRRef, mult);

#pragma omp parallel for
    for (int i = 0; i < (int) images.size(); i++) {
        // Open the image
        Format* f = loadImage(images[i]);
        std::vector<float> errors(ENUMBER_METRIC, -1.f);

        if (f != NULL && f->getWidth() == width && f->getHeight() == height) {
            // FIXME: avoid to do this conversion
            imgRGB imageHDR = std::get<1>(f->pack());
            applyScale(imageHDR, mult);
            // Else, compute the metric
            for (int idError = 0; idError < ENUMBER_METRIC; idError++) {
                // Compute error pixel wise
                std::vector<float> pixelErrors(f->getWidth()*f->getHeight(), 0.0f);
                for (int i = 0; i < width * height; ++i) {
                    pixelErrors[i] = metricPix(imageHDR[i], imgHDRRef[i], (EErrorMetric) idError);
                }

                // Sort the metric vector
                // and compute the metric for the percentage of selected pixels
                std::sort(pixelErrors.begin(), pixelErrors.end());
                float error = 0.f;
                for (int i = 0; i < width * height * percentage; ++i) {
                    error += pixelErrors[i];
                }
                errors[idError] = errorNorm(error, width * height * percentage, (EErrorMetric) idError);
            }
        }

        metrics[i] = errors;

        if(f != NULL) delete f;

        // Display computation results
        std::cout << images[i] << ": ( ";
        for(auto m : errors) {
            std::cout << m << "; ";
        }
        std::cout << ")\n";
    }

    return metrics;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Write Techniques
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/// write HDR files
// Automatically choose the format
bool rgbe_write(const std::string& path, int width, int height,
                const std::vector<std::vector<float>>& image) {
    // Read extension
    auto ext = helper_get_ext(path.c_str());
    if (ext == Extension::Unknown) {
        printf("ERROR: Unkown extension, QUIT: %s\n", path.c_str());
        return false;
    }

    // Write the format by calling the dedicated function
    if (ext == Extension::RGBE) {
        FormatRGBE f(width, height, image);
        f.write(path);
        return true;
    } else if (ext == Extension::PFM) {
        FormatPFM f(width, height, image);
        f.write(path);
        return true;
    } else if (ext == Extension::OpenEXR) {
        FormatOpenEXR f(width, height, image);
        f.write(path);
        return true;
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
readReturn rgbe_read(const std::string& path) {
    auto ext = helper_get_ext(path.c_str());
    if (ext == Extension::Unknown) {
        printf("ERROR: Unkown extension, QUIT: %s\n", path.c_str());
        std::vector<std::tuple<float,float,float>> pixels(0);
        return std::make_tuple(std::make_tuple(0,0), pixels);
    }

    // === Read info
    if (ext == Extension::RGBE) {
        FormatRGBE f(path);
        return f.pack();
    } else if (ext == Extension::PFM) {
        FormatPFM f(path);
        return f.pack();
    } else if (ext == Extension::OpenEXR) {
        FormatOpenEXR f(path);
        return f.pack();
    } else {
        printf("ERROR: No Reader, QUIT: %s\n", path.c_str());
        std::vector<std::tuple<float,float,float>> pixels(0);
        return std::make_tuple(std::make_tuple(0,0), pixels);
    }
}

std::tuple < std::tuple<int, int>, std::vector<std::tuple<int, int, int>> > rgbe_read_tonemap(const std::string &path, float gamma = 2.2f, float exposure = 0.0f)
{
    Format *f = loadImage(path);
    // Header reading
    if (f == NULL)
    {
        std::vector<std::tuple<int, int, int>> pixels(0);
        return std::make_tuple(std::make_tuple(0, 0), pixels);
    }

    // List creation
    std::vector<std::tuple<int, int, int>> pixels(f->getWidth() * f->getHeight());
    exposure = powf(2, exposure);
    float *image = f->getData();
    for (int i = 0; i < f->getWidth() * f->getHeight(); i++)
    {
        pixels[i] = std::make_tuple(
            (int)(powf(image[i * 3] * exposure, 1.f / gamma) * 255),
            (int)(powf(image[i * 3 + 1] * exposure, 1.f / gamma) * 255),
            (int)(powf(image[i * 3 + 2] * exposure, 1.f / gamma) * 255));
    }

    return std::make_tuple(std::make_tuple(f->getWidth(), f->getHeight()),
                           pixels);
}

void rgbe_merge(const std::string& path, int w, int h, int blockSize,
                const std::string& outpath, const std::string& ext = "hdr") {
    // read local ...
    // FIXME: Need to build the list of all the images
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

    // Special read methods
    mIO.def("read_tonemap", &rgbe_read_tonemap, "read and tonemap an hdr image");

    // Write methods
    mIO.def("write", &rgbe_write, "write hdr file (.pfm, .hdr)");

    pybind11::module mMerge = m.def_submodule("merge", "A submodule of 'example'");
    mMerge.doc() = "merge multiple hdr files into a big hdr file";
    mMerge.def("merge", &rgbe_merge, "merge several small hdr image into a big one");

    pybind11::module mUtils = m.def_submodule("utils", "A submodule of 'example'");
    mUtils.doc() = "utility package for hdr common manipulation (i.e. tonemap)";
    mUtils.def("applyExposureGamma", &rgbe_applyExposureGamma, "apply exposure and gamma correction",
               "pixels"_a, "exposure"_a = 1.0, "gamma"_a = 2.2);
}
