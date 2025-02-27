#pragma once
#include "vmec_defines.h"
#include "vmec_types.h"

#include <fstream>


void image_generate(pix_realNumber* pixels_in, pix_RGB* pixels_out);
void write_binary_RBG(pix_RGB* pixels_all, int frame_number);
void write_binary_realNumber(pix_realNumber* pixels_all, int frame_number);
void write_hdr_realNumber(pix_realNumber* pixels_all, int frame_number);
void write_exr_realNumber(pix_realNumber* pixels_all, int frame_number);
