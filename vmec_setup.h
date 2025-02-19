#pragma once
#include "vmec_defines.h"
#include "vmec_types.h"

#include <random>


void generate_noise_points(const int octaves, const int sizeFirstOctave, int &num_disk_points, noise_point* &disk_noise_points);
void generate_lattice_noise_points(const int octaves, const int sizeFirstOctave, int &num_lattice_points, lattice_noise_point* &lattice_noise_points);
void do_setup(pix_realNumber* &pixels_raw, pix_RGB* &pixels_clean, int &num_lattice_points, lattice_noise_point* &lattice_noise_points, int &num_disk_points, noise_point* &disk_noise_points);
void do_cleanup(pix_realNumber* &pixels_raw, pix_RGB* &pixels_clean, lattice_noise_point* &lattice_noise_points, noise_point* &disk_noise_points);
