#pragma once
#include "vmec_defines.h"
#include "vmec_common.h"
#include "vmec_magik.h"

#include "pch_vmec_primary.h"


void trace_all_rays(const realNumber globalTime, const int total_threads, const SensorProperties sensorproperties, RayProperties* all_geodesics, int* geodesic_sizes, realNumber* geodesics_k, noise_point *disk_lattice_noise_points, const int num_lattice_points, lattice_noise_point *FRAME_jet_lattice_noise_points, dynamic_lattice_point *dynamic_noise_points, pix_realNumber *pixels_raw);