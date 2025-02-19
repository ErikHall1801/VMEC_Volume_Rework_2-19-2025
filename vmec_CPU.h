#pragma once
#include "vmec_defines.h"
#include "vmec_common.h"
#include "vmec_magik.h"

#include <mutex>
#include <thread>
#include <random>


void trace_all_rays(const realNumber globalTime, const int total_threads, const SensorProperties sensorproperties, noise_point* disk_noise_points, const int num_lattice_points, lattice_noise_point* lattice_noise_points, pix_realNumber* pixels_raw);