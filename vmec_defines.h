#pragma once
#define _USE_MATH_DEFINES


#include <math.h>
#include <iostream>
#include <cstdlib>
#include <algorithm>

#include "vmec_code_helpers.hpp"

#ifndef VMEC_DEFINES_H
#define VMEC_DEFINES_H

//Set floating point precision
#ifdef USE_SINGLE_PRECISION_FLOATS
typedef float realNumber;
#else
typedef double realNumber;
#endif

realNumber operator ""_real(long double num);

struct vmec_settings
{
    // -- Rendering --
    //Render target
    u32 integration_depth = 156;
    u32 supersamples = 1;
    u32 image_height = 1000;
    u32 image_width = 2000;

    //animation
    bool animate_free_fall = true;
    bool animate_path = false;
    u32 frame_rate = 24;
    u32 number_of_frames = 196;
    u32 f_start = 1; //animation start frame;
    u32 f_end = 196; //animation end frame;
    u32 frame_step = 1; //# frames between each frame;
    realNumber time_ramp = 170847.0_real;//scaling factor for time;

    //scene
    realNumber celestial_sphere_radius = 400.0_real;

    //render mode
    bool render_debug = true;
    bool render_magik = false;
    bool render_magik_sandbox = false;

    // -- debugging --
    //checker disk
    bool checker_disk_toggle = true;
    realNumber checker_disk_major_radius = 56.0_real;
    realNumber checker_disk_minor_radius = 5.0_real;

    //lattice point visualization
    realNumber jet_lattice_point_size = 1.0_real;
    realNumber disk_lattice_point_size = 5.0_real;

    // -- volumetrics --
    //quality disk
    u32 disk_noise_octaves = 1;
    u32 disk_noise_size_first_octave = 24;
    u32 disk_n_points = disk_noise_size_first_octave * std::pow(2,disk_noise_octaves-1);

    //disk lattice points
    realNumber disk_lattice_major_radius = 90.0_real;
    realNumber disk_lattice_minimum_radius = 3.0_real;
    realNumber disk_lattice_vertical_scale = 0.05_real;
    realNumber disk_lattice_velocity_scale = 1.0_real; //MUST BE >= 1

    //bounding volume disk
    bool bv_disk_toggle = true;
    realNumber bv_disk_major_radius = 120.0_real;
    realNumber bv_disk_minor_radius = 0.0_real;
    realNumber bv_disk_mean_slope = 12.0_real;
    realNumber bv_disk_radius_offset = 24.0_real;

    //quality jet
    u32 jet_noise_octaves = 2;
    u32 jet_noise_first_octave = 12;
    u32 jet_n_points = jet_noise_first_octave * std::pow(2,jet_noise_octaves-1);

    //bounding volume jet
    bool bv_jet_toggle = true;
    realNumber bv_jet_height = 10.0_real;
    realNumber bv_jet_deflection = 22.0_real;
    realNumber bv_jet_minor_radius = 16.0_real;

    //size jet
    realNumber jet_minor_radius = 5.0_real;
    realNumber jet_radius_scale = 10.0_real;

    //bounding volume ambient medium
    bool bv_ambient_medium_toggle = true;
    realNumber bv_ambient_medium_radius = 120.0_real;

    // -- Physics Related --
    //Black Hole
    realNumber a = -0.999_real; //Spin parameter
    realNumber rs_scale = 1.01;

    //Astrophysical Jet
    realNumber vy = 6.0_real;//Jet Vertical Velocity 40.5
    realNumber w = 196.0_real; //Jet phi velocity scale 196
    realNumber vr = 24.5_real; //Jet "Spread Factor 212
    realNumber max_t = 8.0_real; //Max life time for jet particle 2
    realNumber alpha = 2.0_real; //Jet phi velocity decay power. Keep it at 2

};

extern vmec_settings settings;

#endif // VMEC_DEFINES_H
