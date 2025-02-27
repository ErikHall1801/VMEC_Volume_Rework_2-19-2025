#pragma once
#define _USE_MATH_DEFINES

#ifndef VMEC_DEFINES_H
#define VMEC_DEFINES_H

#include "pch_vmec_primary.h"

//Set floating point precision
#ifdef USE_SINGLE_PRECISION_FLOATS
typedef float realNumber;
#else
typedef double realNumber;
#endif

using std::vector;

realNumber operator ""_real(long double num);

struct global_database
{
    vector<fstring> args;
};

struct vmec_settings
{
    // -- General Stuff --
    bool arg_help = false;
    bool arg_skip_render = false;
    bool arg_verbose = false;
    u64 session_timestamp = 0; 
    fstring session_dir;
    fstring render_dir = "Render";
    fstring render_binary_dir = "bin";
    fstring render_exr_dir = "EXR";
    fstring render_png_dir = "PNG";

    // -- Rendering --
    //Render target 
    u32 integration_depth = 156;
    u32 supersamples = 1;
    u32 image_height = 250;
    u32 image_width = 500; 

    //animation
    bool animate_free_fall = true;
    bool animate_path = false;
    u32 frame_rate = 24;
    u32 number_of_frames = 1;
    u32 f_start = 0; //animation start frame;
    u32 f_end = 1; //animation end frame;
    u32 frame_step = 1; //# frames between each frame;
    realNumber time_ramp = 2000.0_real;//scaling factor for time;

    //scene
    realNumber celestial_sphere_radius = 400.0_real;

    //render mode
    bool render_debug = true;
    bool render_magik = false;
    bool render_magik_sandbox = true;

    // -- magik sandbox --
    //trivial
    bool magik_trivial_scene = false;

    //single
    bool magik_single_scene = false;
    bool magik_single_MFP_scene = true;

    //multiple
    bool magik_muliple_scene = true;

    // -- debugging --
    //test scenes
    bool test_scene = false;
    bool lattice_scene = true;

    //common
    bool colored_background = false;
    bool integration_depth_overlay = false;
    bool coordinate_overlay = false;
    
    //checker disk
    bool checker_disk_toggle = true; 
    bool checker_disk_Bi_sexual_color = true;
    bool checker_disk_Z_factor = false;
    bool checker_disk_velocity_gradient = false;
    realNumber checker_disk_major_radius = 56.0_real;
    realNumber checker_disk_minor_radius = 5.0_real;

    //lattice point visualization
    bool toggle_disk_lattice_points = false; 
    bool toggle_jet_lattice_points = true;
    realNumber jet_lattice_point_size = 2.0_real;
    realNumber disk_lattice_point_size = 4.0_real;

    // -- volumetrics --
    //quality disk
    u32 disk_noise_octaves = 1;
    u32 disk_noise_size_first_octave = 24;
    u32 disk_n_points = disk_noise_size_first_octave * std::pow(2,disk_noise_octaves-1);

    //disk lattice points
    realNumber disk_lattice_major_radius = 90.0_real;
    realNumber disk_lattice_minimum_radius = 3.0_real;
    realNumber disk_lattice_vertical_scale = 0.05_real;
    realNumber disk_lattice_velocity_scale = 5.0_real; //MUST BE >= 1

    //bounding volume disk
    bool bv_disk_toggle = false;
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
    realNumber bv_jet_height = 12.0_real;
    realNumber bv_jet_deflection = 22.0_real;
    realNumber bv_jet_minor_radius = 24.0_real;

    //size jet
    realNumber jet_minor_radius = 5.0_real;
    realNumber jet_radius_scale = 10.0_real;

    //bounding volume ambient medium
    bool bv_ambient_medium_toggle = false;
    realNumber bv_ambient_medium_radius = 120.0_real;

    // -- Physics Related --
    //Black Hole
    realNumber a = -0.999_real; //Spin parameter
    realNumber rs_scale = 1.01;
    realNumber rs = (1.0+std::sqrt(1.0-(std::abs(a)*std::abs(a))));

    //Astrophysical Jet
    realNumber vy = 6.0_real;//Jet Vertical Velocity 40.5
    realNumber w = 196.0_real; //Jet phi velocity scale 196
    realNumber vr = 24.5_real; //Jet "Spread Factor 212
    realNumber max_t = 8.0_real; //Max life time for jet particle 2
    realNumber alpha = 2.0_real; //Jet phi velocity decay power. Keep it at 2

};

extern global_database global;
extern vmec_settings settings;

#endif // VMEC_DEFINES_H
