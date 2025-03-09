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
    u32 image_height = 1000;  
    u32 image_width = 1660;

    //animation
    bool static_camera = true;
    bool animate_free_fall = false;
    bool animate_path = false;
    u32 frame_rate = 24; 
    realNumber frame_time_step = 1.0/realNumber(frame_rate);
    u32 number_of_frames = 1;
    u32 f_start = 0; //animation start frame;
    u32 f_end = 1; //animation end frame;
    u32 frame_step = 1; //# frames between each frame;
    realNumber time_ramp = 240000.0_real;//scaling factor for time;
    realNumber camera_proper_time = 1.0/frame_rate;

    //scene
    realNumber celestial_sphere_radius = 400.0_real;

    //render mode
    bool render_debug = true;
    bool render_magik = false;
    bool render_magik_sandbox = false;

    // -- magik sandbox --
    //trivial
    bool magik_sandbox_trivial_scene = false;

    //single
    bool magik_sandbox_single_scene = false;

    //multiple
    bool magik_sandbox_muliple_scene = true;
    int magik_sandbox_monte_carlo_samples = 5120;
    int magik_sandbox_multiple_scatter_events = 12;
 
    // -- magik debugging --
    //3D Rasterize Volume Disk
    //3D Rasterize Volume Jet
    //3D Rasterize Volume Ambient

    // -- magik --
    //trivial

    //multiple

    // -- debugging --
    //test scenes
    bool test_scene = true;
    bool lattice_scene = false;

    //common
    bool colored_background = false;
    bool integration_depth_overlay = false;
    bool coordinate_overlay = false;
    bool Z_depth_overlay = false;
    bool world_grid = true; 
    bool world_gird_gravitational_lensing = false;
    bool world_gird_ignore_shadow = false; //true = world-grid is drawn over the shadow, false = it is not
    realNumber grid_spacing_major = 20.0_real;
    realNumber grid_spacing_minor = 5.0_real;

    //checker disk
    bool checker_disk_toggle = true; 
    bool checker_disk_Bi_sexual_color = false;
    bool checker_disk_Z_factor = false; 
    bool checker_disk_velocity_gradient = false;
    bool checker_disk_vortex_noise = true;
    bool checker_disk_vortex_noise_light_travel_delay = true;
    realNumber checker_disk_major_radius = 100.0_real;
    realNumber checker_disk_minor_radius = 5.0_real;

    //lattice point visualization
    bool toggle_disk_lattice_points = true; 
    bool toggle_jet_lattice_points = false;
    realNumber jet_lattice_point_size = 2.0_real;
    realNumber disk_lattice_point_size = 5.0_real;

    //BV Visualization
    bool vis_disk_BV = false;
    bool vis_disk_BV_hit = false;

    bool vis_jet_BV = false;
    bool vis_jet_BV_hit = false;

    bool vis_ambient_BV = false;
    bool vis_ambient_BV_hit = false;

    //Geodesic Export
    bool geodesic_export = false;
    bool geo_ex_UV = false;
    bool geo_ex_ARBITRARY = true;
    bool geo_ex_export_interpolated = false;
    realNumber geo_ex_interpolation_step = 0.1;
    realNumber geo_ex_MU = 0; //Massless vs Massive particle
    realNumber geo_ex_Orig_X = 10.0;
    realNumber geo_ex_Orig_Y = 0.0;
    realNumber geo_ex_Orig_Z = 0.0;
    realNumber geo_ex_Dir_X = -1.0;
    realNumber geo_ex_Dir_Y = 0.0;
    realNumber geo_ex_Dir_Z = 0.0;
    realNumber geo_ex_U1 = 0.0;
    realNumber geo_ex_U2 = 0.0;
    realNumber geo_ex_U3 = 0.0;

    // -- volumetrics --
    //quality disk
    u32 disk_noise_octaves = 8;
    u32 disk_noise_size_first_octave = 24;
    u32 disk_n_points = disk_noise_size_first_octave * std::pow(2,disk_noise_octaves-1);

    //disk lattice points
    realNumber disk_lattice_major_radius = 100.0_real; //90
    realNumber disk_lattice_minimum_radius = 5.0_real; //3
    realNumber disk_lattice_vertical_scale = 0.01_real; //0.05
    realNumber disk_lattice_velocity_scale = 1.0_real; //MUST BE >= 1

    //bounding volume disk
    bool bv_disk_toggle = false;
    realNumber bv_disk_major_radius = 110.0_real;
    realNumber bv_disk_minor_radius = 0.0_real;
    realNumber bv_disk_mean_slope = 12.0_real;
    realNumber bv_disk_radius_offset = 24.0_real;

    //quality jet
    u32 jet_noise_octaves = 2;
    u32 jet_noise_first_octave = 12;
    u32 jet_n_points = jet_noise_first_octave * std::pow(2,jet_noise_octaves-1);

    //bounding volume jet
    bool bv_jet_toggle = false;
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