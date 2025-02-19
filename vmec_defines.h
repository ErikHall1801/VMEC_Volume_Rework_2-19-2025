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

//Rendering 
    //Render target
    #define integrationDepth 3600
    #define SUPERSAMPLES 1
    #define IMAGE_HEIGHT 1000
    #define IMAGE_WIDTH 2000

    //Animation
    #define ANIMATE_FREE_FALL true
    #define ANIMATE_PATH false
    #define FRAME_RATE 24
    #define NUMBER_OF_FRAMES 100
    #define F_START 1 //Animation Start frame
    #define F_END 1 //Animation End frame
    #define FRAME_STEP 1 //# Frames between each Frame
    #define TIME_RAMP 2000.0 //Scaling factor for time. There is no right value for this because VMEC is written in natural units. It can just be interpreted as the unit of time, defined in terms of GM/c², being multipled by this value 

    //Scene 
    #define CELESTIAL_SPHERE_RADIUS 400.0

    //Render Mode
    #define RENDER_DEBUG true 
    #define RENDER_MAGIK false 
    #define RENDER_MAGIK_SANDBOX false 

//Debugging
    //Checker Disk
    #define CHECKER_DISK_TOGGLE true
    #define CHECKER_DISK_MAJOR_RADIUS 56.0
    #define CHECKER_DISK_MINOR_RADIUS 5.0

    //Lattice Point Visualization
        //Disk
        #define DISK_LATTICE_POINT_SIZE 1.0

        //Jet
        #define JET_LATTICE_POINT_SIZE 1.0

//Volumetrics
    //Quality Disk
    #define DISK_NOISE_OCTAVES 5
    #define DISK_NOISE_SIZE_FIRST_OCTAVE 24
    #define DISK_N_POINTS DISK_NOISE_SIZE_FIRST_OCTAVE * std::pow(2,DISK_NOISE_OCTAVES-1)

    //Disk lattice points
    #define DISK_LATTICE_MAJOR_RADIUS 90.0
    #define DISK_LATTICE_MINIMUM_RADIUS 3.0
    #define DISK_LATTICE_VERTICAL_SCALE 0.05
    #define DISK_LATTICE_VELOCITY_SCALE 5.0 //MUST BE >= 1

    //Bouding Volume Disk
    #define BV_DISK_TOGGLE true
    #define BV_DISK_MAJOR_RADIUS 120.0
    #define BV_DISK_MINOR_RADIUS 0.0
    #define BV_DISK_MEAN_SLOPE 12.0
    #define BV_DISK_RADIUS_OFFSET 24.0

    //Quality Jet
    #define JET_NOISE_OCTAVES 2
    #define JET_NOISE_FIRST_OCTAVE 12
    #define JET_N_POINTS JET_NOISE_FIRST_OCTAVE * std::pow(2,JET_NOISE_OCTAVES-1)

    //Bounding Volume Jet
    #define BV_JET_TOGGLE true
    #define BV_JET_HEIGHT 10.0
    #define BV_JET_DEFLECTION 22.0
    #define BV_JET_MINOR_RADIUS 16.0

    //Jet lattice points
    #define JET_MINOR_RADIUS 5.0
    #define JET_RADIUS_SCALE 10.0

    //Bounding Volume Ambient Medium
    #define BV_AMBIENT_MEDIUM_TOGGLE true
    #define BV_AMBIENT_MEDIUM_RADIUS 120.0

//Physics Related
    //Black Hole
    const realNumber a = -0.999_real; //Spin parameter

    //Astrophysical Jet
    const realNumber vy = 6.0_real;//Jet Vertical Velocity 40.5
    const realNumber w = 196.0_real; //Jet phi velocity scale 196
    const realNumber vr = 24.5_real; //Jet "Spread Factor 212
    const realNumber maxT = 8.0_real; //Max life time for jet particle 2
    const realNumber alpha = 2.0_real; //Jet phi velocity decay power. Keep it at 2

// -- Globals --

/*struct global_database
{
    // -- Rendering --
    //Render target
    u32 integration_depth = integrationDepth;
    u32 supersamples = SUPERSAMPLES;
    u32 image_height = IMAGE_HEIGHT;
    u32 image_width = IMAGE_WIDTH;

    //animation
    bool animate_free_fall = ANIMATE_FREE_FALL;
    bool animate_path = ANIMATE_PATH;
    u32 frame_rate = FRAME_RATE;
    u32 number_of_frames = NUMBER_OF_FRAMES;
    u32 f_start = F_START; //animation start frame;
    u32 f_end = F_END; //animation end frame;
    u32 frame_step = FRAME_STEP; //# frames between each frame;
    u32 time_ramp = TIME_RAMP;//scaling factor for time;

    //scene
    realNumber celestial_sphere_radius = CELESTIAL_SPHERE_RADIUS;

    //render mode
    bool render_debug = RENDER_DEBUG;
    bool render_magik = RENDER_MAGIK;

    // -- debugging --
    //checker disk
    bool checker_disk_toggle = CHECKER_DISK_TOGGLE;
    realNumber checker_disk_major_radius = CHECKER_DISK_MAJOR_RADIUS;
    realNumber checker_disk_minor_radius = CHECKER_DISK_MINOR_RADIUS;

    //lattice pou32 visualization
    realNumber jet_lattice_pou32_size = JET_LATTICE_POINT_SIZE;
    realNumber disk_lattice_pou32_size = DISK_LATTICE_POINT_SIZE;

    // -- volumetrics --
    //quality disk
    u32 disk_noise_octaves = DISK_NOISE_OCTAVES;
    u32 disk_noise_size_first_octave = DISK_NOISE_SIZE_FIRST_OCTAVE;

    //disk lattice pou32s
    realNumber disk_lattice_major_radius = DISK_LATTICE_MAJOR_RADIUS;
    realNumber disk_lattice_minimum_radius = DISK_LATTICE_MINIMUM_RADIUS;
    realNumber disk_lattice_vertical_scale = DISK_LATTICE_VERTICAL_SCALE;

    //bouding volume disk
    bool bv_disk_toggle = BV_DISK_TOGGLE;
    realNumber bv_disk_major_radius = BV_DISK_MAJOR_RADIUS;
    realNumber bv_disk_minor_radius = BV_DISK_MINOR_RADIUS;
    realNumber bv_disk_mean_slope = BV_DISK_MEAN_SLOPE;
    realNumber bv_disk_radius_offset = BV_DISK_RADIUS_OFFSET;

    //quality jet
    u32 jet_noise_octaves = JET_NOISE_OCTAVES;
    u32 jet_noise_first_octave = JET_NOISE_FIRST_OCTAVE;

    //bounding volume jet
    bool bv_jet_toggle = BV_JET_TOGGLE;
    realNumber bv_jet_height = BV_JET_HEIGHT;
    realNumber bv_jet_deflection = BV_JET_DEFLECTION;
    realNumber bv_jet_minor_radius = BV_JET_MINOR_RADIUS;

    //size jet
    realNumber jet_minor_radius = JET_MINOR_RADIUS;
    realNumber jet_radius_scale = JET_RADIUS_SCALE;

    //bounding volume ambient medium
    bool bv_ambient_medium_toggle = BV_AMBIENT_MEDIUM_TOGGLE;
    realNumber bv_ambient_medium_radius = BV_AMBIENT_MEDIUM_RADIUS;

//Physics Related
    //Black Hole
    realNumber a = -0.999_real; //Spin parameter

    //Astrophysical Jet
    realNumber vy = 6.0_real;//Jet Vertical Velocity 40.5
    realNumber w = 196.0_real; //Jet phi velocity scale 196
    realNumber vr = 24.5_real; //Jet "Spread Factor 212
    realNumber max_t = 8.0_real; //Max life time for jet particle 2
    realNumber alpha = 2.0_real; //Jet phi velocity decay power. Keep it at 2

};

extern global_database global;*/
#endif // VMEC_DEFINES_H