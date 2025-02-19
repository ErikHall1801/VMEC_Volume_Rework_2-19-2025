#define _CRT_SECURE_NO_DEPRECATE
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>

#include "vmec_defines.h"
#include "vmec_types.h"
#include "vmec_common.h"
#include "vmec_setup.h"
#include "vmec_CPU.h"
#include "vmec_output.h"


int main(){

    pix_realNumber *pixels_raw;
    pix_RGB *pixels_clean;
    int num_disk_lattice_points, num_jet_lattice_points;
    lattice_noise_point *FRAME_jet_lattice_noise_points;
    noise_point *disk_lattice_noise_points;

    vmec_log("Beginning Setup");
    do_setup(pixels_raw, pixels_clean, num_jet_lattice_points, FRAME_jet_lattice_noise_points, num_disk_lattice_points, disk_lattice_noise_points);

    //Set up the initial sensor: camera view and properties
    SensorProperties sensorproperties;
    sensorproperties.maxWidth = IMAGE_WIDTH;
    sensorproperties.maxHeight = IMAGE_HEIGHT;
    sensorproperties.SensorSize = 35.0; //Physical sensor size in mm, here Super 35
    sensorproperties.FocalLength = 24.0; //Camera FOV (focal length) in mm
    sensorproperties.FocalLength /= 1000.0; //Convert to m
    sensorproperties.camPos.x = 0.005; //Cartesian camera position
    sensorproperties.camPos.y = 16.5023;
    sensorproperties.camPos.z = 56.23; 
    sensorproperties.camMomentum.x = -0.0; //Read r, theta, phi COORDINATE velocity
    sensorproperties.camMomentum.y = -0.0;
    sensorproperties.camMomentum.z = 0.0; 
    sensorproperties.AngleX = -16.6123; //Pitch
    sensorproperties.AngleY = 0.064; //Yaw
    sensorproperties.AngleZ = 0.023; //Roll

    //Camera Data for Free fall
    localProperties cameraProperties;
    RayProperties pointProperties;
    Vector3 camPos(sensorproperties.camPos.x,sensorproperties.camPos.y,sensorproperties.camPos.z), camMomentum(sensorproperties.camMomentum.x,sensorproperties.camMomentum.y,sensorproperties.camMomentum.z), camDir;
    BVObjects emptyBV;

    emptyBV.disk = false; emptyBV.jet = false; emptyBV.ambient = false;

    cameraProperties = initializeCamera(camPos,camMomentum);
    pointProperties.t = cameraProperties.t;
    pointProperties.r = cameraProperties.r;
    pointProperties.theta = cameraProperties.theta;
    pointProperties.phi = cameraProperties.phi;
    pointProperties.u0 = cameraProperties.u0;
    pointProperties.u1 = cameraProperties.u1;
    pointProperties.u2 = cameraProperties.u2;
    pointProperties.u3 = cameraProperties.u3;
    pointProperties.hStep = 1.0e-15;

    vmec_log("Finished Setup");

    //Start timing
    auto start = std::chrono::high_resolution_clock::now();

    //Multithreading
    int total_threads=32;

    //Camera Animation variables
    Vector3 XYZCam;
    realNumber hStepSum = 0.0, frameTime = 0.0, frameTimeStep = 1.0/realNumber(FRAME_RATE);

    //Time at infinity
    realNumber globalTime = 0.0;
    realNumber globalTimeOffset = 1000.0;

    //Animation Loop
    vmec_log("Entering Main Program Loop");
    log_flush();
    for(int frameNumber = 0; frameNumber < NUMBER_OF_FRAMES; frameNumber++)
    {
        log_flush();
        vmec_log("Starting frame: ", frameNumber);
        if(ANIMATE_PATH)
        {
            //F
        }

        if(ANIMATE_FREE_FALL)
        {
            while(frameTimeStep*TIME_RAMP > hStepSum)
            {
                pointProperties.hStep = 1e-6*TIME_RAMP;
                pointProperties = rayRKF45(pointProperties,emptyBV,0.0);
                hStepSum += pointProperties.hStep;
            }

            globalTime = pointProperties.t;
            XYZCam = BoyerLindquistToCartesian(pointProperties.r,pointProperties.theta,pointProperties.phi);

            sensorproperties.camPos.x = XYZCam.x;
            sensorproperties.camPos.y = XYZCam.y;
            sensorproperties.camPos.z = XYZCam.z;
            sensorproperties.camMomentum.x = pointProperties.u1;
            sensorproperties.camMomentum.y = pointProperties.u2;
            sensorproperties.camMomentum.z = pointProperties.u3;
        }

        if(BV_JET_TOGGLE)
        {
            //F
        }

        if(frameNumber >= F_START && frameNumber <= F_END && frameNumber % FRAME_STEP == 0)
        {
            trace_all_rays(globalTime+globalTimeOffset, total_threads, sensorproperties, disk_lattice_noise_points, num_jet_lattice_points, FRAME_jet_lattice_noise_points, pixels_raw);

            //Convert raw (doubles) pixels to proper RGB
            write_binary_realNumber(pixels_raw,frameNumber/FRAME_STEP);
        }

        if(frameNumber >= F_END)
        {
            break;
        }

        hStepSum = 0.0;
        frameTime += frameTimeStep;
    }

    vmec_log("End of animation loop");
    //Stop timing
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Time [ms / min]: " << duration.count() / 1000 << " / " << duration.count() / 60000000 << std::endl;

    vmec_log("Starting Program Cleanup");
    do_cleanup(pixels_raw, pixels_clean, FRAME_jet_lattice_noise_points, disk_lattice_noise_points);
    vmec_log("Program Exiting Cleanly");
}