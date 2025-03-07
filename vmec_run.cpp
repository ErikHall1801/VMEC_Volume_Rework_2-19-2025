#define _CRT_SECURE_NO_DEPRECATE

#include "vmec_defines.h"
#include "vmec_types.h"
#include "vmec_common.h"
#include "vmec_setup.h"
#include "vmec_CPU.h"
#include "vmec_output.h"

namespace filesys = std::filesystem;
using std::to_string;

int main(int argc, char** argv) try
{
    global.args.reserve( argc );
    for (int i=0; i<argc; ++i)
    {
        static fstring x_arg;
        x_arg = argv[i];
        global.args.push_back( x_arg );
        if (x_arg == "--help")
        { settings.arg_help = true; vmec_log("arg: help"); }
        else if (x_arg == "--verbose")
        { settings.arg_verbose = true; vmec_log("arg: verbose"); }
        else if (x_arg == "--skip-render")
        { settings.arg_skip_render = true; vmec_log("arg: skip_render"); }
    }

    pix_realNumber *pixels_raw;
    pix_RGB *pixels_clean;
    RayProperties *all_geodesics;
    int *geodesic_sizes;
    realNumber *geodesics_k;	//Path length along geodesic
    int num_disk_lattice_points, num_jet_lattice_points;
    lattice_noise_point *FRAME_jet_lattice_noise_points;
    dynamic_lattice_point *dynamic_noise_points;
    noise_point *disk_lattice_noise_points;

    vmec_log("Beginning Setup");

    // -- Setup filesystem --
    settings.session_timestamp = std::chrono::system_clock::now().time_since_epoch().count();
    // Setup a directory for just this run of a render, might change if multiple runs are done
    settings.session_dir = settings.render_dir + "/version_0_" + to_string(settings.session_timestamp);

    vmec_log("Start Time: ", settings.session_timestamp);
    settings.render_binary_dir =  settings.session_dir + "/" + settings.render_binary_dir;
    settings.render_exr_dir = settings.session_dir + "/" + settings.render_exr_dir;
    settings.render_png_dir = settings.session_dir + "/" + settings.render_png_dir;
    filesys::create_directory(settings.render_dir);
    filesys::create_directory(settings.session_dir);
    filesys::create_directory(settings.render_binary_dir);
    filesys::create_directory(settings.render_exr_dir);
    filesys::create_directory(settings.render_png_dir);

    int total_threads=32;
    do_setup(total_threads, pixels_raw, pixels_clean, all_geodesics, geodesic_sizes, geodesics_k, num_jet_lattice_points, FRAME_jet_lattice_noise_points, dynamic_noise_points, num_disk_lattice_points, disk_lattice_noise_points);

    //Set up the initial sensor: camera view and properties
    SensorProperties sensorproperties;
    sensorproperties.maxWidth = settings.image_width;
    sensorproperties.maxHeight = settings.image_height;
    sensorproperties.SensorSize = 35.0; //Physical sensor size in mm, here Super 35
    sensorproperties.FocalLength = 70.0; //Camera FOV (focal length) in mm
    sensorproperties.FocalLength /= 1000.0; //Convert to m
    sensorproperties.camPos.x = 0.005; //Cartesian camera position 0.005, 0.0023, 86.23
    sensorproperties.camPos.y = 0.0023;
    sensorproperties.camPos.z = 86.23; 
    sensorproperties.camMomentum.x = -0.0; //Read r, theta, phi COORDINATE velocity
    sensorproperties.camMomentum.y = 0.0002;
    sensorproperties.camMomentum.z = 0.0; 
    sensorproperties.AngleX = -0.0123; //Pitch -0.0123, 0.064, 0.523
    sensorproperties.AngleY = 0.064; //Yaw
    sensorproperties.AngleZ = 0.823; //Roll

    //Camera Data for Free fall
    localProperties cameraProperties;
    RayProperties pointProperties, lastProperties;
    Vector3 camPos(sensorproperties.camPos.x,sensorproperties.camPos.y,sensorproperties.camPos.z), camMomentum(sensorproperties.camMomentum.x,sensorproperties.camMomentum.y,sensorproperties.camMomentum.z), camDir;

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
    bool rs_hit = false;

    //Camera Animation variables
    Vector3 XYZCam;
    realNumber hStepSum = 0.0, frameTime = 0.0, frameTimeStep = 1.0/realNumber(settings.frame_rate);

    //Time at infinity
    realNumber globalTime = 0.0;
    realNumber globalTimeOffset = 1000.0;

    try {
    //Animation Loop
    vmec_log("Entering Main Program Loop");
    log_flush();
    int frameNumber;
    for(frameNumber = 0; frameNumber < settings.number_of_frames; frameNumber++)
    {
        log_flush();
        vmec_log("Starting frame: ", frameNumber);
        if(settings.animate_path)
        {
            //F
        }

        if(settings.animate_free_fall && !settings.static_camera)
        {
            //Advance time for the first frame even if animation is static
            lastProperties = pointProperties;

            while(frameTimeStep*settings.time_ramp > hStepSum && !rs_hit && !settings.static_camera)
            {
                if(pointProperties.hStep > 1e-5) {pointProperties.hStep = 1e-5;}
                pointProperties = rayRKF45(pointProperties);

                if(isinf(pointProperties.hStep))
                {
                    pointProperties.hStep = 1e-6;
                }

                hStepSum += pointProperties.hStep;

                if(pointProperties.r < settings.rs*settings.rs_scale*2.0)
                {
                    rs_hit = true;
                }
            }

            vmec_log("hStep" , 1.0 / (hStepSum/settings.time_ramp));

            globalTime = pointProperties.t;
            XYZCam = BoyerLindquistToCartesian(pointProperties.r,pointProperties.theta,pointProperties.phi);

            sensorproperties.camPos.x = XYZCam.x;
            sensorproperties.camPos.y = XYZCam.y;
            sensorproperties.camPos.z = XYZCam.z;
            sensorproperties.camMomentum.x = pointProperties.u1;
            sensorproperties.camMomentum.y = pointProperties.u2;
            sensorproperties.camMomentum.z = pointProperties.u3;
        }

        if(rs_hit)
        {
            vmec_log("Geodesic hit event horizon");
            break;
        }

        if(settings.bv_jet_toggle)
        {
            //F
        }

        if(frameNumber >= settings.f_start && frameNumber <= settings.f_end && frameNumber % settings.frame_step == 0)
        {
            if (settings.arg_skip_render == false)
            {
                trace_all_rays(globalTime+globalTimeOffset, total_threads, sensorproperties, all_geodesics, geodesic_sizes, geodesics_k, disk_lattice_noise_points, num_jet_lattice_points, FRAME_jet_lattice_noise_points, dynamic_noise_points, pixels_raw);
            }

            //Convert raw (doubles) pixels to proper RGB
            write_binary_realNumber( pixels_raw, frameNumber/settings.frame_step );
            write_hdr_realNumber( pixels_raw, frameNumber/settings.frame_step );
            write_exr_realNumber( pixels_raw, frameNumber/settings.frame_step );
        }

        if(frameNumber >= settings.f_end)
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

    } catch (std::exception& e)
    { vmec_log("EXCEPTION: ", e.what()); }

    vmec_log("Starting Program Cleanup");
    do_cleanup(pixels_raw, pixels_clean, all_geodesics, geodesic_sizes, geodesics_k, FRAME_jet_lattice_noise_points, dynamic_noise_points, disk_lattice_noise_points);
    vmec_log("Program Exiting Cleanly");
}
catch (std::exception& e)
{ vmec_log("UNCAUGHT EXCEPTION: ", e.what(), "\nExiting abruptly, might leak or corrupt"); }
