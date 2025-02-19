#include "vmec_CPU.h"


//RNG Points Jet Lattice
static void initialize_per_RAY_lattice_points(const int num_lattice_points, dynamic_lattice_point* RAY_jet_lattice_noise_points, lattice_noise_point* FRAME_jet_lattice_noise_points){

	for(int idx=0;idx<num_lattice_points;++idx)
	{
		RAY_jet_lattice_noise_points[idx].r = FRAME_jet_lattice_noise_points[idx].r;
		RAY_jet_lattice_noise_points[idx].theta = FRAME_jet_lattice_noise_points[idx].t;
		RAY_jet_lattice_noise_points[idx].phi = FRAME_jet_lattice_noise_points[idx].p;;
		RAY_jet_lattice_noise_points[idx].tau = 0.0;
	}
}


//Search for start and end points in sorted array
static void traverseSearch(int startIndex, int endIndex, float target_low, float target_high, int* out_lower, int* out_upper, noise_point* array) {
    int i = startIndex;
    while (i < endIndex && array[i].r < target_low)
    {
        ++i;
    }
    *out_lower = std::max(i, 0);

    while (i < endIndex && array[i].r < target_high)
    {
        ++i;
    }
    *out_upper = std::min(i, endIndex - 1);
}


//Distance between Sample and disk lattice point
static Vector3 evaluateDiskLatticePoint(const realNumber globalTime, const int i, Vector3 sPos, noise_point* disk_lattice_noise_points, const realNumber t, const realNumber r, const realNumber theta)
{
    Vector4 diskLatticePoint;
    diskLatticePoint.x = disk_lattice_noise_points[i].r;
    diskLatticePoint.y = disk_lattice_noise_points[i].t;
    diskLatticePoint.z = disk_lattice_noise_points[i].p;
    diskLatticePoint.t = disk_lattice_noise_points[i].v;

    Vector3 uVelPoint, xPos;
    realNumber gamma, angularDisplacement, dr, dtheta, dphi, sampleDistance;        
    dr = diskLatticePoint.x * DISK_LATTICE_MAJOR_RADIUS + DISK_LATTICE_MINIMUM_RADIUS;
            
    uVelPoint.set(0.0,0.0,uPhi(dr));
    gamma = generalizedTimeDilationFactor(r, theta, uVelPoint, -1.0);
    angularDisplacement = (globalTime - (t/gamma)) * -uVelPoint.z * sign(a);        
    dtheta = (diskLatticePoint.y * DISK_LATTICE_VERTICAL_SCALE) + (M_PI/2.0);
    dphi = diskLatticePoint.z + angularDisplacement;
            
    return sphericalToCartesian(dr, dtheta, dphi);
}


//Vortex Noise algorithm
static Vector3 vortexNoise(const realNumber globalTime, const Vector3& rayOrig, int octaves, int sizeFirstOctave, noise_point* disk_lattice_noise_points, realNumber previousValue, realNumber t, realNumber r, realNumber theta, realNumber phi, realNumber dx)
{
    Vector4 pointData;
    Vector3 xPos, xPosDot, sPos, sPosDot, uVelPoints, out;
    realNumber gamma, initialAngularDisplacement, constWeight = 0.0, variableWeight = 0.0, vWeight = 0.0, cWeight = 0.0, vValue = 0.0, cValue = 0.0, vortexNoiseValue = 0.0, constantVortexNoiseValue = 0.0, vortexNoiseDerivative = 0.0, NthHarmonic = 0.0, angularDisplacement = 0.0, dr, dtheta, dphi;


    //SPT
    sPos = rayOrig;
    sPos.y = sPos.y * 3.0;
    initialAngularDisplacement = std::sqrt(1.0 / ((r+3.0) * (r+3.0))) * -sign(a) * 156.0;
    sPos = rotate2D(sPos, initialAngularDisplacement);
    sPosDot.x = sPos.x;
    sPosDot.y = 0.0;
    sPosDot.z = sPos.z;

    //Temp
    realNumber diskSize = 93.0;
    realNumber animationScale = 10.0;
    realNumber lightTravelDelayExaggeration = 1.0;
    realNumber weightThreshold = 0.000000000001;
    realNumber maxWidth = 0.0;
    realNumber lowerBound, upperBound;
    int octaveIndices = 0, octaveStart = 0, octaveEnd = 0, i_start = octaveStart, i_end = octaveEnd;

    //For # octaves
    for(int o = 0; o < octaves; o++)
    {
        vWeight = 0.0; cWeight = 0.0; vValue = 0.0; cValue = 0.0;
        octaveIndices = sizeFirstOctave*std::pow(2,o);
        octaveStart = octaveEnd;
        octaveEnd = octaveStart+octaveIndices;

        maxWidth = std::pow(1.0/weightThreshold,1.0/(2.0+o)) / 120.0;
        lowerBound = std::max(((r-1.0)/120.0)-maxWidth,0.0);
        upperBound = std::min(((r-1.0)/120.0)+maxWidth,1.0);

        traverseSearch(octaveStart, octaveEnd, lowerBound, upperBound, &i_start, &i_end, disk_lattice_noise_points);

        //for # points in octave within the lower and upper bound
        for(int g = i_start; g < i_end; g++)
        {
            pointData.x = disk_lattice_noise_points[g].r;
            pointData.y = disk_lattice_noise_points[g].t;
            pointData.z = disk_lattice_noise_points[g].p;
            pointData.t = disk_lattice_noise_points[g].v;

            dr = pointData.x * diskSize + 3.0;

            uVelPoints.x = 0.0;
            uVelPoints.y = 0.0;
            uVelPoints.z = uPhi(dr)* lightTravelDelayExaggeration;
            gamma = generalizedTimeDilationFactor(r, theta,uVelPoints, -1.0);
            angularDisplacement = ((globalTime*animationScale) - (t/gamma)) * -uVelPoints.z * sign(a);

            dtheta = (pointData.y * 0.267) + (M_PI/2.0); //0.267 is a scale such that all points are within the boudning region, this should be functionalized
            dphi = pointData.z + angularDisplacement;

            xPos = sphericalToCartesian(dr, dtheta, dphi); //XYZ with Y up
            xPosDot.x = xPos.x;
            xPosDot.y = 0.0;
            xPosDot.z = xPos.z;

            variableWeight = 1.0 / std::pow(euclidianDistance2(xPos,sPos),(2.0+o)/2.0);

            if(variableWeight > weightThreshold);
            {
                constWeight = 1.0 / std::pow(euclidianDistance2(xPosDot,sPosDot),(2.0+o)/2.0);

                vWeight += variableWeight;
                cWeight += constWeight;
                vValue += pointData.t*variableWeight;
                cValue += pointData.t*constWeight;
            }
        }

        vortexNoiseValue += (vValue/vWeight)*(1.0/static_cast<realNumber>(o+1.0));
        constantVortexNoiseValue += (cValue/cWeight)*(1.0f/static_cast<realNumber>(o+1.0));
        NthHarmonic += 1.0/static_cast<realNumber>(o+1);
    }

    vortexNoiseValue /= NthHarmonic;
    constantVortexNoiseValue /= NthHarmonic;
    vortexNoiseDerivative = (vortexNoiseValue - previousValue) / dx;
    out.x = vortexNoiseValue;
    out.y = constantVortexNoiseValue;
    out.z = vortexNoiseDerivative;
    return out;
}


//Jet Noise field generation
static Vector3 maelstromNoise(const Vector3& rayOrig, const Vector3& Vel, realNumber dx, realNumber previousValue, dynamic_lattice_point* RAY_jet_lattice_noise_points, lattice_noise_point* FRAME_jet_lattice_noise_points)
{
    Vector3 pointBltoCartesian, xPos, sPos, out;
    realNumber mappingFactor = 2.0, polarPhi = 0.0, polarR = 0.0, X = 0.0, r = 0.0, theta = 0.0, phi = 0.0, distY = 0.0, mWeight = 0.0, mValue = 0.0, NthHarmonic = 0.0, variableWeight = 0.0, maelstromNoiseValue = 0.0, maelstromNoiseDerivative = 0.0, factor = 0.0;
    int octaveIndices, octaveStart = 0, octaveEnd = 0;
    out.x = 0.0;
    out.y = 0.0;
    out.z = 0.0;

    //Bokeh Noise Algorithm

    //SPT
    distY = std::sqrt(vr*(std::abs(rayOrig.y)/vy) + 4.0); //Min y distance
    factor = (1.0 / std::sqrt( rayOrig.x*rayOrig.x + rayOrig.z*rayOrig.z )) * 0.25;
    sPos = rayOrig;
    sPos = rotate2D(sPos,w*factor);

    for (int o = 0; o < JET_NOISE_OCTAVES; o++)
    {
        mWeight = 0, mValue = 0;
        octaveIndices = JET_NOISE_FIRST_OCTAVE*std::pow(2,o);
        octaveStart = octaveEnd;
        octaveEnd = octaveStart+octaveIndices;

        for(int i = octaveStart; i < octaveEnd; i++)
        {
            //Extract point Position and Value
            X = (FRAME_jet_lattice_noise_points[i].v+1.0)*0.5;
            r = RAY_jet_lattice_noise_points[i].r;
            theta = RAY_jet_lattice_noise_points[i].theta;
            phi = RAY_jet_lattice_noise_points[i].phi;

            //Bl to Cartesian
            xPos = BoyerLindquistToCartesian(r,theta,phi);

            //Weight
            variableWeight = 1.0 / std::pow(euclidianDistance2(xPos,sPos),(2.0+o)/2.0);
            mWeight += variableWeight;
            mValue += X*variableWeight;
        }

        maelstromNoiseValue += (mValue/mWeight)*(1.0/static_cast<realNumber>(o+1.0));
        NthHarmonic += 1.0/static_cast<realNumber>(o+1);
    }

    maelstromNoiseValue /= NthHarmonic;
    maelstromNoiseDerivative = (maelstromNoiseValue - previousValue) / dx;

    //Out
    out.x = maelstromNoiseValue;
    out.y = maelstromNoiseDerivative;
    out.z = 0.0;
    return out;
}


//Advanced Numerical Approximation
static void ANA(realNumber globalTime, realNumber dt, lattice_noise_point* FRAME_jet_lattice_noise_points, dynamic_lattice_point* RAY_jet_lattice_noise_points, int mPoints)
{
    for(int k = 0; k < mPoints; k++)
    {
        //Initials & Wrapping
        realNumber verticalDirection = (FRAME_jet_lattice_noise_points[k].s);
        realNumber vtau = RAY_jet_lattice_noise_points[k].tau;
        realNumber phi0 = FRAME_jet_lattice_noise_points[k].p;
        realNumber r0 = (FRAME_jet_lattice_noise_points[k].r * JET_RADIUS_SCALE) + JET_MINOR_RADIUS;

        //Wrapping
        realNumber Iter = std::trunc(vtau / maxT);
        realNumber truncTime = maxT*Iter;
        vtau -= truncTime;

        //Position and Velocity
        realNumber ftau = std::sqrt(vr*vtau + r0*r0);
        realNumber rtauSqrtTerm1 = std::sqrt( ((vy*vy*vtau*vtau - a*a + ftau*ftau)*(vy*vy*vtau*vtau - a*a + ftau*ftau)) + 4.0*vy*vy*a*a*vtau*vtau );
        realNumber rtauSqrtTerm2 = std::sqrt( vy*vy*vtau*vtau - a*a + ftau*ftau +  rtauSqrtTerm1 );

        //Position
        realNumber rtau = (1.0/std::sqrt(2.0)) * rtauSqrtTerm2;
        realNumber thetatau = std::acos((vy*vtau/rtau*sign(verticalDirection)));

        Vector3 u = uATJ(vtau,verticalDirection,rtau,thetatau);

        //Time Dilation Factor
        realNumber Gamma = generalizedTimeDilationFactor(rtau,thetatau,u,-1.0);

        //Update point
        RAY_jet_lattice_noise_points[k].tau -= dt/Gamma;
        RAY_jet_lattice_noise_points[k].r = rtau;
        RAY_jet_lattice_noise_points[k].theta = thetatau;
        RAY_jet_lattice_noise_points[k].phi -= u.z*dt/Gamma;
    }
}


//Geodesic Integration - does one integration pass to identify if the ray hits disk/jet/wind,
//then second pass to obtain the emission at every step (much more time-consuming) from the identified feature
static Vector4 traceRay_twopass(const realNumber globalTime, const Vector3& rayOrig, const Vector3& rayDir, const Vector3& camPos, const Vector3& camMomentum, noise_point* disk_lattice_noise_points, int num_lattice_points, lattice_noise_point* FRAME_jet_lattice_noise_points, dynamic_lattice_point* RAY_jet_lattice_noise_points)
{
    int i, hit_disk_ID = 0, hit_jet_ID = 0;
    bool hit_disk = false, hit_jet = false, hit_ambient = false, hit_disk_lattice = false, hit_jet_lattice = false, hit_disk_BV = false, hit_jet_BV = false, hit_ambient_BV = false;

    realNumber
        hStep = 1e-3,
        dx = 1.0,
        dotCheck = 0.0,
        checker = 0.0,
        tempt = 0.0,
        Medium = 1.0,
        dt = 0.0;

    Vector3
        dir(rayDir.x,rayDir.y,rayDir.z), 
        XYZpos(rayOrig.x,rayOrig.y,rayOrig.z),
        previousXYZpos(0.0,0.0,0.0),
        backgroundCol(0.0,0.0,0.0),
        ZAMODir(0.0,0.0,0.0), 
        checkerDiskNormal(0.0,1.0,0.0),
        checkerDiskp0(0.0,0.0,0.0),
        pointPosition(0.0,0.0,0.0);

    localProperties
        cameraProperties = initializeCamera(camPos,camMomentum);

    RayProperties
        rayProperties = initializeRay(rayOrig,rayDir,cameraProperties);
        rayProperties.hStep = hStep;

    BVObjects
        BoundingVolumes;

    Vector4
        initialU(rayProperties.u0,rayProperties.u1,rayProperties.u2,rayProperties.u3), 
        RGBA_out(0.0,0.0,0.0,0.0), 
        ZAMOCheck(0.0,0.0,0.0,0.0);

    //Render Mode
    bool toggleDebugView = RENDER_DEBUG;
    bool toggleMagik = RENDER_MAGIK;

        //Render Debug
            //Scene Selection
            bool toggleTestScene = true;
            bool toggleLatticeScene = false;

            //Relative to Both Scenes
            bool toggleColorBackground = true;

                //Relative to Test scene
                bool toogleIntegrationDepthOverlay = false; //Pixel color = Number of steps taken by a ray
                bool toogleCoordinateOverlay = false; //Pixel color = coordinates of the rays final position (R = r, G = theta, B = phi)
                bool toggleAberrationDot = false; //Pixel color = dot product between the pre and post aberration ray direction vectors
                bool toggleDiskVisualization = true; //Renders a 2D checker board disk

                //Relative to Lattice Points
                bool toggleDiskBV = true;
                bool toggleJetBV = false;
                bool toggleAmbientBV = false;

                    //Disk
                    bool toggleDiskLatticePoints = false;

                    //Jet
                    bool toggleJetLatticePoints = false;

                        //Realtive to Both
                        bool toggleLocalZAMOVelocity = false;

        //Render Magik
            //Method Selection

    //First integration pass - detect if disk/wind/jet is hit
    for(i = 0; i <= integrationDepth; i++)
    {
           
        if(rayProperties.r < (1.0+std::sqrt(1.0-(std::abs(a)*std::abs(a))))*1.0001) //Event Horizon
        {
            break;
        }
        else if(rayProperties.r > CELESTIAL_SPHERE_RADIUS || i+1 == integrationDepth) //Celestial Sphere
        {
            break;
        }
        else if(hit_disk && hit_jet && hit_ambient) //All features hit; might as well stop this pass
        {
            break;
        }

        //In-bounding region check
        pointInBounds(BoundingVolumes,rayProperties.r,rayProperties.theta);

        if(BoundingVolumes.disk)
        {
            hit_disk = true;
        }
        if(BoundingVolumes.jet)
        {
            hit_jet = true;
        }
        if(BoundingVolumes.ambient)
        {
            hit_ambient = true;
        }

        //Store current ray time coordinate
        tempt = rayProperties.t;

        //March step
        rayProperties = rayRKF45(rayProperties,BoundingVolumes,Medium);

        //Calculate dt
        dt = rayProperties.t-tempt;

        //Calculate dx
        previousXYZpos = XYZpos;
        XYZpos = BoyerLindquistToCartesian(rayProperties.r,rayProperties.theta,rayProperties.phi);
        dx = euclidianDistance(previousXYZpos,XYZpos);
        dir = XYZpos - previousXYZpos;
        dir = dir / VectorLength(dir);
    }

    //VMEC Debug Viewer
    if(toggleDebugView)
    {
        cameraProperties = initializeCamera(camPos,camMomentum);
        rayProperties = initializeRay(rayOrig,rayDir,cameraProperties);
        rayProperties.hStep = hStep;

        initialU.set(rayProperties.u0,rayProperties.u1,rayProperties.u2,rayProperties.u3);
        dir.set(rayDir.x,rayDir.y,rayDir.z);
        XYZpos.set(rayOrig.x,rayOrig.y,rayOrig.z);

        //Integrate geodesic - Test Scene
        if(toggleTestScene)
        {
            for(i = 0; i <= integrationDepth; i++)
            {
                if(toggleAberrationDot) //Absolute Visualization
                {
                    ZAMOCheck = restframeToZAMO(rayDir,cameraProperties);
                    ZAMODir.set(ZAMOCheck.x,ZAMOCheck.y,ZAMOCheck.z);
                    ZAMODir = ZAMODir / VectorLength(ZAMODir);
                    dotCheck = dotProduct(rayDir,ZAMODir);
                    RGBA_out.set(dotCheck,dotCheck,dotCheck,0.0);
                    break;
                }

                if(rayProperties.r < (1.0+std::sqrt(1.0-(std::abs(a)*std::abs(a))))*1.01) //Event Horizon
                {
                    RGBA_out.x += 0.0;
                    RGBA_out.y += 0.0;
                    RGBA_out.z += 0.0;
                    RGBA_out.t += 0.0;
    
                    if(toogleIntegrationDepthOverlay) //Absolute Visualization
                    {
                        RGBA_out.set(i,i,i,0.0);
                    }
    
                    if(toogleCoordinateOverlay) //Absolute Visualization
                    {
                        RGBA_out.set(rayProperties.r,rayProperties.theta,rayProperties.phi,0.0);
                    }
    
                    break;
                }

                if(rayProperties.r > CELESTIAL_SPHERE_RADIUS) //Celestial Sphere
                {
                    if(toggleColorBackground)
                    {
                        XYZpos = XYZpos / VectorLength(XYZpos);
                        backgroundCol.set((XYZpos.x + 1.0) / 2.0, (XYZpos.y + 1.0) / 2.0, (XYZpos.z + 1.0) / 2.0);
                    }
                    
                    RGBA_out.x += backgroundCol.x;
                    RGBA_out.y += backgroundCol.y;
                    RGBA_out.z += backgroundCol.z;
                    RGBA_out.t = 0.0;
    
                    if(toogleIntegrationDepthOverlay) //Absolute Visualization
                    {
                        RGBA_out.set(i,i,i,0.0);
                    }

                    if(toogleCoordinateOverlay) //Absolute Visualization
                    {
                        RGBA_out.set(rayProperties.r,rayProperties.theta,rayProperties.phi,0.0);
                    }

                    break;
                }

                if(i+1 == integrationDepth) //Insufficient depth //Absolute Visualization
                {
                    RGBA_out.set(1.0,0.0,0.0,0.0);
                    break;
                }

                if(rayDiskIntersection(checkerDiskNormal,checkerDiskp0,previousXYZpos,dir,dx) && toggleDiskVisualization)// Checker disk
                {
                    checker = radialChecker(rayProperties.r, rayProperties.phi, 56.0, 0.1, 5.0);
                    XYZpos.y = 0.0;
                    XYZpos = XYZpos / VectorLength(XYZpos);
                    RGBA_out.x += checker * ((XYZpos.x + 1.0) / 2.0);
                    RGBA_out.y += checker * ((XYZpos.y));
                    RGBA_out.z += checker * ((XYZpos.z));
                    RGBA_out.t = 0.0;

                    if(toogleIntegrationDepthOverlay) //Absolute Visualization
                    {
                        RGBA_out.set(i,i,i,0.0);
                    }

                    if(toogleCoordinateOverlay) //Absolute Visualization
                    {
                        RGBA_out.set(rayProperties.r,rayProperties.theta,rayProperties.phi,0.0);
                    }

                    break;
                }

                //March step
                rayProperties = rayRKF45(rayProperties,BoundingVolumes,Medium);

                //Calculate dx
                previousXYZpos = XYZpos;
                XYZpos = BoyerLindquistToCartesian(rayProperties.r,rayProperties.theta,rayProperties.phi);
                dx = euclidianDistance(previousXYZpos,XYZpos);
                dir = XYZpos - previousXYZpos;
                dir = dir / VectorLength(dir);
            }
        }

        if(toggleLatticeScene)
        {   
            if(hit_jet)
            {
                //Create Per ray copy of FRAME array
                initialize_per_RAY_lattice_points(num_lattice_points, RAY_jet_lattice_noise_points, FRAME_jet_lattice_noise_points);
            }

            for(i = 0; i <= integrationDepth; i++)
            {
                if(hit_disk || hit_jet || hit_ambient)
                {
                    pointInBounds(BoundingVolumes,rayProperties.r,rayProperties.theta);

                    if(toggleDiskBV && hit_disk_BV != true && BoundingVolumes.disk)
                    {
                        RGBA_out.x += 0.05;
                        hit_disk_BV = true;
                    }

                    if(toggleJetBV && hit_jet_BV != true && BoundingVolumes.jet)
                    {
                        RGBA_out.y += 0.05;
                        hit_jet_BV = true;
                    }

                    if(toggleAmbientBV && hit_ambient_BV != true && BoundingVolumes.ambient)
                    {
                        RGBA_out.z += 0.05;
                        hit_ambient_BV = true;
                    }
                }

                if(toggleDiskLatticePoints && hit_disk) //Disk Lattice points
                {
                    if(BoundingVolumes.disk)
                    {
                        for(int g = 0; g < DISK_N_POINTS; g++)
                        {
                            pointPosition = evaluateDiskLatticePoint(globalTime, g, XYZpos, disk_lattice_noise_points, rayProperties.t, rayProperties.r, rayProperties.theta);

                            if(euclidianDistance(XYZpos,pointPosition) < DISK_LATTICE_POINT_SIZE)
                            {
                                hit_disk_lattice = true;
                                hit_disk_ID = g;
                                break;
                            }
                        }
                    }
                }

                if(hit_disk_lattice)
                {
                    RGBA_out.x += ndetermin(hit_disk_ID);
                    RGBA_out.y += ndetermin(hit_disk_ID+465);
                    RGBA_out.z += ndetermin(hit_disk_ID+5742);
                    RGBA_out.t = 0.0;
                    RGBA_out = RGBA_out/VectorLength(RGBA_out);
                    break;
                }

                if(toggleJetLatticePoints && hit_jet) //Jet Lattice points
                {
                    if(BoundingVolumes.jet)
                    {
                        for(int g = 0; g < JET_N_POINTS; g++)
                        {
                            if(euclidianDistance(XYZpos,BoyerLindquistToCartesian(RAY_jet_lattice_noise_points[g].r,RAY_jet_lattice_noise_points[g].theta,RAY_jet_lattice_noise_points[g].phi)) < JET_LATTICE_POINT_SIZE)
                            {
                                hit_jet_lattice = true;
                                hit_jet_ID = g;
                                break;
                            }
                        }
                    }
                }

                if(hit_jet_lattice)
                {
                    RGBA_out.x += ndetermin(hit_jet_ID+5647);
                    RGBA_out.y += ndetermin(hit_jet_ID+68934);
                    RGBA_out.z += ndetermin(hit_jet_ID+923845);
                    RGBA_out.t = 0.0;
                    RGBA_out = RGBA_out/VectorLength(RGBA_out);

                    if(toggleLocalZAMOVelocity)
                    {
                        //F
                    }

                    break;
                }
                
                if(rayProperties.r < (1.0+std::sqrt(1.0-(std::abs(a)*std::abs(a))))*1.01) //Event Horizon
                {
                    RGBA_out.x += 0.0;
                    RGBA_out.y += 0.0;
                    RGBA_out.z += 0.0;
                    RGBA_out.t = 0.0;
                    break;
                }
    
                if(rayProperties.r > CELESTIAL_SPHERE_RADIUS) //Celestial Sphere
                {                    
                    RGBA_out.x += backgroundCol.x;
                    RGBA_out.y += backgroundCol.y;
                    RGBA_out.z += backgroundCol.z;
                    RGBA_out.t += 0.0;
                    break;
                }
    
                if(i+1 == integrationDepth) //Insufficient depth
                {
                    RGBA_out.set(1.0,0.0,0.0,0.0);
                    break;
                }

                //Store current ray time coordinate
                tempt = rayProperties.t;

                //March step
                rayProperties = rayRKF45(rayProperties,BoundingVolumes,Medium);

                //Calculate dt
                dt = rayProperties.t-tempt;
    
                //Calculate dx
                previousXYZpos = XYZpos;
                XYZpos = BoyerLindquistToCartesian(rayProperties.r,rayProperties.theta,rayProperties.phi);
                dx = euclidianDistance(previousXYZpos,XYZpos);
                dir = XYZpos - previousXYZpos;
                dir = dir / VectorLength(dir);
            }
        }
    }

    return RGBA_out;
}


//Magik Sandbox
static Vector4 traceRay_sandbox(const Vector3 rayOrig, const Vector3 rayDir)
{    
    Vector3
        AABBmax(10.0,10.0,10.0),
        AABBmin(-10.0,-10.0,-10.0),
        hit1(0.0,0.0,0.0),
        hit0(0.0,0.0,0.0),
        backgroundCol(1.0,1.0,1.0);

    Vector4
        RGBA_out(0.0,0.0,0.0,0.0);

    if(intersectAABB(rayOrig,rayDir,AABBmax,AABBmin,hit1,hit0))
    {
        realNumber dx = euclidianDistance(hit1,hit0);
        realNumber T = std::exp(-dx*0.25);
        backgroundCol = backgroundCol*T;
        //RGBA_out.set((hit0.x + AABBmax.x) / (2.0*AABBmax.x),(hit0.y + AABBmax.y) / (2.0*AABBmax.x),(hit0.z + AABBmax.z) / (2.0*AABBmax.x),0.0);
    }

    RGBA_out.set(backgroundCol.x,backgroundCol.y,backgroundCol.z,0.0);
    return RGBA_out;
}

//Main Pixel Shader
static pix_realNumber pixel_function(const realNumber globalTime, int width, int height, const SensorProperties sensorproperties, noise_point *disk_lattice_noise_points, int num_lattice_points, lattice_noise_point *FRAME_jet_lattice_noise_points,  dynamic_lattice_point* RAY_jet_lattice_noise_points)
{
    pix_realNumber pixel;

    Vector3 fragmentCoords;
    Vector3 rayOrig, rayDir(0.0,0.0,0.0);
    Vector4 mRender(0.0,0.0,0.0,0.0);
    Matrix3 M;

    //1. Convert 1D pixel array to 2D U,V coordinates to work for arbitary aspect ratios
    realNumber ratio = realNumber(sensorproperties.maxHeight) / realNumber(sensorproperties.maxWidth), u, v;
    u = ((static_cast<realNumber>(width) / static_cast<realNumber>(sensorproperties.maxWidth)) - 0.5);
    v = ((static_cast<realNumber>((sensorproperties.maxHeight-height)) / static_cast<realNumber>(sensorproperties.maxHeight)) - 0.5);
    
    if(ratio < 1.0)
    {
        v = v * ratio;
    }
    else
    {
        u = u / ratio;
    }

    //For each Supersample, sends multiple rays per pixel with random offsets for Anti-Aliasing
    for(int j = 0; j < SUPERSAMPLES; j++)
    {
        //1.1 Randomize UV space
        realNumber pixelJitterScale = 0.00025, rU, rV;
        rU = u + (nrandom() -0.5) * 2.0 * pixelJitterScale;
        rV = v + (nrandom() -0.5) * 2.0 * pixelJitterScale;

        //2. Scale UV´s to sensor
        fragmentCoords.set(-rU,rV,0.0);
        fragmentCoords = fragmentCoords * 0.001; //Scale sensor to mm range
        fragmentCoords = fragmentCoords * sensorproperties.SensorSize; //Scale to fit common sensor size

        //3. UV´s to ray origin and ray direction
        rayOrig.set(0.0,0.0,sensorproperties.FocalLength);
        rayOrig = rayOrig + sensorproperties.camPos;
        rayDir = (fragmentCoords + sensorproperties.camPos) - rayOrig;
        rayDir = rayDir / VectorLength(rayDir);

        //Orient camera
        M = rotateRay(sensorproperties.AngleX, sensorproperties.AngleY, sensorproperties.AngleZ);
        rayDir = M*rayDir;

        //Accumulate results for Debug and Magik render
        if((RENDER_DEBUG || RENDER_MAGIK) && !RENDER_MAGIK_SANDBOX)
        {
            mRender = mRender + traceRay_twopass(globalTime,rayOrig,rayDir,sensorproperties.camPos,sensorproperties.camMomentum,disk_lattice_noise_points,num_lattice_points,FRAME_jet_lattice_noise_points,RAY_jet_lattice_noise_points);
        }

        //Accumlate results for Sandbox
        if(RENDER_MAGIK_SANDBOX)
        {
            mRender = mRender + traceRay_sandbox(rayOrig,rayDir);
        }
    }

    //"Normalize" colors and return
    mRender = mRender / realNumber(SUPERSAMPLES);
    pixel.R = mRender.x;
    pixel.G = mRender.y;
    pixel.B = mRender.z;
    return pixel;
}


static void thread_trace_ray(const realNumber globalTime, const SensorProperties sensorproperties, noise_point *disk_lattice_noise_points, const int num_lattice_points, lattice_noise_point *FRAME_jet_lattice_noise_points, int* idx_ray, std::mutex* thread_lock, pix_realNumber *pixels_raw){
    int x,y;
    int total_rays = sensorproperties.maxWidth*sensorproperties.maxHeight;
    dynamic_lattice_point *RAY_jet_lattice_noise_points = (dynamic_lattice_point*)malloc(sizeof(dynamic_lattice_point)*num_lattice_points);
    while(idx_ray[0]<total_rays){
        thread_lock[0].lock();
        x = idx_ray[0] % sensorproperties.maxWidth;
        y = idx_ray[0] / sensorproperties.maxWidth;
        ++idx_ray[0];
        thread_lock[0].unlock();
        pixels_raw[x+y*sensorproperties.maxWidth]=pixel_function(globalTime,x,y,sensorproperties,disk_lattice_noise_points,num_lattice_points,FRAME_jet_lattice_noise_points,RAY_jet_lattice_noise_points);
    }
    free(RAY_jet_lattice_noise_points);
}


void trace_all_rays(const realNumber globalTime, const int total_threads, const SensorProperties sensorproperties, noise_point *disk_lattice_noise_points, const int num_lattice_points, lattice_noise_point *FRAME_jet_lattice_noise_points, pix_realNumber *pixels_raw){
    int idx_ray = 0;
    std::thread *ThreadPool = new std::thread[total_threads];
    std::mutex thread_lock;

    //Write to each pixel (Array Index)
    for(int i = 0; i<total_threads; ++i){
        ThreadPool[i] = std::thread(thread_trace_ray,globalTime,sensorproperties,disk_lattice_noise_points,num_lattice_points,FRAME_jet_lattice_noise_points, &idx_ray, &thread_lock, pixels_raw);
    }

    for(int i = 0; i<total_threads; ++i){
        ThreadPool[i].join();
    }

    delete [] ThreadPool;
}