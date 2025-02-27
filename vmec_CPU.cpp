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
static Vector3 evaluateDiskLatticePoint(const realNumber globalTime, const int i, const Vector3 &sPos, noise_point* disk_lattice_noise_points, const realNumber &t, const realNumber &r, const realNumber &theta)
{
    Vector4 diskLatticePoint;
    diskLatticePoint.x = disk_lattice_noise_points[i].r;
    diskLatticePoint.y = disk_lattice_noise_points[i].t;
    diskLatticePoint.z = disk_lattice_noise_points[i].p;
    diskLatticePoint.t = disk_lattice_noise_points[i].v;

    Vector3 uVelPoint, xPos;
    realNumber gamma, angularDisplacement, dr, dtheta, dphi, sampleDistance;        
    dr = diskLatticePoint.x * settings.disk_lattice_major_radius + settings.disk_lattice_minimum_radius;
            
    uVelPoint.set(0.0,0.0,uPhi(dr));
    gamma = generalizedTimeDilationFactor(r, theta, uVelPoint, -1.0);
    angularDisplacement = (globalTime - (t/gamma)) * -uVelPoint.z * sign(settings.a);        
    dtheta = (diskLatticePoint.y * settings.disk_lattice_vertical_scale) + (M_PI/2.0);
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
    initialAngularDisplacement = std::sqrt(1.0 / ((r+3.0) * (r+3.0))) * -sign(settings.a) * 156.0;
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
            angularDisplacement = ((globalTime*animationScale) - (t/gamma)) * -uVelPoints.z * sign(settings.a);

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
    distY = std::sqrt(settings.vr*(std::abs(rayOrig.y)/settings.vy) + 4.0); //Min y distance
    factor = (1.0 / std::sqrt( rayOrig.x*rayOrig.x + rayOrig.z*rayOrig.z )) * 0.25;
    sPos = rayOrig;
    sPos = rotate2D(sPos,settings.w*factor);

    for (int o = 0; o < settings.jet_noise_octaves; o++)
    {
        mWeight = 0, mValue = 0;
        octaveIndices = settings.jet_noise_first_octave*std::pow(2,o);
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
        realNumber r0 = (FRAME_jet_lattice_noise_points[k].r * settings.jet_radius_scale) + settings.jet_minor_radius;

        //Wrapping
        realNumber Iter = std::trunc(vtau / settings.max_t);
        realNumber truncTime = settings.max_t*Iter;
        vtau -= truncTime;

        //Position and Velocity
        realNumber ftau = std::sqrt(settings.vr*vtau + r0*r0);
        realNumber rtauSqrtTerm1 = std::sqrt( ((settings.vy*settings.vy*vtau*vtau - settings.a*settings.a + ftau*ftau)*(settings.vy*settings.vy*vtau*vtau - settings.a*settings.a + ftau*ftau)) + 4.0*settings.vy*settings.vy*settings.a*settings.a*vtau*vtau );
        realNumber rtauSqrtTerm2 = std::sqrt( settings.vy*settings.vy*vtau*vtau - settings.a*settings.a + ftau*ftau +  rtauSqrtTerm1 );

        //Position
        realNumber rtau = (1.0/std::sqrt(2.0)) * rtauSqrtTerm2;
        realNumber thetatau = std::acos((settings.vy*vtau/rtau*sign(verticalDirection)));

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


static int interp_find_index(realNumber target, const realNumber *s, const int max_steps, const int prev_idx)
{
    int new_idx = std::max(std::min(prev_idx,max_steps-1),0);
    if(s[new_idx]<=target && s[new_idx+1]>target){
    }
    else if(s[new_idx+1]<=target){
        while(new_idx<max_steps-1 && s[new_idx+1]<=target){
            ++new_idx;
        }
    }
    else{
        while(new_idx>0 && s[new_idx]>target){
            --new_idx;
        }
    }
    return new_idx;
}


// Linear interpolation
static RayProperties interp_along_geodesic(realNumber target, const RayProperties *ray_points, const int num_k, const realNumber *k, int &interp_idx)
{
    RayProperties result;

    interp_idx = interp_find_index(target, k, num_k, interp_idx);
    realNumber coeff = (target - k[interp_idx])/(k[interp_idx+1]-k[interp_idx]); 

    result = ray_points[interp_idx] + coeff*(ray_points[interp_idx+1]-ray_points[interp_idx]);
    return result;
}


//Trace Geodesic
static void trace_geodesic(RayProperties* geodesic, int &geodesic_size, realNumber* geodesic_k, localProperties &positionProperties, const Vector3 &rayOrig, const Vector3 &rayDir, const Vector3 &rayMomentum, bool &ray_convergence, bool &hit_disk, bool &hit_jet, bool &hit_ambient, int &i_rs, int &i_celestial_sphere, int &i_hit_disk, int &i_hit_jet, int &i_hit_ambient) 
{
    //Variable declaration
    bool 
        break_tag = false;
    
    int
        i;

    realNumber
        hStep = 1e-3;

    localProperties
        pointProperties = initializeCamera(rayOrig,rayMomentum);
        positionProperties = pointProperties;

    RayProperties
        pathProperties = initializeRay(rayOrig,rayDir,pointProperties);
        pathProperties.hStep = hStep;

    BVObjects 
        BoundingVolumes;

    realNumber
        k_current = 0.0_real;

    geodesic_size = 0;

    //Geodesic Integration
    for(i = 0; i < settings.integration_depth; i++)
    {
           
        if(pathProperties.r < settings.rs*settings.rs_scale) //Event Horizon
        {
            i_rs = i;
            break_tag = true;
        }
        else if(pathProperties.r > settings.celestial_sphere_radius) //Celestial Sphere
        {
            i_celestial_sphere = i;
            break_tag = true;
        }
        else if(i+1 == settings.integration_depth) //Ray did not converge
        {
            ray_convergence = false;
            break_tag = true;
        }

        //In-bounding region check
        pointInBounds(BoundingVolumes,pathProperties.r,pathProperties.theta);

        if(BoundingVolumes.disk && !hit_disk)
        {
            i_hit_disk = i;
            hit_disk = true;
        }
        if(BoundingVolumes.jet && !hit_jet)
        {
            i_hit_jet = i;
            hit_jet = true;
        }
        if(BoundingVolumes.ambient && !hit_ambient)
        {
            i_hit_ambient = i;
            hit_ambient = true;
        }

        //Update Array
        geodesic[i] = pathProperties;
        geodesic_k[i] = k_current;
        ++geodesic_size;

        //March step
        pathProperties = rayRKF45(pathProperties);
        k_current += euclidianDistance(BoyerLindquistToCartesian(geodesic[i].r,geodesic[i].theta,geodesic[i].phi),BoyerLindquistToCartesian(pathProperties.r,pathProperties.theta,pathProperties.phi));

        if(break_tag)
        {
            break;
        }
    }
}


//Debugging viewer
static Vector4 vmec_debugger(const realNumber globalTime, const Vector3& rayOrig, const Vector3& rayDir, const Vector3& camPos, const Vector3& camMomentum, RayProperties *geodesic, int &geodesic_size, realNumber* geodesic_k, noise_point* disk_lattice_noise_points, int num_lattice_points, lattice_noise_point* FRAME_jet_lattice_noise_points, dynamic_lattice_point* RAY_jet_lattice_noise_points)
{
    if(settings.test_scene && settings.lattice_scene)
    {
        //return 0; //Error
    }

    //Precompute Geodesic
        //Geodesic tags and bools
        bool 
            ray_convergence = true, 
            hit_disk = false,
            hit_jet = false,
            hit_ambient = false;

        int
            i_rs = settings.integration_depth - 1,
            i_celestial_sphere = settings.integration_depth - 1,
            i_hit_disk = settings.integration_depth - 1,
            i_hit_jet = settings.integration_depth - 1,
            i_hit_ambient = settings.integration_depth - 1;

        localProperties
            camProperties;

        //Compute Geodesic
        trace_geodesic(geodesic,geodesic_size,geodesic_k,camProperties,rayOrig,rayDir,camMomentum,ray_convergence,hit_disk,hit_jet,hit_ambient,i_rs,i_celestial_sphere,i_hit_disk,i_hit_jet,i_hit_ambient);

    //Debugger Variables
        bool 
            hit_disk_lattice = false,
            hit_jet_lattice = false,
            break_tag = false;
    
        int
            i = 0,
            hit_disk_ID = 0,
            hit_jet_ID = 0,
            idx_interp = 0,
            interp_idx = 0;

        realNumber
            hStep = 1e-3,
            dx = 1.0,
            checker = 0.0,
            redshift_factor = 0.0,
            k_intersect = 0.0,
            geodesic_length = 0.0,
            k = 0.0,
            dk = 0.5;

        Vector3
            dir(rayDir.x,rayDir.y,rayDir.z), 
            p1(rayOrig.x,rayOrig.y,rayOrig.z),
            p0(0.0,0.0,0.0),
            backgroundCol(0.0,0.0,0.0),
            checkerDiskNormal(0.0,1.0,0.0),
            checkerDiskp0(0.0,0.0,0.0),
            pointPosition(0.0,0.0,0.0),
            hitDisk(0.0,0.0,0.0),
            localVelocity(0.0,0.0,0.0);
    
        Vector4
            initialU(geodesic[0].u0,geodesic[0].u1,geodesic[0].u2,geodesic[0].u3), 
            RGBA_out(0.0,0.0,0.0,0.0);

        RayProperties
            rayProperties;

    //Checker Disk test scene
    if(settings.test_scene)
    {
        for(i = 0; i < geodesic_size; i++)
        {
            rayProperties = geodesic[i];

            //Event Horizon
            if(i == i_rs && !break_tag)
            {
                break_tag = true;
            }

            //Celestial Sphere
            if(i == i_celestial_sphere && !break_tag)
            {
                if(settings.colored_background)
                {
                    p1 = p1 / VectorLength(p1);
                    backgroundCol.set((p1.x + 1.0) / 2.0, (p1.y + 1.0) / 2.0, (p1.z + 1.0) / 2.0);
                }

                RGBA_out.set(RGBA_out.x + backgroundCol.x , RGBA_out.y + backgroundCol.y , RGBA_out.z + backgroundCol.z , 0.0);
            }

            //Ray did not hit anything
            if(!ray_convergence && !break_tag)
            {
                RGBA_out.set(1.0,1.0,1.0,1.0);
                break_tag = true;
            }

            //Checker Disk
            if(rayDiskIntersection(checkerDiskNormal,checkerDiskp0,p0,dir,dx,hitDisk) && settings.checker_disk_toggle && !break_tag)
            {
                k_intersect = geodesic_k[i-1] + euclidianDistance(p0,hitDisk);
                
                if(settings.checker_disk_Bi_sexual_color)
                {
                    p1 = hitDisk;
                    hitDisk = cartesianToBL(hitDisk);
                    checker = radialChecker(hitDisk.x, hitDisk.z, settings.checker_disk_major_radius, 0.1, settings.checker_disk_minor_radius);
                    p1 = p1 / VectorLength(p1);
                    RGBA_out.set(RGBA_out.x + checker*((p1.x + 1.0) / 2.0) , RGBA_out.y + checker*p1.y , RGBA_out.z + checker*p1.z, 0.0);
                }
                else if(settings.checker_disk_Z_factor)
                {
                    int idx_interp = 0;
                    rayProperties = interp_along_geodesic(k_intersect, geodesic, geodesic_size, geodesic_k, idx_interp);
                    hitDisk = cartesianToBL(hitDisk);
                    localVelocity.set(0.0,0.0,(uPhi(hitDisk.x)));
                    redshift_factor = generalizedRedshiftFactor(rayProperties,camProperties,initialU,localVelocity);
                    RGBA_out.set(redshift_factor,redshift_factor,redshift_factor,0.0);
                }

                break_tag = true;
            }

            //Break Condition
            if(break_tag)
            {
                if(settings.integration_depth_overlay) //Absolute Visualization
                {
                    RGBA_out.set(i,i,i,0.0);
                }

                if(settings.coordinate_overlay) //Absolute Visualization
                {
                    RGBA_out.set(rayProperties.r,rayProperties.theta,rayProperties.phi,0.0);
                }

                break;
            }

            //Calculate dx
            p0 = p1;
            p1 = BoyerLindquistToCartesian(rayProperties.r,rayProperties.theta,rayProperties.phi);
            dx = euclidianDistance(p0,p1);
            dir = p1 - p0;
            dir = dir / VectorLength(dir);
        }
    }

    //Lattice Point test scene
    if(settings.lattice_scene)
    {
        if(hit_jet)
        {
            //Create Per ray copy of FRAME array
            initialize_per_RAY_lattice_points(num_lattice_points, RAY_jet_lattice_noise_points, FRAME_jet_lattice_noise_points);
        }

        if(hit_disk && settings.bv_disk_toggle)
        {
            RGBA_out.x += 0.05;
        }

        if(hit_jet && settings.bv_jet_toggle)
        {
            RGBA_out.y += 0.05;
        }

        if(hit_ambient && settings.bv_ambient_medium_toggle)
        {
            RGBA_out.z += 0.05;
        }

        geodesic_length = geodesic_k[geodesic_size-1];

        while(geodesic_length > k)
        {
            //Lerp along stored geodesic and derive iterator
            rayProperties = interp_along_geodesic(k, geodesic, geodesic_size, geodesic_k, idx_interp);
            interp_idx = interp_find_index(k, geodesic_k, geodesic_size, interp_idx);
            i = interp_idx;

            //Event Horizon
            if(i >= i_rs && !break_tag)
            {
                break_tag = true;
            }

            //Celestial Sphere
            if(i >= i_celestial_sphere && !break_tag)
            {
                if(settings.colored_background)
                {
                    p1 = p1 / VectorLength(p1);
                    backgroundCol.set((p1.x + 1.0) / 2.0, (p1.y + 1.0) / 2.0, (p1.z + 1.0) / 2.0);
                }

                RGBA_out.set(RGBA_out.x + backgroundCol.x , RGBA_out.y + backgroundCol.y , RGBA_out.z + backgroundCol.z , 0.0);
            }

            //Ray did not hit anything
            if(!ray_convergence && !break_tag)
            {
                RGBA_out.set(1.0,1.0,1.0,1.0);
                break_tag = true;
            }

            //Visualize Disk Points
            if(settings.toggle_disk_lattice_points && hit_disk && i >= i_hit_disk && !break_tag)
            {
                for(int g = 0; g < settings.disk_n_points; g++ && !hit_disk_lattice)
                {
                    pointPosition = evaluateDiskLatticePoint(globalTime, g, p1, disk_lattice_noise_points, rayProperties.t, rayProperties.r, rayProperties.theta);
                    if(euclidianDistance(pointPosition,BoyerLindquistToCartesian(rayProperties.r, rayProperties.theta, rayProperties.phi)) < settings.disk_lattice_point_size)
                    {
                        hit_disk_lattice = true;
                        hit_disk_ID = g;
                        break_tag = true;
                    }
                }

                if(hit_disk_lattice)
                {
                    RGBA_out.x += ndetermin(hit_disk_ID+275);
                    RGBA_out.y += ndetermin(hit_disk_ID+4653);
                    RGBA_out.z += ndetermin(hit_disk_ID+57542);
                    RGBA_out.t = 0.0;
                    RGBA_out = RGBA_out/VectorLength(RGBA_out);
                }
            }

            //Visualize Jet Points
            if(settings.toggle_jet_lattice_points && hit_jet && i >= i_hit_jet && !break_tag)
            {
                for(int g = 0; g < settings.jet_n_points; g++)
                {
                    if(euclidianDistance(p1,BoyerLindquistToCartesian(RAY_jet_lattice_noise_points[g].r,RAY_jet_lattice_noise_points[g].theta,RAY_jet_lattice_noise_points[g].phi)) < settings.jet_lattice_point_size)
                    {
                        hit_jet_lattice = true;
                        hit_jet_ID = g;
                        break_tag = true;
                    }
                }

                if(hit_jet_lattice)
                {
                    RGBA_out.x += ndetermin(hit_jet_ID+5647);
                    RGBA_out.y += ndetermin(hit_jet_ID+68934);
                    RGBA_out.z += ndetermin(hit_jet_ID+923845);
                    RGBA_out.t = 0.0;
                    RGBA_out = RGBA_out/VectorLength(RGBA_out);
                }
            }

            //Break Condition & Absolute Visualizations
            if(break_tag)
            {
                if(settings.integration_depth_overlay)
                {
                    RGBA_out.set(i,i,i,0.0);
                }

                if(settings.coordinate_overlay)
                {
                    RGBA_out.set(rayProperties.r,rayProperties.theta,rayProperties.phi,0.0);
                }

                break;
            }

            //Calculate dx
            p0 = p1;
            p1 = BoyerLindquistToCartesian(rayProperties.r,rayProperties.theta,rayProperties.phi);
            dx = euclidianDistance(p0,p1);
            dir = p1 - p0;
            dir = dir / VectorLength(dir);

            //Parameter along geodesic
            k += dk;
        }
    }
    return RGBA_out;
}


//Magik Sandbox
static Vector4 vmec_magik_sandbox(const Vector3 rayOrig, const Vector3 rayDir)
{
    Vector4
        RGBA_out(0.0,0.0,0.0,0.0);
    
    if(settings.magik_trivial_scene)
    {
        Magik::trivial(RGBA_out,rayOrig,rayDir);
    }
    else if(settings.magik_single_scene)
    {
        Magik::single(RGBA_out,rayOrig,rayDir);
    }
    else if(settings.magik_single_MFP_scene)
    {
        Magik::single_MFP(RGBA_out,rayOrig,rayDir);
    }
    else if(settings.magik_muliple_scene)
    {
        Magik::multiple();
    }

    return RGBA_out;
}


//Main Pixel Shader
static pix_realNumber pixel_function(const realNumber globalTime, int width, int height, const SensorProperties sensorproperties, RayProperties* geodesic, int &geodesic_size, realNumber* geodesic_k, noise_point *disk_lattice_noise_points, int num_lattice_points, lattice_noise_point *FRAME_jet_lattice_noise_points,  dynamic_lattice_point* RAY_jet_lattice_noise_points)
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
    for(int j = 0; j < settings.supersamples; j++)
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

        //Accumulate results for Debugger
        if(!settings.render_magik_sandbox && settings.render_debug)
        {
            mRender = mRender + vmec_debugger(globalTime, rayOrig, rayDir, sensorproperties.camPos, sensorproperties.camMomentum, geodesic, geodesic_size, geodesic_k, disk_lattice_noise_points, num_lattice_points, FRAME_jet_lattice_noise_points, RAY_jet_lattice_noise_points);
        }

        //Accumlate results for Sandbox
        if(settings.render_magik_sandbox)
        {
            mRender = mRender + vmec_magik_sandbox(rayOrig,rayDir);
        }
    }

    //"Normalize" colors and return
    mRender = mRender / realNumber(settings.supersamples);
    pixel.R = mRender.x;
    pixel.G = mRender.y;
    pixel.B = mRender.z;
    return pixel;
}


static void thread_trace_ray(int thread_idx, const realNumber globalTime, const SensorProperties sensorproperties, RayProperties* all_geodesics, int* geodesic_sizes, realNumber* geodesics_k, noise_point *disk_lattice_noise_points, const int num_lattice_points, lattice_noise_point *FRAME_jet_lattice_noise_points, dynamic_lattice_point *dynamic_noise_points, int* idx_ray, std::mutex* thread_lock, pix_realNumber *pixels_raw){
    int x,y;
    int total_rays = sensorproperties.maxWidth*sensorproperties.maxHeight;
    RayProperties* geodesic = &all_geodesics[thread_idx * settings.integration_depth];
    int geodesic_size = geodesic_sizes[thread_idx];
    realNumber* geodesic_k = &geodesics_k[thread_idx * settings.integration_depth];
    dynamic_lattice_point *RAY_jet_lattice_noise_points = &dynamic_noise_points[num_lattice_points*thread_idx];
    while(idx_ray[0]<total_rays){
        thread_lock[0].lock();
        x = idx_ray[0] % sensorproperties.maxWidth;
        y = idx_ray[0] / sensorproperties.maxWidth;
        ++idx_ray[0];
        thread_lock[0].unlock();
        pixels_raw[x+y*sensorproperties.maxWidth]=pixel_function(globalTime, x, y, sensorproperties, geodesic, geodesic_size, geodesic_k, disk_lattice_noise_points, num_lattice_points, FRAME_jet_lattice_noise_points, RAY_jet_lattice_noise_points);
    }
}


void trace_all_rays(const realNumber globalTime, const int total_threads, const SensorProperties sensorproperties, RayProperties* all_geodesics, int* geodesic_sizes, realNumber* geodesics_k, noise_point *disk_lattice_noise_points, const int num_lattice_points, lattice_noise_point *FRAME_jet_lattice_noise_points, dynamic_lattice_point *dynamic_noise_points, pix_realNumber *pixels_raw){
    int idx_ray = 0;
    std::thread *ThreadPool = new std::thread[total_threads];
    std::mutex thread_lock;

    //Write to each pixel (Array Index)
    for(int i = 0; i<total_threads; ++i){
        ThreadPool[i] = std::thread(thread_trace_ray, i, globalTime,sensorproperties, all_geodesics, geodesic_sizes, geodesics_k, disk_lattice_noise_points, num_lattice_points, FRAME_jet_lattice_noise_points, dynamic_noise_points, &idx_ray, &thread_lock, pixels_raw);
    }

    for(int i = 0; i<total_threads; ++i){
        ThreadPool[i].join();
    }

    delete [] ThreadPool;
}