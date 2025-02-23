#include "vmec_setup.h"
#include "vmec_common.h"


//RNG Points
void generate_noise_points(int &num_disk_points, noise_point* &disk_noise_points){

	int octaveSize = 0;
	num_disk_points = 0;
	noise_point noise_point_curr;
	std::default_random_engine generator;
	std::uniform_real_distribution<realNumber> r_dist(0.0,1.0);
	std::uniform_real_distribution<realNumber> t_dist(-45.0*(M_PI/180.0),45.0*(M_PI/180.0));
	std::uniform_real_distribution<realNumber> p_dist(0.0,2.0*M_PI);
	std::uniform_real_distribution<realNumber> v_dist(-1.0,1.0);

	std::vector<noise_point> noise_points_vec;
	for(int o = 0; o < settings.disk_noise_octaves; o++)
	{
		octaveSize = settings.disk_noise_size_first_octave*std::pow(2,o);
		for(int idx = 0; idx < octaveSize; idx++){
			noise_point_curr.r = std::sqrt(r_dist(generator));
			noise_point_curr.t = t_dist(generator);
			noise_point_curr.p = p_dist(generator);
			noise_point_curr.v = v_dist(generator);
			noise_points_vec.push_back(noise_point_curr);
		}
		std::sort(noise_points_vec.begin() + num_disk_points, noise_points_vec.begin() + num_disk_points + octaveSize, [](const noise_point &a, const noise_point &b){return a.r < b.r;});
		num_disk_points+=octaveSize;
	}

	disk_noise_points = allocate<noise_point>(num_disk_points);
	std::copy(noise_points_vec.begin(),noise_points_vec.end(),disk_noise_points);

}


//RNG Points Jet Lattice
void generate_lattice_noise_points(int &num_lattice_points, lattice_noise_point* &lattice_noise_points){

	int octaveIndices = 0;
	num_lattice_points = 0;
	num_lattice_points = settings.jet_noise_first_octave * std::pow(2,settings.jet_noise_octaves-1);

   lattice_noise_points = allocate<lattice_noise_point>( num_lattice_points );

	std::default_random_engine generator;
	std::uniform_real_distribution<realNumber> r_dist(0.0,1.0);
	std::uniform_real_distribution<realNumber> p_dist(0.0,2.0*M_PI);
	std::uniform_real_distribution<realNumber> v_dist(-1.0,1.0);
	std::uniform_real_distribution<realNumber> W_dist(0,M_PI_2);
    std::uniform_real_distribution<realNumber> s_dist(1,-1);

	for(int idx=0;idx<num_lattice_points;++idx)
    {
		lattice_noise_points[idx].r = (std::sqrt(r_dist(generator)) * settings.jet_radius_scale) + settings.jet_minor_radius;
		lattice_noise_points[idx].p = p_dist(generator);
		lattice_noise_points[idx].t = M_PI_2;
		lattice_noise_points[idx].v = v_dist(generator);
		lattice_noise_points[idx].W = W_dist(generator);
        lattice_noise_points[idx].s = s_dist(generator);
	}
}


void do_setup(pix_realNumber* &pixels_raw, pix_RGB* &pixels_clean, int &num_lattice_points, lattice_noise_point* &lattice_noise_points, int &num_disk_points, noise_point* &disk_noise_points){

	pixels_raw = allocate<pix_realNumber>( settings.image_height * settings.image_width );
	pixels_clean = allocate<pix_RGB>( settings.image_height * settings.image_width );

	//Jet Lattice points
	generate_lattice_noise_points(num_lattice_points, lattice_noise_points);

	//noise points - sorted by r within each octave
	generate_noise_points(num_disk_points, disk_noise_points);

}


void do_cleanup(pix_realNumber* &pixels_raw, pix_RGB* &pixels_clean, lattice_noise_point* &lattice_noise_points, noise_point* &disk_noise_points){

	//Free Memory
	free(pixels_raw);
	free(pixels_clean);
	free(lattice_noise_points);
	free(disk_noise_points);

}
