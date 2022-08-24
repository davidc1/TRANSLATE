#include <array>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>

#include "Electron.h"

/*
 * Converts electronvolts to joules
 * 
 * @param eV Energy (in eV)
 * @return Energy (in J)
 */
inline double eV_to_J(double eV)
{
    return eV * 1.60218e-19;
}

/*
 * Converts joules to electronvolts
 * 
 * @param J Energy (in J)
 * @return Energy (in eV)
 */
inline double J_to_eV(double J)
{
    return J * 6.242e18;
}

/*
 * Given an electron's kinetic energy, returns its velocity (classically)
 * 
 * @param eV The electron's kinetic energy (in eV)
 * @return The velocity of the electron (in m/s)
 */
inline double eV_to_v(double eV)
{
    double J = eV_to_J(eV);
    return sqrt(2 * J / m_e);
}

/*
 * Produces a random unit vector
 * 
 * @param gen The random number generator to be used
 * @return Random unit vector
 */
Vec random_unit_vector(std::mt19937& gen)
{
    std::normal_distribution<double> dist(0.0, 1.0);

    Vec v(dist(gen), dist(gen), dist(gen));

    if (norm(v) != 0)
        v /= norm(v);

    return v;
}

/*
 * Return a velocity whose magnitude is chosen from the Maxwell-Boltzmann distribution for
 * argon / the electron. The velocity's direction is randomly chosen
 * 
 * @param gen The random number generator to be used
 * @param electron Whether the velocity should be chosen from the electron distribution. If false, it draws from the argon distribution
 * @return The velocity vector (in m/s)
 */
Vec random_velocity(std::mt19937& gen, bool electron=false)
{
    if (!interactions)
        return Vec(0, 0, 0);

    double a;

    if (!electron)
        a = sqrt(k_b * T / M_a);
    else
        a = sqrt(k_b * T / m_e);

    std::gamma_distribution<double> maxwell(1.5, 1 / (2 * a * a));
    Vec v = random_unit_vector(gen) * sqrt(maxwell(gen));

    return v;
}

/*
 * Chooses the starting position of simulated electron randomly from a disc at the
 * top of the simulated volume. If running the simulation in a uniform field, this returns
 * the origin.
 * 
 * @param gen The random number generator to be used
 * @return The position vector (in m)
 */
Vec starting_pos(std::mt19937& gen)
{
    if (uniform_field)
        return Vec(0, 0, 0);

    std::normal_distribution<double> dist(0.0, 1.0);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    Vec random_disc(0, dist(gen), dist(gen));

    if (norm(random_disc) != 0)
        random_disc /= norm(random_disc);

    random_disc *= uniform(gen);
    random_disc *= (r_max - r_min) * spawn_width_scale;

    return Vec((z_max - z_min) * spawn_height_scale + z_min, random_disc.y, random_disc.z);
}

/*
 * Given a parallel coordinate to the tip array plane (y or z), this shifts it to 
 * within the bounds of a single tip volume. Effectively, this tiles the simulated volume
 * like a chessboard.
 * 
 * @param x The coordinate to be shifted (in m)
 * @return The shifted coordinate (in m)
 */
double shift_coord(double x)
{
    bool negative = std::signbit(x);
    double width = r_max - r_min;
    
    if (negative) {
        while ((x + 2 * width) < r_max)
            x += 2 * width;
    } else {
        while ((x - 2 * width) > -r_max)
            x -= 2 * width;
    }

    return x;
}

/*
 * Shifts the y and z coordinates of a position vector to lie within the bounds of
 * a single tip, for the purposes of calculating the force. This function is used
 * only when simulating a tip array.
 * 
 * @param pos The position vector to be altered (in m)
 * @return The shifted position vector
 */
Vec local_coords(Vec pos)
{
    return Vec(pos.x, shift_coord(pos.y), shift_coord(pos.z));
}

/*
 * Gives the radial coordinate of the simulated electron. If the electron lies outside
 * the simulated volume, its radial coordinate is chosen to be the maximum of the volume.
 * 
 * @param local_pos The position vector of the electron relative to the nearest tip (in the single tip case, this is the origin) (in m)
 * @return The electron's radial coordinate
 */
double r(Vec local_pos)
{
    double r_val = sqrt(local_pos.y * local_pos.y + local_pos.z * local_pos.z);

    if (r_val > r_max)
        r_val = r_max;

    return r_val;
}

/*
 * Gives the height cordinate of the simulated electron. This is just x, so this function
 * effectively acts like a clamp.
 * 
 * @param local_pos The position vector of the electron relative to the nearest tip (in m) (in the single tip case, this is the origin)
 * @return The electron's height coordinate
 */
double z(Vec local_pos)
{
    double z_val = local_pos.x;

    if (z_val > z_max)
        z_val = z_max;

    if (z_val < 0)
        z_val = 0;

    return z_val;
}

/*
 * Given a position vector, returns whether it lies within or on the surface of the
 * simulated tip (or tip array).
 * 
 * @param pos The position vector the check (in m)
 * @return Whether the vector touches the tip
 */
bool hit_check(Vec pos)
{
    Vec local = local_coords(pos);

    double r_val;
    double z_val;
    
    if (!single_tip) {
        r_val = r(local);
        z_val = z(local);
    } else {
        r_val = r(pos);
        z_val = z(pos);
    }

    if (z_val <= z_min) 
        return true;

    // For rounded tips
    if (rounded_tip) {
        // For cylindrical sides
        if (top_radius == bot_radius) {
            if ( ((z_val - top_z) * (z_val - top_z) + r_val * r_val <= top_radius * top_radius) || (z_val <= top_z && r_val <= top_radius) )
                return true;
        // For tapered sides
        } else {
            // For cylindrical sides after a tapered section
            if (z_val < bot_z) {
                if (r_val < bot_radius)
                    return true;
                else
                    return false;
            }
            if ( ((z_val - top_z) * (z_val - top_z) + r_val * r_val <= top_radius * top_radius) || (z_val <= top_z && (z_val <= (bot_z - top_z) / (bot_radius - top_radius) * (r_val - bot_radius) + bot_z)) )
                return true;
        }
    // For flat tips
    } else {
        // For cylindrical sides
        if (top_radius == bot_radius) {
            if (z_val <= top_z && r_val <= top_radius)
                return true;
        // For tapered sides
        } else {
            if (z_val <= top_z && (z_val <= (bot_z - top_z) / (bot_radius - top_radius) * (r_val - bot_radius) + bot_z))
                return true;
        }
    }

    return false;
}

/*
 * Gives the radial index to access the electric field lookup table.
 * 
 * @param r_val The radial coordinate (in m)
 * @return The associated index
 */
int r_index(double r_val)
{
    return floor(r_conversion_factor * (r_val - r_min));
}


/*
 * Gives the height index to access the electric field lookup table.
 * 
 * @param z_val The height coordinate (in m)
 * @return The associated index
 */
int z_index(double z_val)
{
    return floor(z_conversion_factor * (z_val - z_min));
}

/*
 * Gives the acceleration felt by the electron in a non-uniform field.
 * 
 * @param pos The position vector of the electron (in m)
 * @param volts The bulk field (in V/cm)
 * @return The acceleration of the electron (in m/s^2)
 */
Vec accel_from_E(Vec pos, double volts)
{
    Vec local = local_coords(pos);

    double r_val;
    double z_val;
    
    if (!single_tip) {
        r_val = r(local);
        z_val = z(local);
    } else {
        r_val = r(pos);
        z_val = z(pos);
    }

    int i = r_index(r_val);
    int j = z_index(z_val);
    
    double E_r_val = E_r[num_z_steps * i + j];
    double E_z_val = E_z[num_z_steps * i + j];

    Vec inward_vec;
    
    if (!single_tip)
        inward_vec = Vec(0, -local.y, -local.z);
    else
        inward_vec = Vec(0, -pos.y, -pos.z);

    if (norm(inward_vec) != 0)
        inward_vec /= norm(inward_vec);

    inward_vec *= E_r_val;

    Vec e_field(-E_z_val, inward_vec.y, inward_vec.z);

    return (volts / bulk_field) * (e / m_e) * e_field * 1e2;
}

/*
 * The constructor for the electron class. Sets up the basic initial conditions.
 * 
 * @param ions A parameter passed to child electrons to keep track of the total number of ionizations
 * @param volts The bulk field (in V/cm)
 * @param position The initial position vector (in m)
 * @param velocity The initial velocity vector (in m)
 * @param gen The random number generator to be used
 */
Electron::Electron(double initial_time, double volts, Vec position, Vec velocity, std::mt19937& gen, int debug):
  child_ions(0), x(position), v(velocity), accel((e / m_e) * volts * 1e2, 0, 0), volts_per_cm(volts), total_time(initial_time), generator(gen), debug(debug)
{
    if (!uniform_field)
        accel = accel_from_E(x, volts_per_cm);

    energy = J_to_eV(0.5 * m_e * dot(v, v));

    scatter = 0;
}


/*
 * Gives the index to be used when accessing cross-section lookup tables.
 * 
 * @param v The velocity of the electron (in m/s) (this is sometimes the velocity of the electron relative to a chosen argon atom)
 * @return The associated index
 */
int Electron::index(double v) const
{
    double E = J_to_eV(0.5 * m_e * v * v);

    if (E < eV_min)
        E = eV_min;

    int ind = round(((eV_steps - 1) / log10(eV_max / eV_min)) * (log10(E) - log10(eV_min)));
    assert(ind < 1000);

    if (debug) { std::cout << std::scientific << "[ DEBUG ] velocity " << v << " corresponds to energy " << E << std::endl; }

    return ind;
}


/*
 * Gives the index to be used when accessing cross-section lookup tables.
 * 
 * @param v The velocity of the electron (in m/s) (this is sometimes the velocity of the electron relative to a chosen argon atom)
 * @return The associated index
 */
int Electron::index_diff(double v) const
{
  double E = J_to_eV(0.5 * m_e * v * v);

  int ind = 0;
  
  if (E < eV_min_diff){
    E = eV_min_diff;
    ind = 0;
  }
  else if (E > eV_max_diff) {
    E = eV_max_diff;
    ind = eV_steps_diff - 1;
  }
  else {
    while (E <= energy_vals_diff[ind])
      ind++;
  }

  // account for binning taking energy and angle bins into account
  ind = 1 + ind * (angle_steps_diff + 2);

  if (debug) {
    std::cout << std::scientific;
    std::cout << "[ DEBUG ] Energy " << E << " -> index " << ind << std::endl;
    std::cout << "[ DEBUG ] final index " << ind << std::endl;
  }
  
  return ind;
}


/*
 * Removes an amount of energy from the electron, then updates its velocity to match.
 * 
 * @param eV The energy to be removed (in eV)
 */
void Electron::remove_energy(double eV)
{
    energy -= eV;

    if (energy < 0)
        energy = 0;
    
    v = sqrt((2 * eV_to_J(energy)) / m_e) * (v / norm(v));

    return;
}

/*
 * Removes the ionization energy of argon from the electron, then, if using the full
 * algorithm, spawns a new electron at the current one's position. The ionization counter is
 * increased by one.
 * 
 * @param electron_list The vector to which child electrons should be pushed
 * @param total_ionizations The global ionization counter to increase upon ionization
 */
void Electron::ionization(std::vector<Electron*> &electron_list, int& total_ionizations)
{
    if (track_child_ions) {
        std::uniform_real_distribution<double> rand_eV(1.0, 5.0);
        Vec near_therm_vel = random_unit_vector(generator) * eV_to_v(rand_eV(generator));
        electron_list.push_back(new Electron(total_time, volts_per_cm, x, near_therm_vel, generator,debug));
    }

    remove_energy(15.76);
    total_ionizations++;
    child_ions++;

    return;
}

/*
 * Updates the electron's velocity, assuming an elastic collision with an argon
 * atom.
 * 
 * @param u The magnitude of the electron's velocity relative to a chosen argon atom (in m/s)
 * @param vm The chosen argon atom's velocity vector (in m/s)
 */
void Electron::elastic_collision(double u, Vec vm)
{
    Vec unit_vec = random_unit_vector(generator);

    Vec a = (M_a * u) * unit_vec;
    Vec b = (m_e) * v;
    Vec c = (M_a) * vm;

    Vec v1 = (a + b + c) / (m_e + M_a);

    double r = x_sec_p(u) / std::max(x_sec_e(u), x_sec_p(u));

    std::uniform_real_distribution<double> rand_momentum_change(0.0, 1.0);

    if (rand_momentum_change(generator) < r)
        v = v1;
    else
        v = norm(v1) / norm(v) * v;

    return;
}


/*
 * Updates the electron's velocity, assuming an elastic collision with an argon
 * atom accounting for angle-depenence in xsec.
 * reference: http://homepages.engineering.auckland.ac.nz/~pkel015/SolidMechanicsBooks/Part_III/Chapter_1_Vectors_Tensors/Vectors_Tensors_05_Coordinate_Transformation_Vectors.pdf
 * reference: https://mathworld.wolfram.com/SphericalCoordinates.html
 * 
 * @param u The magnitude of the electron's velocity relative to a chosen argon atom (in m/s)
 * @param vm The chosen argon atom's velocity vector (in m/s)
 */
void Electron::elastic_collision_diff(double u, Vec vm)
{

  auto scatterinfo = x_sec_e_angle_diff(u);

  double theta = scatterinfo.second * 3.1415 / 180.; // scattering angle in radians

  scatter = scatterinfo.second;

  if (debug)
  std::cout << std::scientific << " [ DEBUG ] " << " simulated scattering angle in rad is " << theta << std::endl;

  // calculate cosine of scatter -> this is the component w.r.t. the
  // incoming electron's direction
  double uz = cos(theta);

  if (debug)
    std::cout << std::scientific << " [ DEBUG ] " << " simulated scattering angle in rad is " << theta << " with cos(theta) = " << uz << std::endl;

  // random unit vector
  
  // random phi angle to complete full determination of new direction vector
  std::uniform_real_distribution<double> uniform(0.0, 2 * 3.1415);
  double phi = uniform(generator);

  if (debug)
    std::cout << std::scientific << " [ DEBUG ] " << " simulated angle phi in rad is " << phi << std::endl;

  double ux = cos(phi) * sin(theta);
  double uy = sin(phi) * sin(theta);

  if (debug) {
    std::cout << std::scientific << " [ DEBUG ] " << " scatter angle     vector ( " << ux << ", " << uy << ", " << uz << ") with magnitude " << norm(v) << std::endl;
    std::cout << std::scientific << " [ DEBUG ] " << " original velocity vector ( " << v.x << ", " << v.y << ", " << v.z << ") with magnitude " << norm(v) << std::endl;
  }
  
  // compute "velocity" vector phi and theta
  double phi_vel = atan2(v.y,v.x);
  double theta_vel = acos(v.z / norm(v));
  // define coordinate vectors of phi_vel / theta_vel coordinate system
  // in the standard cartesian coordinate system
  //Vec r3(cos(phi_vel)*sin(theta_vel), sin(phi_vel)*sin(theta_vel), cos(theta_vel));
  //Vec r2(-sin(phi_vel),cos(phi_vel),0.);
  //Vec r1(cos(phi_vel)*cos(theta_vel),sin(phi_vel)*cos(theta_vel),-sin(theta_vel));

  // apply rotation
  //Vec e1(1,0,0);
  //Vec e2(0,1,0);
  //Vec e3(0,0,1);
  
  //double ux_rot = dot(r1,e1) * ux + dot(r1,e2) * uy + dot(r1,e3) * uz;
  //double uy_rot = dot(r2,e1) * ux + dot(r2,e2) * uy + dot(r2,e3) * uz;
  //double uz_rot = dot(r3,e1) * ux + dot(r3,e2) * uy + dot(r3,e3) * uz;

  double ux_rot = cos(phi_vel) * cos(theta_vel) * ux - sin(phi_vel) * uy + cos(phi_vel) * sin(theta_vel) * uz;
  double uy_rot = sin(phi_vel) * cos(theta_vel) * ux + cos(phi_vel) * uy + sin(phi_vel) * sin(theta_vel) * uz;
  double uz_rot = -sin(theta_vel) * ux + cos(theta_vel) * uz;

  Vec newdir(ux_rot,uy_rot,uz_rot);

  if (debug) 
    std::cout << std::scientific << " [ DEBUG ] " << " direction unit vector norm is " << norm(newdir) << std::endl;

  if (debug) {
    // dot product between old and new direction vector
    double dirdot = dot(newdir,v) / (norm(v) * norm(newdir));
    std::cout << std::scientific << " [ DEBUG ] scatter dot product is " << dirdot << " and simulated angle is " << cos(theta) << std::endl;
  }

  v = newdir * norm(v);

  if (debug)
    std::cout << std::scientific << " [ DEBUG ] " << " output velocity vector ( " << v.x << ", " << v.y << ", " << v.z << ") with magnitude " << norm(v) << std::endl;  
  
  return;
}


/*
 * Classically updates the electron's position and velocity.
 */
void Electron::update_pos_vel()
{
    Vec x_old = x;
    
    if (!uniform_field) {
        x += v * time_to_collision + 0.5 * accel_from_E(x, volts_per_cm) * time_to_collision * time_to_collision;
        v += 0.5 * (accel_from_E(x_old, volts_per_cm) + accel_from_E(x, volts_per_cm)) * time_to_collision;
    } else {
        x += v * time_to_collision + 0.5 * accel * time_to_collision * time_to_collision;
        v += accel * time_to_collision;
    }

    return;
}

/*
 * Draws the next timestep to use.
 * 
 * @param u The magnitude of the electron's velocity relative to a chosen argon atom (in m/s)
 * @return The drawn timestep (in s)
 */
double Electron::next_collision(double u)
{
    std::exponential_distribution<double> exponential(lambda);
    double time_step = exponential(generator);

    if (time_step > 3 * beta)
        time_step = 3 * beta;

    return time_step;
}

/*
 * Gives the probability for an interaction.
 * 
 * @param u The magnitude of the electron's velocity relative to a chosen argon atom (in m/s)
 * @param x_sec The value of the interaction's cross-section at the kinetic energy associated with the above relative velocity (in cm^2)
 * @return The associated interaction probability
 */
double Electron::probability(double u, double x_sec)
{
    double prob = (u * 1e2 * x_sec) / K_max;

    if (debug) {
      std::cout << std::scientific;
      std::cout << "[ DEBUG ] K max : " << K_max << std::endl;
      std::cout << "[ DEBUG ] Prob  : " << prob << std::endl;
      std::cout << "[ DEBUG ] reached a simulation time of " << total_time << std::endl;
    }
    
    assert (prob < 1);
    return prob;
}

/*
 * Advances the simulation by one step. In this order, this function...
 *      1. Draws a timestep
 *      2. Updates the parent electron's position and velocity
 *      3. Determines the probability of all interactions
 *      4. Simulates the interaction that passes the probability check
 *      5. Updates the electron's energy and total simulated time 
 * 
 * @param electron_list The vector to which child electrons should be pushed
 * @param total_ionizations The global ionization counter to increase upon ionization
 */
void Electron::update(std::vector<Electron*> &electron_list, int& total_ionizations)
{   
    // Draw a timestep and update the electron's position and velocity
    time_to_collision = next_collision(norm(v));
    update_pos_vel();

    scatter = -1;

    // If interactions are turned off, implement a terminal velocity
    if (!interactions) {
        if (norm(v) > 1000) {
            v /= norm(v);
            v *= 1000;
        }
    }

    // Draw an argon vector and determine interaction probabilities
    energy = J_to_eV(0.5 * m_e * dot(v, v));

    Vec vm = random_velocity(generator);
    double u = norm(v - vm);

    if (debug) { std::cout << "[ DEBUG ] elastic collision..." << std::endl; }
    
    double col_prob = probability(u, x_sec_e_diff(u));

    if (debug)
      std::cout << std::scientific << "[ DEBUG ] velocity "  << u << " and prob " << col_prob << std::endl;

    if (debug) { std::cout << "[ DEBUG ] ionization..." << std::endl; }
    double ion_prob = probability(u, x_sec_i(u));
    if (debug) { std::cout << "[ DEBUG ] excitation 11..." << std::endl; }
    double ex_11_prob = probability(u, x_sec_ex(u, 11));
    if (debug) { std::cout << "[ DEBUG ] excitation 13..." << std::endl; }
    double ex_13_prob = probability(u, x_sec_ex(u, 13));
    if (debug) { std::cout << "[ DEBUG ] excitation 14..." << std::endl; }
    double ex_14_prob = probability(u, x_sec_ex(u, 14));
    if (debug) { std::cout << "[ DEBUG ] excication 15..." << std::endl; }
    double ex_15_prob = probability(u, x_sec_ex(u, 15));

    assert(col_prob + ion_prob + ex_11_prob + ex_13_prob + ex_14_prob + ex_15_prob <= 1);

    // If interactions are turned off, set all probabilities to 0
    if (!interactions) {
        col_prob = 0;
        ion_prob = 0;
        ex_11_prob = 0;
        ex_13_prob = 0;
        ex_14_prob = 0;
        ex_15_prob = 0;
    }
    
    std::uniform_real_distribution<double> rand_event(0.0, 1.0);
    double prob = rand_event(generator);
    
    // Simulate the interaction that passes the check
    if (prob < ion_prob) {
        ionization(electron_list, total_ionizations);
    } else if (prob < (ion_prob + ex_11_prob)) {
        remove_energy(0.00);
    } else if (prob < (ion_prob + ex_11_prob + ex_13_prob)) {
        remove_energy(11.68);
    } else if (prob < (ion_prob + ex_11_prob + ex_13_prob + ex_14_prob)) {
        remove_energy(13.33);
    } else if (prob < (ion_prob + ex_11_prob + ex_13_prob + ex_14_prob + ex_15_prob)) {
        remove_energy(14.22);
    } else if (prob < (ion_prob + ex_11_prob + ex_13_prob + ex_14_prob + ex_15_prob + col_prob)) {
      //elastic_collision(u, vm);
        elastic_collision_diff(u, vm);
    }
    
    // Update the electron's energy and total time simulated
    energy = J_to_eV(0.5 * m_e * dot(v, v));
    total_time += time_to_collision;
}

/*
 * Generates a set of simulations, saving the files to py/studies/simulation-runs/
 * The files generated are plain text, with each row being a simulation step. The columns
 * are...
 * 
 * time (in ns), x-position (in um), y-position (in um), z-position (in um), energy (in eV), drift velocity (in m/s), number of ionizations
 * 
 * @param volts The bulk field (in V/cm)
 * @param cutoff The cutoff time (in s) (this is superseded by hit checks when using a non-uniform field)
 * @param cores The number of available processor cores (used in naming files)
 * @param write_every Save one simulated step for every write_every steps
 * @param k The core that the simulation runs on (used in naming files)
 * @param batches The number of batches to be generated
 * @param bar The progress bar to be used
 */
void generate_plot(int volts, double cutoff, int cores, int write_every, int k, int batches, int debug, ProgressBar& bar)
{
    for (int i = 0; i < batches; i++) {
        // Setup a progress bar and create the random number generator for the thread
        bar.new_file(cores * i + (k + 1), k);
        std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());

        // Setup the array of electrons to be simulated
        std::vector<Electron*> electron_list;
        electron_list.push_back(new Electron(0, volts, starting_pos(generator), random_velocity(generator, true), generator, debug));

        // Open the file to write to
        std::ofstream file("../py/simulation-runs/" + std::to_string(volts) + "V_" + std::to_string(cores * i + (k + 1)) + ".txt");
        assert(file.is_open());
        
        // Grab initial conditions
        double starting_x = 0;
        if (!uniform_field)
            starting_x = electron_list[0]->position().x;
        double t = electron_list.front()->elapsed_time();
        double x = electron_list.front()->position().x;
        double y = electron_list.front()->position().y;
        double z = electron_list.front()->position().z;
        double ke = electron_list.front()->ke();
        double vd = (x - starting_x) / t;
	double s = electron_list.front()->angle();
        int total_ionizations = 0;

        // Create a counter to check against when skipping steps to be written
        int simulation_step = 1;

        while (electron_list[0]->elapsed_time() < cutoff) {
            // Update electrons, adding child electrons if necessary
            std::vector<Electron*> new_electrons;
            for (auto elec : electron_list)
                elec->update(new_electrons, total_ionizations);
            electron_list.insert(electron_list.end(), new_electrons.begin(), new_electrons.end());

            // Get the relevant information of the primary electron (saved here in case this loop ends)
            t = electron_list.front()->elapsed_time();
            x = electron_list.front()->position().x;
            y = electron_list.front()->position().y;
            z = electron_list.front()->position().z;
            ke = electron_list.front()->ke();
            vd = (x - starting_x) / t;
	    s  = electron_list.front()->angle();

	    //if (s < 0) continue;

            // Remove electrons that have reached the anode
            if (!uniform_field){
	      electron_list.erase(std::remove_if(electron_list.begin(), electron_list.end(), [](decltype(*begin(electron_list)) ii){ return hit_check(ii->position()); }), electron_list.end());
	      //electron_list.erase(std::remove_if(electron_list.begin(), electron_list.end(), [](auto const& i){ return hit_check(i->position()); }), electron_list.end());
	    }

            // If all electrons have reached the anode, end the simulation
            if ((!uniform_field) && (electron_list.size() == 0))
                break;

            // Only save every write_every steps
            if (simulation_step++ % write_every != 0)
                continue;

            // Update the progress bar
            if (uniform_field)
                bar.update(electron_list.front()->elapsed_time() / cutoff, k);
            else
                bar.update(1 - (electron_list.front()->position().x - z_min) / ((z_max - z_min) * spawn_height_scale), k);

            // Display the progress bar
            if (bar.min_prog(k))
                bar.display();

            // Write the primary electron's information to a file (if using the full ionization algorithm, this may cause only the time and total ionizations to be accurate)
            file << t * 1e9 << "," << x * 1e6 << "," << y * 1e6 << "," << z * 1e6 << "," << ke << "," << vd << ","<< s << "," << total_ionizations << "\n";
        }

        // Update the progress bar one last time
        bar.update(1, k);

        // Write information to the file one last time (to catch all final ionizations)
        file << t * 1e9 << "," << x * 1e6 << "," << y * 1e6 << "," << z * 1e6 << "," << ke << "," << vd << ","<< s << "," << total_ionizations << "\n";
        file.close();
        assert(!file.is_open());
    }

    return;
}
