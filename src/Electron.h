#ifndef ELECTRON_H
#define ELECTRON_H

#include <algorithm>
#include <assert.h>
#include <random>
#include <vector>
#include <utility>

#include "Constants.h"
#include "LUTs/CrossSections.h"
#include "LUTs/DiffCrossSections.h"
#include "LUTs/DiffCrossSectionsLAr.h"
#include "LUTs/Tip.h"
#include "ProgressBar.h"
#include "Vec.h"

class Electron
{
private:
  
  Vec _x, _v, _accel;
  Vec _lastinteractionposition;
  double _scatteringangle; // scattering angle
  double _dist;
  
  // variable constants used to dynamically determine simulation step-size
  double _K_max_var;
  double _lambda_var;
  double _beta_var;
  
  int _debug;
  int _status;
  int _interaction;
  
  double _volts_per_cm;
  double _time_to_collision;
  double _total_time;
  double _angle_param;
  double _density_param;
  
  double _energy;
  int _child_ions;
  
  std::mt19937& generator;
  
  int index(double v) const;
  int index_diff(double v) const;
  
  /*
   * Given a velocity, returns the effective scattering cross-section at the 
   * kinetic energy associated with the input velocity.
   * 
   * @param v The input velocity (in m/s)
   * @return The effective scattering cross-section (in cm^2)
   */
  inline double x_sec_e(double v) const { return (is_gas) ? momentum_xsec_gas[index(v)] : effective_xsec_liquid[index(v)]; }
  
  /*
   * Given a velocity, returns the elastic cross-section at the 
   * kinetic energy associated with the input velocity.
   * 
   * @param v The input velocity (in m/s)
   * @return The momentum-transfer cross-section (in cm^2) and scattering angle [deg]
   */
  double x_sec_p_diff(double v) const
  {
    
    auto idx = index_diff(v);
    double xsec = 0.;
    if (is_gas == false)
      xsec = momentum_xsec_liq_diff[idx];
    else
      xsec = momentum_xsec_gas_diff[idx];
    
    //if (status) { std::cout << std::scientific << "IDX = " << idx << " and XSEC = " << xsec << std::endl; }
    
    if (_debug) {
      std::cout << std::scientific;
      if (is_gas == false)
	std::cout << "[ DEBUG ] xsec @ idx " << 1 << " is " << momentum_xsec_liq_diff[0] << std::endl;
      else
	std::cout << "[ DEBUG ] xsec @ idx " << 1 << " is " << momentum_xsec_gas_diff[0] << std::endl;
      std::cout << "[ DEBUG ] xsec @ idx " << idx << " is " << xsec << std::endl;
    }
    
    return xsec;
  }
  
  /*                                                                                                                     
   * Given a velocity, returns the elastic cross-section at the                                                                  
   * kinetic energy associated with the input velocity.                                                                            
   *                                                                                                                                
   * @param v The input velocity (in m/s)                                                                                 
   * @return The momentum-transfer cross-section (in cm^2) and scattering angle [deg]                                 
   */
  double x_sec_e_diff(double v) const
  {
    
    auto idx = index_diff(v);
    
    double xsec = 0;

    if (is_gas == false)
      xsec = energy_xsec_liq_diff[idx];
    else
      xsec = energy_xsec_gas_diff[idx];

    //if (status) { std::cout << std::scientific << "IDX = " << idx << " and XSEC = " << xsec << std::endl; }
    
    if (_debug) {
      std::cout << std::scientific;
      if (is_gas == false)
	std::cout << "[ DEBUG ] xsec @ idx " << 1 << " is " << energy_xsec_liq_diff[0] << std::endl;
      else
	std::cout << "[ DEBUG ] xsec @ idx " << 1 << " is " << energy_xsec_gas_diff[0] << std::endl;
      
      std::cout << "[ DEBUG ] xsec @ idx " << idx << " is " << xsec << std::endl;
    }
    return xsec;
  }

    /*
     * Given a velocity, returns the momentum-transfer
     * cross-section and randomly sampled scattering angle
     * 
     * @param v The input velocity (in m/s)
     * @return The momentum-transfer cross-section (in cm^2) and scattering angle [deg]
     */
    std::pair<double,double> x_sec_p_angle_diff(double v) const
    {

      //std::cout << std::scientific;
      
      auto idx = index_diff(v);

      std::uniform_real_distribution<double> uniform(0.0, 1.0);

      double xsec = momentum_xsec_gas_diff[idx];

      // relative differential max of xsec to make sampling more efficient
      double relmax = 0.; 
      if (is_gas == false) { relmax = *(std::max_element(std::begin(momentum_xsec_liq_diff) + idx, std::begin(momentum_xsec_liq_diff) + idx + 180)); }
      else { relmax = *(std::max_element(std::begin(momentum_xsec_gas_diff) + idx, std::begin(momentum_xsec_gas_diff) + idx + 180)); }

      if (_debug) {
	std::cout << "[ DEBUG ] maxium relative cross-section for cross-section at index " << idx << " is " << relmax << std::endl;
	std::cout << "[ DEBUG ] xsec max for velocity " << v << " and bin " << idx << "is" << relmax << std::endl;
      }

      // now find the direction of the scatter
      bool direction = false;
      double angle = 0;
      while (direction == false) {
	angle = uniform(generator) * 180.;
	double prob  = uniform(generator) * relmax;
	double relxsec = 0.;
	if (is_gas == false) { relxsec = momentum_xsec_liq_diff[idx + int(angle)]; }
	else { relxsec = momentum_xsec_gas_diff[idx + int(angle)]; }
	if (_debug) std::cout << "[ DEBUG ] for angle " << angle << " we simulated a probability of " << prob << " and the local rel xsec prob is " << relxsec << std::endl;
	if (prob < relxsec)
	  direction = true;
      }
      return std::make_pair(xsec,angle);
    }

    /*                                                                                            
     * Given a velocity, returns the energy-transfer
     * cross-section and randomly sampled scattering angle
     *                                                                                                           
     * @param v The input velocity (in m/s)                                                                     
     * @return The momentum-transfer cross-section (in cm^2) and scattering angle [deg]                                
     */
    std::pair<double,double> x_sec_e_angle_diff(double v) const
    {

      auto idx = index_diff(v);
      std::uniform_real_distribution<double> uniform(0.0, 1.0);      
      double xsec = energy_xsec_gas_diff[idx];
      // relative differential max of xsec to make sampling more efficient

      double relmax = 0.;
      if (is_gas == false) { relmax = *(std::max_element(std::begin(energy_xsec_liq_diff) + idx, std::begin(energy_xsec_liq_diff) + idx + 180)); }
      else { relmax = *(std::max_element(std::begin(energy_xsec_gas_diff) + idx, std::begin(energy_xsec_gas_diff) + idx + 180)); }
      
      if (_debug) {
        std::cout << "[ DEBUG ] maxium relative cross-section for cross-section at index " << idx << " is " << relmax << std::endl;
        std::cout << "[ DEBUG ] xsec max for velocity " << v << " and bin " << idx << "is" << relmax << std::endl;
      }
      
      // now find the direction of the scatter                                                                                                               
      bool direction = false;
      double angle = 0;
      while (direction == false) {
        angle = uniform(generator) * 180.;
        double prob  = uniform(generator) * relmax;
        double relxsec = 0.;
	if (is_gas == false) { relxsec = energy_xsec_liq_diff[idx + int(angle)]; }
	else { relxsec = energy_xsec_gas_diff[idx + int(angle)]; }
        if (_debug) std::cout << "[ DEBUG ] for angle " << angle << " we simulated a probability of " << prob << " and the local rel xsec prob is " << relxsec << std::endl;
        if (prob < relxsec)
          direction = true;
      }
      return std::make_pair(xsec,angle);
    }
    
    /*
     * Given a velocity, returns the momentum-transfer cross-section at the 
     * kinetic energy associated with the input velocity.
     * 
     * @param v The input velocity (in m/s)
     * @return The momentum-transfer cross-section (in cm^2)
     */
    inline double x_sec_p(double v) const { return (is_gas) ? momentum_xsec_gas[index(v)] : momentum_xsec_liquid[index(v)]; }

    /*
     * Given a velocity, returns the ionization cross-section at the 
     * kinetic energy associated with the input velocity.
     * 
     * @param v The input velocity (in m/s)
     * @return The ionization cross-section (in cm^2)
     */
    inline double x_sec_i(double v) const { return ionization_xsec[index(v)]; }
     

    /*
     * Given a velocity, returns the excitation cross-section at the 
     * kinetic energy associated with the input velocity.
     * 
     * @param v The input velocity (in m/s)
     * @param level The excitation energy 'band'
     * @return The excitation cross-section (in cm^2)
     */
    inline double x_sec_ex(double v, int level) const
    {

      auto idx = index(v);

      if (_debug) {
	std::cout << std::scientific
		  << "[ DEBUG ] velocity " << v << " -> index " << idx << std::endl;
      }
      
        switch (level) {
	  
            case 11:
                return excite_xsec_11[idx];
            case 13:
                return excite_xsec_13[idx];
            case 14:
                return excite_xsec_14[idx];
	    case 15:
	      return excite_xsec_15[idx];
            default:
                assert(false);
        }
    }

  double probability(double u, double x_sec);
  double next_collision(double u);
  
  void update_pos_vel();
  void remove_energy(double eV);
  void ionization(std::vector<Electron*>& electron_list, int& total_ionizations);
  void elastic_collision(double u, Vec vm);
  void angle_scatter(double angle);
  void momentum_collision(double u, Vec vm);

public:

  Electron(double initial_time, double volts, Vec position, Vec velocity, std::mt19937& gen, int debug, int status, double angle, double density);
  
  inline Vec position() const { return _x; }
  inline Vec velocity() const { return _v; }
  inline double elapsed_time() const { return _total_time; }
  inline double timestep() const { return _time_to_collision; }
  inline double ke() const { return _energy; }
  inline double angle() const { return _scatteringangle; }
  inline double distancesincelastinteraction() const { return _dist; }
  
  void update(std::vector<Electron*>& electron_list, int& total_ionizations);
};

void generate_plot(int volts, double elec_energy, double angle, double density, double cutoff, int cores, int write_every, int k, int batches, int debug, int status, ProgressBar& bar);

#endif
