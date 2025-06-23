#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "LUTs/CrossSections.h"

// Liquid number density [cm^-3]
const double n = 2.11e22;

// Gas number densities
//const double n = 8.37e17;                       // at 10 mbar (0.01 atm) (7.5 Torr)
//const double n = 1.15e19;                       // at 2 PSI (103 Torr)
// const double n = 4.02e19;                       // at 7 PSI (362 Torr) (4.06 cm^-3 according to Van der Waals equation)
// const double n = 8.61e19;                       // at 15 PSI (776 Torr) (8.81 cm^-3 according to Van der Waals equation)

// Global constants
const double k_b = 1.38e-23;                    // Boltzmann constant (in J/K)
const double T = 87;                            // The argon temperature (in K)
const double M_a = 6.6e-26;                     // Mass of the argon atom (in kg)
const double m_e = 9.1e-31;                     // Mass of the electron (in kg)
const double e = 1.6e-19;                       // Charge of the electron (in C)
const double lambda = n * K_max;                // Inverse mean of the exponential distribution for drawing timesteps
const double beta = 1 / lambda;                 // Mean of the exponential distribution for drawing timesteps

// Editable parameters
const bool is_diff = false;                      // Use differential cross-sections?
const bool is_gas = false;                       // Is this in gas? Or liquid? Make sure to update the number density
const bool track_child_ions = false;            // Keep track of child electrons and their ionizations? Or only ionizations of the main, simulated electron?
const bool uniform_field = true;                // Is this in a uniform field? Or are we simulating a tip / tip array?
const bool interactions = true;                 // Simulate interactions?
const bool single_tip = true;                   // Is this is a non-uniform field, simulate a single tip? Or simulate an array of tips?
const double spawn_height_scale = 1.0;          // Within the simulated volume (when using a non-uniform field), where should the spawn disc be located? 1.0 is at the top, 0.0 is at the bottom.
const double spawn_width_scale = 1.0;           // Within the simulated volume (when using a non-uniform field), how should the radius of the spawn disc be scaled? 1.0 is the width of the simulated volume, 0.0 is a point.

#endif
