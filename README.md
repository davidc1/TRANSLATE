```
mkdir build
cd build
cmake ..
make
```

### Setup
To run the code type `./LArCADe` and follow the prompts. To change running configuration edit the `srcs/Constants.h` file and specifically update the value of `n`, the density of argon atoms, `is_diff` and `is_gas` to choose the simulation mode.

Configurable constants:

```
// Liquid number density [cm^-3]                                                                                                                       
const double n = 2.11e22;       // Number density of argon atoms

// Editable parameters                                                                                                                                 
const bool is_diff = false;     // Use differential cross-sections?                                                                    
const bool is_gas = false;      // Is this in gas? Or liquid? Make sure to update the number density 
```

Default gas mode (integrated xsec): `n = 1.15e19`, `is_diff = false`, and `is_gas = true`.

Default liquid mode (integrated xsec): `n = 2.11e22`, `is_diff = false`, and `is_gas = false`.

### Simulation Output
Information for each simulated electron is stored in text files. Each iteration of the simulation for which output is stored is saved as a new row, and information in each row is arranged as described below, separated by commas:

time (ns), time-step (ns), x coordinate (um), y coordainte (um), z coordainte (um), kinetic energy (eV), drift velocity (m/s), scattering angle (degrees), distance since last interaction (meters), interaction category, total ionizations initiated

### Notebooks
Validations of swarm parameters are carried out through the notebooks below:

```
lar-diffusion-simple.ipynb
gar-diffusion.ipynb
lar-drift-velocity.ipynb
gar-drift-velocity.ipynb
lar-amplification-townsend.ipynb
gar-amplification-townsend.ipynb
```
