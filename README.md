# UCNP-RET-MC
A small utility for simulating resonance energy transfer from upconverting nanoparticles (UCNPs) to organic dyes, using a Monte Carlo approach. If you are using this code, please cite Dukhno et al. (2017), DOI:10.1039/C6NR09706E, and the original article containing the algorithm, Corry et al. (2005) DOI:10.1529/biophysj.105.069351.

### Concept
Resonance energy transfer in a multi-donor-multi-acceptor system can be complicated to express in a closed-form solution, especially in relatively saturated conditions where excitation and emission might be hindered due to the donors or acceptors spending a significant time in an excited state. Additionally, a closed-form solution wouldn't allow to estimate differences between systems with slightly different geometries. Meanwhile, a naive Monte Carlo model with a fixed time step could work when the lifetimes of donors and acceptors are of the same order of magnitude, but would be expensive to scale to long timeframes and high time granularities.

The algorithm described in Corry et al. (2005) allows to calculate only the "interesting" events while keeping track of possible loss of excitation and energy transfer due to donors/acceptors being in excited states.
This is an implementation of this algorithm to the particular case of Yb-Er doped UCNPs. More details are provided in Dukhno et al. (2017).

### Usage
Run UCNP_RET_montecarlo.py from command line or UCNP_RET_montecarlo.ipynb from JupyterLab. Follow the prompts and fill them in with your corresponding parameters. 

The script allows working with several particle geometries: uniformly doped particle, core-doped particle and shell-doped particle. 

The script saves your configuration to an XML file (keep a copy in case you want to repeat a calculation later). Default values for  a run are provided in respective XML files for different particle geometries.

In case you want to visualise particle geometry, use the `plot_geometry` function.

### Notes, precautions and caveats 
- Excitation flux: currently the script does not calculate excitation flux for the particle. An example of this calculation is provided in the original paper.
- The calculation of donor and acceptor quantities is hard-coded for beta-NaYF4 only. Crystal structure parameters for different kinds of particles need to be appropriately parametrized.

### Dependencies
Numpy (works with 1.17.0, haven't tested other versions).
Matplotlib (works with 3.1.1, haven't tested other versions).
