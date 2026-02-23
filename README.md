# Kettera - A numerical solver for the time-dependent Schrödinger equation

Kettera (from the finnish word ketterä, meaning agile or nimble and from 'ket', the name commonly used for quantum state vectors) is 
a numerical solver for the time-dependent Schrödinger equation (TDSE) written in modern Fortran. It also comes equipped with
python scripts for generating plots and animations based on the data produced by the solver. Kettera is intented to be used on
GNU/Linux systems and no guarantees for functionality are provided for users on other operating systems.

## Features

Kettera can handle multiple different initial conditions and potentials to simulate many different quantum systems.
Initial conditons can also be read in from an external file. The user can also choose which numerical method is used
to solve the TDSE. Kettera implements the Crank - Nicolson method as well as the split-step Fourier method. The python
scripts that come equipped with Kettera can be used to plot certain physical quantities as a function of time and for
creating animations that showcase the dynamics 

### Initial conditions

The initial state of the system can be provided in two ways:
    1. As normal input parameters describing the type, initial offset and initial momentum of the wave function
    2. As an external file that contains the shape of the wave function at time t = 0

### Potentials

The user can decide what kind of environment the wavefunction is in by providing the parameters for a potential function.
The user can control the type, offset and strength of the potential function or choose to set the potential to zero
to simulate a free particle

### Visualization

After a simulation is done, the user can plot many physical quantities describing the system (eg. total energy, wavefunction normalization, etc.)
as a function of time using the python script 'ketplot.py'. The user can also generate high quality animations of the simulation using the
python script 'ketanim.py'. They can choose which quantities to animate and a gif of the dynamics will be produced

## External Dependencies

Kettera uses the [FFTW](https://www.fftw.org/) library to calculate Fourier transforms and the python scripts require
[Numpy](https://numpy.org/) as well as [Matplotlib](https://matplotlib.org/) to function.

## Installation

NOTE: This section assumes you have configured the external dependencies listed in the previous section to work on your system

## Usage

To invoke Kettera, simply type `kettera` into a terminal. This should produce some simple usage information about kettera.
Providing the option `--help` will provide more detailed usafe information about Kettera and providing the option `--version`
will output Kettera's current version

### Running simulations

### Using external input files

### Using the visualization tools

## License

Kettera is provided under [version 3 of the GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Author

Kettera is written by Juuso Kaarela in 2026


