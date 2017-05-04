Creates MeerKAT subarray using simulations. Outputs array layout, uv-coverage and PSF

Requires
-----
[casapy](http://casa.nrao.edu/casa_obtaining.shtml)   
[simms](https://github.com/SpheMakh/simms)   
makems and pyxis from [radio-astro](https://launchpad.net/~radio-astro/+archive/ubuntu/main)

Usage
-----
select_subarray.py -- generating sets of antennas for subarrays (rerun until satisfied)   
build_array.py -- script to show subarray and MeerKAT array layout   
create_simms_inputfile.py -- autogen .cfg file used by simms   
make_cfg.py -- generate measurement set to get ANTENNA table for configuration

