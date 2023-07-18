# DefectNetwork
This repository simulates the formation of defects network near a saturated grain boundary


The project consists of three parts: 
1. Simulating multiple batches of **prismatic loop evolution** using kMC algorithm embedded with Dislocation Dynamics module [dd3d_diffuse.m](dd3d_diffuse.m), all results from different batches are 'parallel' and stiched together using [merge.m](merge.m) and the [dislocation data files](dislocation.txt) are exported from it.
2. Generating **cubic volume tessellation** with the [tessellation_orthogonal.m](tessellation_orthogonal.m). The [output](orthogonal20/n20_vorvx0.txt) contains all the vertices of each volume element. 
3. Conducting **GND signal analysis** with [main.py](GND rendering/main.py), which reads both the [dislocation data](dislocation.txt) and [volume mesh data](n20_vorvx0.txt). It outputs information about dislocation segments after being truncated by volume cells, the volume, vertice and GND content of each volume cell.
4. Generating **GND signal 3D map** with [voro_plot.m](GND_3d_analysis/voro_plot.m) and [result](orthogonal_merged20_voro_color.txt) is rendered in [GND rendering/3d_rendering.py](GND rendering/3d_rendering.py) for final 3d image of GND signal near the GB.


## Table of Contents

- [Install](#install)
- [Usage](#usage)
- [Examples](#example)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)


## Install

This project consists of code developed in [Python](https://www.python.org/) and [MATLAB](https://www.mathworks.com/products/matlab.html). Commercial Software [OVITO](https://www.ovito.org/) is used for data processing. 


## Usage

### Dislocation extraction
Starting from the [sample MD result]() of Fe grain boundary with a prismatic loop, we use OVITO to export XYZ file that contains coordinates of all particles. The set up for OVITO and the dislocation extraction image are as shown below:

<img src="ovito_setup.png" width="300" height="700">     <img src="dislocation.png" width="300" height="700"> 

The output file should look like the [dislocation.XYZ](). Then the python script [dxa_analysis.py](dxa_analysis.py), which is built upon [OVITO's Python interface](https://docs.ovito.org/python/) extracts dislocation information from the XYZ file and outputs dislocation data for the types '100', '110', '111' and other into txt files. 

### Spatial tessellation
There are two types of meshes available: [Voronoi Tessellation](tessellation_voronoi.m) and [regular Hexahedron Tesselation](tessellation_cubic.m). To construct spatial tessellation, run the corresponding MATLAB code and the coordinates of all volume mesh vertices will be saved into text files in designated directories.

### Gaussian GND calculation
The [main.py](main.py) calculates the GND density by combining the dislocation information and meshed volumes. The output from this script includes: 
1. Simulation parameters like the mesh size and the dimension of the simulation volume; 
2. Volume of each meshed volume element; 
3. GND signal intensity of each meshed volume element; 
4. Vertices of dislocation segments with the truncating effect of meshed volume elements; 
5. The correspondence between each dislocation segment and the meshed volume element that fully contains it.

### Analyze GND calculation

## Example



## Maintainers

[@AlanHe](https://github.com/hsc1993).

## Contributing


### Contributors

This project is supported by research group at UCLA, Johns Hopkins University, Hongkong City University and Pennsylvania University.


## License

[MIT](LICENSE) © Sicong He





