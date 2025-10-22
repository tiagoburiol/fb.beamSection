## Mesh generation with GiD and MAT-Fem

---

### GiD

GiD is a user-friendly pre and post-processing software developed by the International Center for Numerical Methods in Engineering [(CIMNE)](www.cimne.com) that is easy to use and customize.  The software handle meshes with 1D, 2D and 3D elements and the free version can export meshes with up to 10,000 nodes. More information, manuals and download is available at [GiD Simulation](https://www.gidsimulation.com/gid-for-science/). 

#### Download and Installation

If you're having issues installing GiD, please refer to [the official GID installation page](https://www.gidsimulation.com/gid-for-science/support/installation/).

### MAT-Fem

Along with GiD, we use the [MAT-Fem](https://www.gidsimulation.com/downloads/educational-finite-element-codes-matfem/) plugin to assign material properties and export the mesh as a Python-ready file. The export is done with a `MAT-Fem.bas` edited file that can be found [in this folder](https://github.com/tiagoburiol/napy.torsion/tree/main/mesh_files). MAT-Fem is a GiD plugin developed by researchers at [CIMNE](https://cimne.com/) and made to be a Finite Element utility to generate meshes to be used with MATLAB codes.

#### Download and Installation

To download MAT-Fem visit the [GiD downloads page](https://www.gidsimulation.com/downloads/educational-finite-element-codes-matfem/).

* You can install MAT-Fem by simply locating the `problem types` folder on GiD's installation folder (default: C:\Program Files\GiD\GiD <version>\problemtypes)

* As a plugin MAT-Fem, can be also be used from __*Data > Problem Type > Load*__  tab on GiD. On the *Read Problem Type* window navigate to the folder where `MAT-Fem.gid` is downloaded and click **Open**. 

#### .bas File

* You can also place the `MAT-fem_Python_Multimat_dens.bas` file from this repository in the Tamplates folder on GiD's installation folder (default: C:\Program Files\GiD\GiD <version>\templates) if you wish for it to be more easily accessed inside GiD as one of the options at ***Files > Export > Using template .bas (only mesh)***. 

* Other-wise you'll have to go to t ***Files > Export > Using template .bas (only mesh) >  Others...*** and navigate to where the .bas file is located at every time you wish to export a mesh. 

### Example

 After opening GiD, If you have MAT-Fem already installed in the`problem types` folder go to ***Data > Problem Type > MAT-Fem***, or go to __*Data > Problem Type > Load*__  and find where the MAT-Fem.gid file is. A welcome window (Fig. 1) wil appear and a new sidebar (Fig. 2) will be available.

| ![]() |     |
| ----- | --- |


