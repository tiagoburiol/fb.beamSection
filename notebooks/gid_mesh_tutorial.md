# Mesh generation with GiD and MAT-Fem

---

### GiD

GiD is a user-friendly pre and post-processing software developed by the International Center for Numerical Methods in Engineering [(CIMNE)](www.cimne.com) that is easy to use and customize.  The software handle meshes with 1D, 2D and 3D elements and the free version can export meshes with up to 10,000 nodes. More information, manuals and download is available at [GiD Simulation](https://www.gidsimulation.com/gid-for-science/). 

#### Download and Installation

If you're having issues installing GiD, please refer to [the official GID installation page](https://www.gidsimulation.com/gid-for-science/support/installation/).

### MAT-Fem

Along with GiD, we use the [MAT-Fem](https://www.gidsimulation.com/downloads/educational-finite-element-codes-matfem/) plugin to assign material properties and export the mesh as a Python-ready file. The export is done with a `MAT-Fem.bas` edited file that can be found [in this folder](https://github.com/tiagoburiol/fb.beamSection/tree/main/mesh_files). MAT-Fem is a GiD plugin developed by researchers at [CIMNE](https://cimne.com/) and made to be a Finite Element utility to generate meshes to be used with MATLAB codes.

#### Download and Installation

To download MAT-Fem visit the [GiD downloads page](https://www.gidsimulation.com/downloads/educational-finite-element-codes-matfem/).

* You can install MAT-Fem by putting the downloaded file on the `problem types` folder on GiD's installation folder (default: *C:\Program Files\GiD\GiD <version>\problemtypes*)

* MAT-Fem, can also be initialized from __*Data > Problem Type > Load*__  tab on GiD. On the ***Read Problem Type*** window navigate to the folder where `MAT-Fem.gid` is downloaded and click **Open**. 

#### The .bas File

The .bas file is a template for exporting the mesh in a format suitable for use in FlightBEND.

* You can also place the `MAT-fem_Python_Multimat_dens.bas` file from this repository in the Tamplates folder on GiD's installation folder (default: *C:\Program Files\GiD\GiD <version>\templates*) if you wish for it to be more easily accessed inside GiD as one of the options at ***Files > Export > Using template .bas (only mesh)***. 

* Other-wise you'll have to go to t ***Files > Export > Using template .bas (only mesh) >  Others...*** and navigate to where the .bas file is located at every time you wish to export a mesh. 

## Steps

---

 After opening GiD, If you have MAT-Fem already installed in the`problem types` folder go to ***Data > Problem Type > MAT-Fem***, or go to __*Data > Problem Type > Load*__  and find where the MAT-Fem.gid file is. A welcome window (**Fig. 1**) wil appear and a new sidebar (**Fig. 2**) will be available. 

| ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/MAT-Fem_welcome.png "MAT-Fem Welcome") | ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/MAT-Fem_sidebar.png "MAT-Fem Sidebar") |
| ----------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| **Fig. 1** -  MAT-Fem welcome window.                                                                             | **Fig. 2** -  MAT-Fem sidebar.                                                                                    |

### Draw Geometry

Use the ***create line*** ![](https://github.com/tiagoburiol/fb.beamSection/tree/main/notebooks/imgs/icons/create_line.png =16x16 "create line" )tool to draw straight lines; the ***create arc***  ![](https://github.com/tiagoburiol/fb.beamSection/tree/main/notebooks/imgs/icons/create_arc.png =16x16 "create arc")for circular arcs and ***create nurbs line*** ![](https://github.com/tiagoburiol/fb.beamSection/tree/main/notebooks/imgs/icons/create_nurbs_line.png =16x16 "create nurbs line") for curves (splines). More options are also available on the ***Geometry*** tab (**Fig. 3**). GiD also has capabilities to import geometry from CAD editors and various formats on the ***File > Import*** tab (**Fig. 4**). 

> Draw the geometry only on the XY plane on GiD. But the mesh will be exported in the YZ plane.

| ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/GiD_geometry_tab.png "Geometry") | ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/GiD_import.png "Import") |
| ------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------- |
| **Fig. 3** -  Geometry tab.                                                                                   | **Fig. 4** -  Import options.                                                                         |

### Define Surfaces

It is important that you define *surfaces*  in the geometry so GiD undertands that 2D element meshes should be generated. To define a surface click the ***Create NURBS surface*** ![](https://github.com/tiagoburiol/fb.beamSection/tree/main/notebooks/imgs/icons/create_nurbs_surface.png =16x16 "Create nurbs surface") tool, select only the parts of the geometry that is closed by clicking and holding the left mouse button (**Fig. 5**) (make sure that the geometry is *closed*, otherwise the surface can't be defined) then press the ***Esc*** keyboard key to confirm your selection. *<mark>Remember that in GiD, to confirm a action you use the **Esc** key! </mark>* After that, if a valid surface is crated, a smaller copy of the geometry will appear in magenta color inside the original drawing (**Fig. 6**). 

| ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/GiD_surface_select.png "Selecting geometry with NURBS surface tool") | ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/GiD_surface.png "Surface created") |
| ----------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- |
| **Fig. 5** -  Selecting geometry to define a surface.                                                                                           | **Fig. 6** -  Surface created successfully.                                                                   |



### Assign Material

The most convenient way to assign material to geometry is before mesh creation, by assigning the material properties to surfaces. This way, if you wish/need to regenerate the mesh the assigment stays defined. To do that, with the MAT-Fem sidebar, open click the ***Assign Material Properties***  ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/icons/assign_material.png "Assign Material Properties" =16x16)button. The Materials window will appear (**Fig. 7**). You can directly edit the properties currently shown on the window and save you own material definition if you wish. GiD does export measument units, only the values  as they appear in the material window, so <mark>make sure to remember the units of length, stress and mass you used</mark>. 

> Note: the thickness property is not used in our program, so it can be left as is. 

After filling the properties, click the ***Assign*** button then the ***Surfaces*** option (**Fig. 7**). Now you can select multiple surfaces to assign the current properties to. After selection (**Fig. 8**) press the ***Esc*** key to confirm. On the materials window you can alse use the ***Draw*** button to show the assignment of specific materials or all materials (**Figs. 9** and **10**) to verify if the assignment step was done correctly.

| ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/MAT-Fem_material_window.png "Materials window")             | ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/MAT-Fem_material_window_surface_select.png "Surface selected for material assignment") |
| -------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Fig. 7** -  MAT-Fem materials window.                                                                                                | **Fig. 8** -  Surface selected for material assignment.                                                                                                            |
| ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/MAT-Fem_material_window_draw_material.png "Draw materials") | ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/MAT-Fem_material_window_draw_material_show.png "Displaying assigned materials")         |
| **Fig. 9** -  Draw material option.                                                                                                    | **Fig. 10** -  Displaying assigned materials.                                                                                                                      |

### Assign Element Type and Order

 By default, GiD generates 2D meshes with linear triangular elements, if you want other element type go to ***Mesh > Element Type*** and select either Triangle, or Quadrilateral (**Fig. 11**). **Important:** *To correctly generate quad meshes, make sure the surfaces are bounded by only 4 lines (which may be curved or straight), otherwise the mesh can't be generated*. After that you'll be prompted to select the geometry. After selection press the ***Esc*** key to confirm. To change the intepolation order go to **Mesh > Quadratic Type** and select between ***Normal (Linear)***, ***Quadratic*** or ***Quadratic9*** (**Fig. 12**). The last option is for complete quadratic Lagrange elements Tri-6 and Quad-9, the second is for Quad-8 serendipidy elements. The interpolation option is automaticaly applied to the entire mesh by GiD to ensure compatibility.

| ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/GiD_element_type.png "Element type selection") | ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/GiD_element_order.png "Element order selection") |
| ------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| **Fig. 11** -  Selecting element type.                                                                                    | **Fig. 12** -  Selecting interpolation order.                                                                               |

### Assign Mesh Parameters

GiD is able to generate both unstructured (**Fig. 13**) and structured (**Fig. 14**) meshes and you have control over a few parameters so that the resulting mesh is adequate for your application. Under the options for ***structured*** meshes you'll have to select either ***lines*** or ***surfaces***  to assign parameters to. For surfaces you may choose to assign ***number of divisions*** or ***element size*** to any ***sides*** of a surface (**Fig. 14**). The same choices are applicable to the lines option. 



When assigning structured parameters to a surface take notice that the following steps are required: 

1. First a dialogue box will appear for you to type the value of the parameter;

2. After that you'll be prompted to *choose the surfaces*, and only _after pressing the Esc key_ you may *select the sides to apply the paramenter*. 

3. After selecting the sides pressing Esc will apply the parameters and open a new dialogue, returning to Step 1.

4. To exit you can close the dialogue window.
   
   

If you wish to generate a unstructured mesh on a given surface, you may just select ***Mesh > Unstructured > Surfaces*** and select which surfaces may be left unstructured.

| ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/GiD_unstructured_mesh.png "Unstructured mesh") | ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/GiD_structured_mesh_surface.png "Structured surface mesh options") |
| ---------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Fig. 13** -  Unstructured mesh options.                                                                                    | **Fig. 14** -  Structured surface mesh options.                                                                                                  |

### Generate Mesh

To generate the mesh you may use the ***Generate Mesh*** ![](https://github.com/tiagoburiol/fb.beamSection/tree/main/notebooks/imgs/icons/generate_mesh.png "Generate Mesh") button from MAT-Fem, or access the tab ***Mesh > Generate Mesh***, or the keyboard shortcut Ctrl + G. The mesh generation window will appear (**Fig. 15**). The element size option will only apply to geometry that has not been given structured parameters. Type a element size of leave it as-is then click on **Ok**. The ***Progress in meshing*** window will appear and display a graph with the evolution of the meshing process over time. Pay special attention to the contex box above the graph, there will be show the type and number of elements generated. *<mark>Make sure no **line elements** are generated</mark>*. I there are, get rid of any single lines that are not conected to surfaces and mesh again. Also, <mark>*make sure only **a sigle type of element** is present*</mark> -- either triangles or quadrilaterals -- because our routine assumes the mesh is composed only of a single element type. Click on **View mesh** to see the result.

| ![](https://github.com/tiagoburiol/fb.beamSection/blob/main/notebooks/imgs/GiD_mesh_generation_window.png "Generate mesh window") | ![](ttps://github.com/tiagoburiol/napy.torsion/blob/main/notebooks/imgs/GiD_progress_in_meshing.png "Progress in meshing ") |
| --------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| **Fig. 15** - "Generate mesh" window.                                                                                                   | **Fig. 16** -  "Progress in meshing" window.                                                                                |

### Export Mesh

After generation the mesh can be exported as various file types on **Files > Export**. To export a mesh compatible with FlightBEND use the `.bas` file provided in this folder and save it on the `templates` folder on GiD's intallation folder. This way, to export the mesh go to ***File > Export > Using template .bas (only mesh) > MAT-fem_Python_Multimat_dens.bas***. A dialog box will appear; navigate to a desired folder and *<mark>**don't forget to name the file with a ``.py``  extension**</mark>*.

# Structured Mesh Exemple

---

Figs. 17 to XX show the steps to generate a structured mesh with 4 elements across the horizontal lines and 8 elements across the vertical lines of a trapezoidal geometry.

1. Use the line tool to draw a trapezoid like shown in **Fig. 17**. When clicking the last point on the last edge the **Create poin procedure** window may appear (**Fig. 18**); click on **Join** to make sure the geometry will be closed. Press **Esc** to exit the line tool.

2. Click the **Create NURBS surface** button, select the geometry then press **Esc** to confirm. 

3. Assign material properties on the surface with the **Assign Material Properties** button on the MAT-Fem side bar; Enter the desired values; click on the button **Assign > Surfaces**; select the surface and press **Esc** to confirm. 
   
   > If the sidebar is not present, go to ***Data > Problem type > MAT-Fem***, if you have MAT-Fem in the *Problem types* folder on the GiD installation folder; or ***Data > Problem type > Load...*** and navigate to the where *MAT-Fem.gid* is downloaded. 

4. Select element type on **Mesh > Element Type > Quadrilateral** then select the surface and press **Esc** to confirm.

5. Select element order on **Mesh > Quadratic type > Quadratic9**  (to assign Quad-9 elements)
   
   > You may choose to assign element type **before** or **after** material assignment, not necessarily in the order shown in this tutorial.

6. Go to **Mesh > Structured > Surfaces > Assign number of divisions to surface lines**

7. Select the trapezoidal surface.

8. Press **Esc** on the keyboard.

9. Type 4 in the box and click **Assign**. 

10. Select the upper edge. 
    
    > Notice the bottom edge will also be selected, since in a structured mesh opposite sides must have the same number of elements.

11. Press **Esc** on the keyboard. This will assign 4 divisions on the selected lines. After, the Dialog window will reappear.

12. To assign a new size to the vertical lines, type 8 in the box then click **Assign**.

13. Select either the left of the right edge. 
    
    > Notice the other edge will also be selected.

14. Press **Esc** on the keyboard. This will assign 8 divisions on the selected lines. After, the Dialog window will reappear.

15. Close the dialog box.

16. Click the **Generate mesh buttom** then **Ok**. Make sure there are only quadrilateral elements then click **View Mesh**.

17. Export the mesh on ***File > Export > Using template .bas (only mesh) > MAT-fem_Python_Multimat_dens.bas***. making sure to give it a name and a .py extension.

# Useful Tips for Geometry Editing

---

* GiD is a buggy program. Be sure to save your project often as you edit.

* The free version of GiD allows only a limited number of surfaces to be saved on a file, but you can still export the mesh and later delete a few surfaces to be able to save the project.

* At any point you can delete geometry by clicking the **Delete** button an then delecting the type of geometry to delete.

* After generanting the mesh, you may swich between mesh and geometry views using the **toggle geometry/mesh view** ![](https://github.com/tiagoburiol/fb.beamSection/tree/main/notebooks/imgs/icons/toggle_geometry_mesh.png "Toggle geometry/mesh") button.

* You can show node and element labels by clicking the **set all labels on/off** ![](https://github.com/tiagoburiol/fb.beamSection/tree/main/notebooks/imgs/icons/turn_labels.png "Labels on/off") button.

* You can subdivide lines with the **Divide line into a number of segments** ![](https://github.com/tiagoburiol/fb.beamSection/tree/main/notebooks/imgs/icons/divide_line.png "Divide line")tool. Note that the original line will be substituted by the new segments. You can intersect lines with the **Intersect lines**  ![](https://github.com/tiagoburiol/fb.beamSection/tree/main/notebooks/imgs/icons/intersect_line.png "Intersect line") tool. Similar tools are available for surfaces as well.

* Before mesh creation it can be useful to collapse points and lines that are too close together with the **Delete duplicate entities with a tolerance** ![](https://github.com/tiagoburiol/fb.beamSection/tree/main/notebooks/imgs/icons/delete_duplicate "Delete duplicate") tool. 
  
  

|                        |                                               |
| ------------------------------- | ------------------------------------------------ |
| **Fig. 17** - Draw a trapezoid. | **Fig. 18** -  Join the last point to the first. |
