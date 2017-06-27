# VorticellaSimulations
A video of the full simulation is available here: https://www.youtube.com/watch?v=lZbrAN-_Qck

Compiling the code:

$ make all

You can change the compiler in the Makefile. Please check the Makefile before compiling the code.

Output of the code:
Once you run
$make all
the Makefile will create an executable named rod.
The output of rod will be the x,y,z coordinates of each element of the Vorticella stalk at each timepoint of the simulation. These coordinates are best viewed in VMD, which can be downloaded from here:
http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD

The output x,y,z coordinates are saved in gromax format in these files:

animation.gro
spasmo.gro
head.gro

All these files need to be opened in VMD simultaneously to view the full simulation.

The unit.xlsx (and unit.ods) file contains the units of all physical quantities employed in the model.
The main.c file employs slightly different units that the rod.c file. However, all unit conversions are taken care of within the codes.
the unit of distance is 1 micrometer in main.c, while it is 1 nanometer in rod.c.
the unit of time is 1 millisecond in main.c, while it is 1 microsecond in rod.c.
the unit of force is 1 femtonewton in main.c, it is 1 femtonewton in rod.c as well.
