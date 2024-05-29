# TopologyOptimization-FORTRAN

This repository contains a series of FORTRAN routines to perform a topology optimization using FEM in linear elasticity. The code operates by inputting structure and load information into the "DataStructure" folder and defining characteristics in the main.f90 file. It utilizes the MA87 library from HSL (with an adapted interface), as well as the LAPACK, BLAS, Metis, and OpenMP libraries. The code uses either the method of moving asymptotes (MMA) or the optimality criterion (OC) presented by Professors Svanberg (1987) and Bendsoe & Sigmund (2013). This code is for academic purposes only.

A .makefile is included to facilitate compilation and execution. During the post-processing stage, the code generates files in EnSight format (https://dav.lbl.gov/archive/NERSC/Software/ensight/), which can be read using ParaView, a free and open-source tool (https://www.paraview.org/). The code was based on the work of Liu & Tovar (2014) and Andreassen et al. (2011). Any suggestions for improving the code performance or reports of bugs are welcome.

- Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B. S., & Sigmund, O. (2011). Efficient topology optimization in MATLAB using 88 lines of code. Structural and Multidisciplinary Optimization, 43, 1-16.
- Bendsoe, M. P., & Sigmund, O. (2013). Topology optimization: theory, methods, and applications. Springer Science & Business Media.
- Svanberg, K. (1987). The method of moving asymptotes—a new method for structural optimization. International journal for numerical methods in engineering, 24(2), 359-373.
- Liu, K., & Tovar, A. (2014). An efficient 3D topology optimization code written in Matlab. Structural and multidisciplinary optimization, 50, 1175-1196.
