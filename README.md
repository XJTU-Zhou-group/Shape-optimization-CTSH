# Shape-optimization-CTSH
Created by XJTU-Zhou-group, May, 2022


The ABAQUS/Explicit script folder contains calculation scripts for the CTSH folding and unfolding quasi-static process, which can automatically output the maximum unfolding moment peak, mass, maximum storage strain energy, maximum curvature, and maximum failure indices.

The four folders beginning with Optimization contain the codes of constructing the four agent models(ANN, RBF, GPR and Kriging) and the isight integration files that call the surrogate model for optimization.

The sample point data folder contains the original data calculated from the ABAQUS script used for the training and testing of the surrogate models.

The DOE & Calculating data set by isight folder contains the integrated files of isight and the abaqus python script, which could complete DOE sampling calculations automatically.
