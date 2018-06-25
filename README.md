# lidar-optimization

Requirement:
Gurobi solver

Abstract:
LiDARs plays an important role in self-driving cars and its configuration such as the location placement for each LiDAR can influence object detection performance. This paper aims to investigate an optimal configuration that maximizes the utility of on-hand LiDARs. First, a perception model of LiDAR is built based on its physical attributes. Then a generalized optimization model is developed to find the optimal configuration, including the pitch angle, roll angle, and position of LiDARs. In order to fix the optimization issue with off-the-shelf solvers, we proposed a lattice-based approach by segmenting the LiDAR’s range of interest into finite subspaces, thus turning the optimal configuration into a nonlinear optimization problem. A cylinder-based method is also proposed to approximate the objective function, thereby making the nonlinear optimization problem solvable. A series of simulations are conducted to validate our proposed method. This proposed approach to optimal LiDAR configuration can provide a guideline to researchers to maximize the utility of LiDARs.

![image](https://github.com/zhao-lab/mou-lidar-optimization-itsc18/blob/master/case1_2lidar_2laser_822.jpeg.001.jpeg)

In folder "tool", folder "c" and folder "matlab" are the source code offered by Gurobi.
You can't directly git clone these two folders and need to install the source code from Gurobi official website because there's a license issue.

Put the "main.m" in folder "src" into folder "matlab" installed and run it to get a ".lp" file which is the model created by Gurobi.

Then run "mip2_c.c" in folder "c" installed to get the positions of lidars in the model you designed in "main.m".
