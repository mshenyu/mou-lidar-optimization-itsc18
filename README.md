# lidar-optimization

Requirement:
Gurobi solver

Abstract:
LiDARs plays an important role in self-driving cars and its configuration such as the location placement for each LiDAR can influence object detection performance. This paper aims to investigate an optimal configuration that maximizes the utility of on-hand LiDARs. First, a perception model of LiDAR is built based on its physical attributes. Then a generalized optimization model is developed to find the optimal configuration, including the pitch angle, roll angle, and position of LiDARs. In order to fix the optimization issue with off-the-shelf solvers, we proposed a lattice-based approach by segmenting the LiDAR’s range of interest into finite subspaces, thus turning the optimal configuration into a nonlinear optimization problem. A cylinder-based method is also proposed to approximate the objective function, thereby making the nonlinear optimization problem solvable. A series of simulations are conducted to validate our proposed method. This proposed approach to optimal LiDAR configuration can provide a guideline to researchers to maximize the utility of LiDARs.
