function [Q_sol,time,VV_k,U_middle,actuation_left,actuation_right,Y_ref,dY_ref,dt_control,cost]=LoadDataFromGazeboSim_ISMPC()
Q_sol = importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/ISMPC/Q_sol.txt');
time= importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/ISMPC/Time.txt');
actuation_left=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/ISMPC/U_left.txt');
actuation_right=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/ISMPC/U_right.txt');
Y_ref=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/ISMPC/Y_ref.txt');
dY_ref=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/ISMPC/dY_ref.txt');
dt_control=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/ISMPC/dt_control.txt');
cost=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/ISMPC/cost.txt');
U_middle=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/ISMPC/U_middle.txt');
VV_k=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/ISMPC/VV_k.txt');
end
