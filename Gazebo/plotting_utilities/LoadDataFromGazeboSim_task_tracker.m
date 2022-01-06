function [Q_sol,time,X_k,U_middle,actuation_left,actuation_right,Chi,Chi_dot,Chi_ref,Chi_dot_ref,Chi_ddot_ref,dt_control,cost]=LoadDataFromGazeboSim_task_tracker()
Q_sol = importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/task_tracker/Q_sol.txt');
time= importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/task_tracker/Time.txt');
actuation_left=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/task_tracker/U_left.txt');
actuation_right=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/task_tracker/U_right.txt');
Chi=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/task_tracker/Chi.txt');
Chi_dot=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/task_tracker/Chi_dot.txt');
Chi_ref=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/task_tracker/Chi_ref.txt');
Chi_dot_ref=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/task_tracker/Chi_dot_ref.txt');
Chi_ddot_ref=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/task_tracker/Chi_ddot_ref.txt');
dt_control=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/task_tracker/dt_control.txt');
cost=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/task_tracker/cost.txt');
U_middle=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/task_tracker/U_middle.txt');
X_k=importdata('/home/andrea/Desktop/università/magistrale/tesi/ROS_Gazebo/Gazebo/sim_results/task_tracker/X_k.txt');
end
