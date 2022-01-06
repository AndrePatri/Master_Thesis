% USED TO CALL gazebo plotter function and save results
close all
clear all
[figures_handles,graphs_handles]=graph_plot_Gazebo();
savefig(figures_handles,'/home/andrea/Desktop/università/magistrale/tesi/simulations/thesis_sim/Task_space_tracker/simple_relative_CoM_tracking/gazebo/sim_n4/matlab_sim.fig'); %saves graphs figures to file
savefig(figures_handles,'/home/andrea/Desktop/università/magistrale/tesi/simulations/thesis_sim/ISMPC_tracking/gazebo/general_traj_track/matlab_sim.fig'); %saves graphs figures to file
