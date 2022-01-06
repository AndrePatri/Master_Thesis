#!/bin/bash
killall -9 gazebo & killall -9 gzserver  & killall -9 gzclient #kills any Gazebo instance (Gazebo will not run while another instance of it is already executing)

export GAZEBO_PLUGIN_PATH=${GAZEBO_PLUGIN_PATH}:$HOME/Gazebo_standalone/plugins/build # exporting Gazebo plugin path
export GAZEBO_MODEL_PATH=$GAZEBO_MODEL_PATH:/usr/share/gazebo-11/models:/usr/share/gazebo-11/models:$HOME/IIT_CRIS_HHCM/Gazebo_standalone/models # exporting Gazebo plugin path

gazebo -u $HOME/Gazebo_standalone/worlds/my_world/my_world.world --verbose

$SHELL #keep the terminal open for debugging purposes