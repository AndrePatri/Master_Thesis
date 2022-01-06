#!/bin/bash

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Modify according to the plugins (by adding/deleting paths) 
echo 'export GAZEBO_PLUGIN_PATH=${GAZEBO_PLUGIN_PATH}:'$SCRIPT_DIR'/plugins/my_proto_knee_ISMPC/build:'$SCRIPT_DIR'/plugins/my_proto_knee_Task_Tracker/build'>> ~/.bashrc
#Model paths
echo 'export GAZEBO_MODEL_PATH=$GAZEBO_MODEL_PATH:/usr/share/gazebo-11/models:'$SCRIPT_DIR'/models'>> ~/.bashrc

# Sourcing setup
echo 'source /usr/share/gazebo/setup.sh'>> ~/.bashrc

$SHELL