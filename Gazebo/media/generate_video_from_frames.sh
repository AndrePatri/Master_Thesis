#!/bin/bash

ffmpeg -r 30 -pattern_type glob -i 'video_frames/default_camera_camera_link_my_camera*.jpg' -c:v libx264 generated_videos/sim_video.mp4