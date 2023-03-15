#!/usr/bin/env python3

from nuscenes.nuscenes import NuScenes

# Load dataset
nusc = NuScenes(version='v1.0-mini', dataroot='./datasets/nuscenes', verbose=True)
nusc.list_scenes()

# Examine a scene
scene = nusc.scene[0]
first_sample = scene['first_sample_token']
nusc.render_sample(first_sample)
