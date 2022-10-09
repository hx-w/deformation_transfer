# -*- coding: utf-8 -*-
'''
show obj file:
    - static/source_mesh/reference.obj
    - static/target_mesh/reference.obj
with markers on.
'''

import trimesh
import numpy as np

def read_markers(filename: str) -> np.ndarray:
    '''
    read markers from file
    '''
    markers = np.loadtxt(filename, delimiter=':', dtype=np.int32)
    return markers

def test_show_markers():
    # load source and target mesh
    source = trimesh.load('../static/source_mesh/reference.obj')
    target = trimesh.load('../static/target_mesh/reference.obj')

    markers = read_markers('../static/markers.txt')

    # colorize markers
    source.visual.vertex_colors[markers[:, 0]] = [255, 0, 0, 255]
    target.visual.vertex_colors[markers[:, 1]] = [255, 0, 0, 255]

    source.show()
    target.show()


if __name__ == '__main__':
    test_show_markers()
