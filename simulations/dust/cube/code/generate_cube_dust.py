import numpy as np
import os
import matplotlib.pyplot as plt
# os.path.dirname(os.path.abspath(__file__)) 
from definitions import PATH_TO_DUST_IMG




def generate_cube_dust():
    dust_cube_test_spitzer = np.load(PATH_TO_DUST_IMG)
    sizes = dust_cube_test_spitzer.shape
    # if len(sizes) == 1:
    #     sizes = np.array(int())
    print(sizes)
    dust_cube_test_spitzer = dust_cube_test_spitzer.reshape(sizes[0],sizes[1])
    third = np.zeros((sizes[0]*sizes[1],sizes[0])).astype(bool)
    for ii,b in enumerate(dust_cube_test_spitzer.flatten()):
        third[ii, :int(b)] = True
    img = third.reshape(sizes[0],sizes[1],sizes[0])
    # plt.imshow(img)
    return img


def generate_cube_dust_random():
    dust_cube_test_spitzer = np.load(PATH_TO_DUST_IMG)
    sizes = dust_cube_test_spitzer.shape
    print(sizes)
    dust_cube_test_spitzer = dust_cube_test_spitzer.reshape(sizes[0],sizes[1])
    third = np.zeros((sizes[0]*sizes[1],sizes[0])).astype(bool)
    indexes = np.arange(0, 44)
    np.random.seed(32)
    for ii,b in enumerate(dust_cube_test_spitzer.flatten()):
        ij = np.random.choice(indexes, size=int(b))
        third[ii, ij] = True
    img = third.reshape(sizes[0],sizes[1],sizes[0])
    return img

def generate_cube_dust_nonbool():
    dust_cube_test_spitzer = np.load(PATH_TO_DUST_IMG)
    sizes = dust_cube_test_spitzer.shape
    dust_cube_test_spitzer = dust_cube_test_spitzer.reshape(sizes[0],sizes[1])
    third = np.zeros((sizes[0]*sizes[1],sizes[0]))
    for ii,b in enumerate(dust_cube_test_spitzer.flatten()):
        third[ii, :int(b)] = b
    img = third.reshape(sizes[0],sizes[1],sizes[0])
    return img