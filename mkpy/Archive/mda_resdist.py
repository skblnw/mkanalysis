import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import os

# MAIN

if __name__ == '__main__':

    os.mkdir('resdist')

    file_list = [ 
        '3.xtc',
        '4.xtc',
        '5.xtc',
        '6.xtc',
        '7.xtc',
        '8.xtc'
        ]
    offset_nframe = 100
    offset_res = 193

    ulist = [ mda.Universe('stripped.pdb',trajectory) for trajectory in file_list]
    outdir = 'resdist'

    for ii, u in enumerate(ulist):
        PROA = u.select_atoms('resid 4-197')
        PROB = u.select_atoms('resid 4-197')

        PROA_com = PROA.center_of_mass(compound='residues')
        PROB_com = PROB.center_of_mass(compound='residues')

        n_PROA = len(PROA_com)
        n_PROB = len(PROB_com)

        print('The trajectory has {} frames'.format(u.trajectory.n_frames))
        print('PROA has {} residues and PROB has {} residues'.format(n_PROA, n_PROB))


        res_dist0 = distances.distance_array(PROA_com, PROB_com,
                                        box=u.dimensions)
        res_dist_md = res_dist0

        for ts in u.trajectory[offset_nframe:]:
            PROA_com = PROA.center_of_mass(compound='residues')
            PROB_com = PROB.center_of_mass(compound='residues')

            res_dist = distances.distance_array(PROA_com, PROB_com,
                                            box=u.dimensions)
            res_dist_md += res_dist

        res_dist_md /= u.trajectory.n_frames - offset_nframe

        fig, ax = plt.subplots()
        im2 = ax.imshow(res_dist_md, origin='upper', cmap='bone')

        data=np.loadtxt('../pairs_draw_cov70')
        data[:,1:2]=data[:,1:2]+offset_res
        data=data[:10,:2].transpose()
        ax.plot(data[0],data[1],'g*')
        data=np.loadtxt('../pairs_draw_cov90')
        data[:,1:2]=data[:,1:2]+offset_res
        data=data[:10,:2].transpose()
        ax.plot(data[0],data[1],'r*')


        # add residue ID labels to axes
        tick_interval = 30
        ax.set_yticks(np.arange(n_PROA)[::tick_interval])
        ax.set_xticks(np.arange(n_PROB)[::tick_interval])
        ax.set_yticklabels(PROA.residues.resids[::tick_interval])
        ax.set_xticklabels(PROB.residues.resids[::tick_interval])

        # Ranges
        plt.xlim(0, 193)
        plt.ylim(194, 387)

        # add figure labels and titles
        plt.ylabel('Subunit B')
        plt.xlabel('Subunit A')
        plt.title('Distance between center-of-mass')

        # colorbar
        cbar2 = fig.colorbar(im2)
        cbar2.ax.set_ylabel('Distance (Angstrom)')
        im2.set_clim(0,25)

        plt.savefig('resdist/'+str(ii+3)+'.png')
