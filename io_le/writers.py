from pathlib import Path
from simulations.infplane.models import SimulationResult, PhaseResult
from simulations.infplane.config import InfPlaneConfig
import numpy as np
import os
from dataclasses import asdict
import yaml
import matplotlib.pyplot as plt


class ResultWriter:
    def __init__(self, output_root: str):
        self.outdir = output_root
        outdir = Path(os.path.join(self.outdir, 'arrays'))
        outdir.mkdir(parents=True, exist_ok=True)
        outdir = Path(os.path.join(self.outdir, 'figures'))
        outdir.mkdir(parents=True, exist_ok=True)
        outdir = Path(os.path.join(self.outdir, 'fits'))
        outdir.mkdir(parents=True, exist_ok=True)


    def write_arrays_and_metadata(self, result: SimulationResult, cfg: InfPlaneConfig):
        np.save(os.path.join(self.outdir, 'arrays/surface.npy'), result.all_surface)
        np.save(os.path.join(self.outdir, 'arrays/surface_values.npy'), result.surface_image)
        np.save(os.path.join(self.outdir, 'arrays/ximg_arcsec.npy'), result.x_img_arc)
        np.save(os.path.join(self.outdir, 'arrays/yimg_arcsec.npy'), result.y_img_arc)
        np.save(os.path.join(self.outdir, 'arrays/zimgly.npy'), result.z_img_ly)
 
        
        np.save(os.path.join(self.outdir, 'arrays/x_ly.npy'), result.all_x_inter_values)
        np.save(os.path.join(self.outdir, 'arrays/y_ly.npy'), result.all_y_inter_values)
        np.save(os.path.join(self.outdir, 'arrays/z_ly.npy'), result.all_z_inter_values)

        with open(os.path.join(self.outdir,"run_params.yml"), "w") as f:
            yaml.safe_dump(asdict(cfg), f, sort_keys=False)

        return [os.path.join(self.outdir, 'arrays/surface.npy'),
                os.path.join(self.outdir, 'arrays/surface_values.npy'),
                os.path.join(self.outdir, 'arrays/ximg_arcsec.npy'), 
                os.path.join(self.outdir, 'arrays/yimg_arcsec.npy'),
                os.path.join(self.outdir, 'arrays/zimgly.npy'),
                os.path.join(self.outdir, 'arrays/x_ly.npy'),
                os.path.join(self.outdir, 'arrays/y_ly.npy'),
                os.path.join(self.outdir, 'arrays/z_ly.npy')]
    
    def plot_save(self, result: SimulationResult, cfg: InfPlaneConfig):
        mags = -2.5*np.log10(result.surface_image)-48.6
        mags = np.nan_to_num(mags, nan=0.0, posinf=0.0, neginf=None)
        fig, ax = plt.subplots(1,1, figsize = (8,8))
        aja = ax.imshow(mags, origin = "lower", cmap="RdPu")
        plt.colorbar(aja)
        ax.set_title(f'surface image at {cfg.ct_years} years')
        plt.savefig(os.path.join(self.outdir, 'figures/surface.png'), dpi=100)
        plt.close()

