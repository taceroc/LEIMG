import logging
import numpy as np
from simulations.base import BaseSimulation, SimulationContext
from .config import InfPlaneConfig
from .models import SimulationResult, PhaseResult
from .surface import SurfaceAssembler
# from . import phase_ops
from dataclasses import dataclass
from .DustShape import InfPlane
import simulations.infplane.LightEcho as LE
from .Source import Source
import simulations.infplane.SurfaceBrightness as sb
from .LE_img import LEImageAnalytical
from utils import fix_constants as fc
from utils import extract_lc
import pandas as pd
from collections.abc import Iterable
import sys
logger = logging.getLogger(__name__)

@dataclass
class LightCurve:
    time_years: np.ndarray
    mag: np.ndarray

@dataclass
class InfPlaneState:
    cfg: InfPlaneConfig
    lc: LightCurve
    source: Source
    plane: InfPlane
    too_big: bool
    x_bins: np.array
    y_bins: np.array
    ranges: list
    pixel_resolution: float
    act_all: list #center for phase 0
    bct_all: list
    z_all_ly: np.array #z for all
    surface_known: list
   
    # add anything phase ops need (constants, masks, cached objects)

class InfPlaneSimulator(BaseSimulation):
    def __init__(self, config: InfPlaneConfig, context: SimulationContext, pixel_resolution: float = 0.2):
        self.not_use = 0
        self.config = config
        self.context = context
        self.pixel_resolution = pixel_resolution
        self.too_big = False
        self.max_size_img = 3000
        self.min_size_img = 500
        self.results_simulation = SimulationResult(surface_image = np.array([]),
                            x_img_arc= np.array([]),
                            y_img_arc= np.array([]),
                            z_img_ly= np.array([]),
                            all_surface= np.array([]),
                            all_x_inter_values= np.array([]),
                            all_y_inter_values= np.array([]),
                            all_z_inter_values= np.array([]),
                            x_known= np.array([]),
                            y_known= np.array([]),
                            act_all= np.array([]),
                            bct_all= np.array([]),
                            all_r_le_in_arc=np.array([]),
                            all_r_le_out_arc=np.array([]),
                            metadata = {}
                            )

        self.state = self._build_state(config)
        # init source/dust/light-curve/state objects


    def _build_state(self, config):
        lc = self._load_lightcurve(file_path="data/lightcurves/LC_sn2011fe.csv")
        source = Source(config.dt0_years, config.so_d_ly, Flmax=0)
        plane = InfPlane([config.a, config.ay, config.az, -config.z0ly], config.dz0_ly)
 
        return InfPlaneState(cfg=config, lc=lc, source=source, plane=plane, 
                             pixel_resolution=self.pixel_resolution, too_big=self.too_big,
                            x_bins=np.array([]), y_bins=np.array([]), ranges=[], act_all=[], bct_all=[], z_all_ly=np.array([]), surface_known=[])
        
    def _load_lightcurve(self, file_path):
        file_path = '/pscratch/sd/t/taceroc/LE_experiments/LE_pkg/data/lightcurves/LC_sn2011fe.csv'
        lc_sn2011fe = extract_lc.read_from_file(file_path)
        lc = {}
        lc['mag'] = lc_sn2011fe['mag'].values - 5
        lc['time'] = ((lc_sn2011fe['time'].values) - lc_sn2011fe['time'].min() )* fc.dtoy
        lc['time'] = lc['time'] - lc['time'][lc['mag'] == lc['mag'].min()]
        return LightCurve(time_years=lc['time'], mag=lc['mag'])

    def calculate_sb(self, tt_years, x_inter_values, y_inter_values, z_inter_values, LE_plane1source1_tt):
        """
            Calculate the surface brightness for the phase
            Arguments:
                tt: phase in years
                x_inter_values, y_inter_values, z_inter_values: values in ly for intersection
                LE_plane1source1_tt: LE_planesource object (LightEcho.LEPlane())
            Return:
                surface: surface brightness in the units ly,cm,s
        """
        # Calculate surface brightness
        ctt = self.state.cfg.ct_years - tt_years
        sb_plane_all = sb.SurfaceBrightnessAnalytical(self.state.cfg.wavel, self.state.source, 
                                                      LE_plane1source1_tt, [x_inter_values, y_inter_values, z_inter_values], 
                                                      self.state.cfg.composition_s_c, self.state.cfg.dust_env)
        
        sb_plane_all.lc = pd.DataFrame({
                    'time': self.state.lc.time_years,
                    'mag': self.state.lc.mag
                })

        cossigma, surface = sb_plane_all.calculate_surface_brightness(tt_years)

        return cossigma, surface
        
    def _calculate_phase_0(self):
        x_inter_values, y_inter_values, z_inter_values, new_xs, new_ys, new_zs, LE_plane1source1_tt1 = self._calculate_xyz_phase_lc_time(self.state.lc.time_years[::40][0],
                                                                                                                                          return_LE_plane1source1=True)

        if np.mean(np.sqrt([x_inter_values**2 + y_inter_values**2]))/self.state.cfg.so_d_ly > 1e-2:
            logger.critical("Elliposoid -> Paraboloid approximation no longer true %s", 
                           np.mean(np.sqrt([x_inter_values**2 + y_inter_values**2]))/self.state.cfg.so_d_ly)
            return None #sys.exit(1)
        else:
            logger.info('approximation holds %s', np.mean(np.sqrt([x_inter_values**2 + y_inter_values**2]))/self.state.cfg.so_d_ly)            
            logger.info('%s', (np.mean(np.sqrt([x_inter_values**2 + y_inter_values**2])), self.state.cfg.so_d_ly))
            

            
        _, surface = self.calculate_sb(self.state.lc.time_years[::40][0], x_inter_values, y_inter_values, z_inter_values, LE_plane1source1_tt1)
        
        le_img1 = LEImageAnalytical(LE_plane1source1_tt1, self.state.plane, surface, pixel_resolution=self.pixel_resolution)
        logger.info("rin and rout for phase 0%s", (le_img1.r_le_in[0], le_img1.r_le_out[0]))
        return le_img1

    def _theta_from_start_span(self, start_deg, span_deg):
        start = float(start_deg) % 360
        span = max(float(span_deg), 0.5)
        end = start + span
        return [np.deg2rad(start), np.deg2rad(end)]

    def _prepare_image_grid(self):
        le_img1 = self._calculate_phase_0()
        if le_img1 is None:
            return None
        x_size_img, y_size_img = le_img1._estimate_image_size()
        logger.info("estimated img pixels %s", (x_size_img, y_size_img))

        # iterate to find a size that can be handled
        if max(x_size_img, y_size_img) > self.max_size_img:
            self.state.too_big = True
            ini_angle = float(self.state.cfg.angles_deg[0]) % 360
            end_angle = float(self.state.cfg.angles_deg[1])
            if end_angle < 0:
                end_norm = end_angle + 360
            else:
                end_norm = end_angle % 360
            requested_span = (end_norm - ini_angle) % 360
            if np.isclose(requested_span, 0):
                requested_span = 360.0

            if requested_span >= 359.0:
                current_span = 90.0
                logger.info("Image too large. Starting adaptive arc at 90 degrees.")
            else:
                current_span = requested_span
                logger.info(f"Image too large. Starting adaptive arc from requested span {requested_span:.2f} degrees.")

            min_span_deg = 0.5
            while True:
                self.state.theta = self._theta_from_start_span(ini_angle, current_span)
                x_inter_values, y_inter_values, z_inter_values, new_xs, new_ys, new_zs, LE_plane1source1_tt1 = self._calculate_xyz_phase_lc_time(self.state.lc.time_years[::40][0],
                                                                                                                                          return_LE_plane1source1=True)
                _, surface1 = self.calculate_sb(self.state.lc.time_years[::40][0], x_inter_values, y_inter_values, z_inter_values, LE_plane1source1_tt1)
                le_img1 = LEImageAnalytical(LE_plane1source1_tt1, self.state.plane, surface1, pixel_resolution=self.pixel_resolution)
                x_size_img, y_size_img = le_img1._estimate_image_size()
                logger.info(f"adaptive arc span {current_span:.4f} deg -> estimated img pixels {x_size_img}, {y_size_img}")

                if (max(x_size_img, y_size_img) <= self.max_size_img) or (min(x_size_img, y_size_img) <= self.min_size_img) :
                    break
                if current_span <= min_span_deg:
                    logger.info(f"Reached minimum arc span ({min_span_deg} deg). Proceeding with the smallest arc.")
                    break
                current_span = max(current_span / 2.0, min_span_deg)

        logger.info("rin and rout for phase 0 after too large%s", (le_img1.r_le_in[0], le_img1.r_le_out[0]))
        self.state.x_bins, self.state.y_bins, x_all, y_all, x_all_ly, y_all_ly, self.state.z_all_ly, range_x, range_y, range_x_p, range_y_p = le_img1.create_image_grid(x_size_img, y_size_img)
        self.state.ranges = [range_x, range_y, range_x_p, range_y_p]

        self.state.act_all.extend([le_img1.act])
        self.state.bct_all.extend([le_img1.bct])

        return True

    def _calculate_phases_geometry(self, tt_years):
        x_inter_values, y_inter_values, z_inter_values, new_xs, new_ys, new_zs, LE_plane1source1_tt = self._calculate_xyz_phase_lc_time(tt_years, return_LE_plane1source1=True)
        bool_mask_sizes = ((x_inter_values <= self.state.ranges[0][0]) & (x_inter_values >= self.state.ranges[0][1]))
        x_inter_values = x_inter_values[bool_mask_sizes]
        y_inter_values = y_inter_values[bool_mask_sizes]
        size_common_all = np.min([len(x_inter_values), len(y_inter_values)])
        x_inter_values = x_inter_values[:size_common_all]
        y_inter_values = y_inter_values[:size_common_all]
        z_inter_values = LE_plane1source1_tt.func_for_z_plane(x_inter_values, y_inter_values)
        bool_mask_sizes = ((LE_plane1source1_tt.x_projected[:,0].flatten() <= self.state.ranges[2][0]) & (LE_plane1source1_tt.x_projected[:,0].flatten() >= self.state.ranges[2][1]))
        x_project_0 = LE_plane1source1_tt.x_projected[:,0, bool_mask_sizes].flatten()
        x_project_1 = LE_plane1source1_tt.x_projected[:,1, bool_mask_sizes].flatten()
        size_common = np.min([x_project_0.shape[-1], x_project_1.shape[-1]])
        
        LE_plane1source1_tt.x_projected = np.concatenate([x_project_0[:size_common], x_project_1[:size_common]]).reshape(1,2,size_common)
        
        y_project_0 = LE_plane1source1_tt.y_projected[:,0, bool_mask_sizes].flatten()
        y_project_1 = LE_plane1source1_tt.y_projected[:,1, bool_mask_sizes].flatten()
        size_common = np.min([y_project_0.shape[-1], y_project_1.shape[-1]])
        LE_plane1source1_tt.y_projected = np.concatenate([y_project_0[:size_common], y_project_1[:size_common]]).reshape(1,2,size_common)

        size_common = np.min([LE_plane1source1_tt.y_projected.shape[-1], LE_plane1source1_tt.x_projected.shape[-1]])
        LE_plane1source1_tt.x_projected = LE_plane1source1_tt.x_projected[:,:,:size_common_all]
        LE_plane1source1_tt.y_projected = LE_plane1source1_tt.y_projected[:,:,:size_common_all]
        LE_plane1source1_tt.z_projected = LE_plane1source1_tt.func_for_z(LE_plane1source1_tt.x_projected, LE_plane1source1_tt.y_projected)

        return x_inter_values, y_inter_values, z_inter_values, LE_plane1source1_tt

    def calculate_le_2d(self, LE_plane1source1_tt, surface):
        """
            Fill the surface_val given the bins
            Arguments:
               LE_plane1source1_tt: LE_planesource object (LightEcho.LEPlane())
               surface: surface values in the units ly,cm,s
            Return:
                le_img: the LE image with each LE by phase, no interpolation here
        """
        # Initialize and calculate the LE image from LE and surface brightness
        le_img = LEImageAnalytical(LE_plane1source1_tt, self.state.plane, surface, pixel_resolution=self.pixel_resolution)
        xs_outer = le_img.new_xs[0, 0, :]
        ys_outer = le_img.new_ys[0, 0, :]
        xs_inner = le_img.new_xs[0, 1, :]
        ys_inner = le_img.new_ys[0, 1, :]

        xs_all = np.concatenate([xs_outer, xs_inner])
        ys_all = np.concatenate([ys_outer, ys_inner])

        inds_x2 = np.digitize(xs_all, self.state.x_bins)
        inds_y2 = np.digitize(ys_all, self.state.y_bins)
        inds_x2 = np.clip(inds_x2, 1, len(self.state.x_bins) - 1)
        inds_y2 = np.clip(inds_y2, 1, len(self.state.y_bins) - 1)

        x_known = self.state.x_bins[inds_x2-1]
        y_known = self.state.y_bins[inds_y2-1]
        # self.state.surface_known = np.concatenate([surface, surface])###############################
        self.state.surface_known.extend(np.concatenate([surface, surface]))###############################
        
    
        # self.surface_known.extend(surface)
    
        # self.act_all.extend([le_img.act])
        # self.bct_all.extend([le_img.bct])

        

        return le_img, x_known, y_known
    
    def _calculate_xyz_phase_lc_time(self, tt_years, return_LE_plane1source1 = False):
        """
            Calculate the intersection paraboloid and ellipsoid for a phase of the LC
            Arguments:
                tt: phase in years
                return_LE_plane1source1: bool, True to return the LE_planesource object (LightEcho.LEPlane())
            Return:
                x,y,z: values in ly for intersection
                new_xs, new_ys, new_zs: x,y,z position in the x-y plane in arcseconds 
                x_inter_values, y_inter_values, z_inter_values, new_xs, new_ys, new_zs, LE_plane1source1_tt
        """
        ctt = self.state.cfg.ct_years - tt_years
   
        if ctt <= 0:
            self.not_use = self.not_use + 1
        else:
            pass
            # self.times_lc_use.append(ctt)
            # self.times_lc.append(tt_years)
        if self.not_use == len(self.state.lc.time_years):
            logger.warning("No valid phase for ct=%s", tt_years)
            raise TypeError("no LE at this phase")
        
        logger.info("TIME LE elliposid in days=%s", ctt / fc.dtoy)
        logger.info("TIME of the SN LC respect to peak %s", tt_years/ fc.dtoy)
        
        # initiate the geometry, source, time    
        LE_plane1source1_tt = LE.LEPlane(ctt, self.state.plane, self.state.source)
        LE_plane1source1_tt.theta = [0, 2*np.pi]
        if self.state.too_big == True:
            LE_plane1source1_tt.theta = self.state.theta
            
        # calculate x,y,z values of LE equation
        x_inter_values, y_inter_values, z_inter_values, new_xs, new_ys, new_zs = LE_plane1source1_tt.run()
        

        if return_LE_plane1source1:
            return x_inter_values, y_inter_values, z_inter_values, new_xs, new_ys, new_zs, LE_plane1source1_tt
        else:
            return x_inter_values, y_inter_values, z_inter_values, new_xs, new_ys, new_zs
    
    def _simulate_phase(self) -> list:
        out = []
        for idx, tt_years in enumerate(self.state.lc.time_years[1::40]):
            x_inter_values, y_inter_values, z_inter_values, LE_plane1source1_tt = self._calculate_phases_geometry(tt_years)
            try:
                print(len(x_inter_values))
            except:
                logger.info('No len for x_inter_values skipping to next phase')
                continue
            _, surface = self.calculate_sb(tt_years, x_inter_values, y_inter_values, z_inter_values, LE_plane1source1_tt)
            le_img, x_known, y_known = self.calculate_le_2d(LE_plane1source1_tt, surface)

            out.append(PhaseResult(
                tt_years=float(tt_years),
                x_inter=x_inter_values,
                y_inter=y_inter_values,
                z_inter=z_inter_values,
                surface=surface,
                r_le_in=LE_plane1source1_tt.r_le_in,
                r_le_out=LE_plane1source1_tt.r_le_out,
                le_img=le_img,
                x_known=x_known,
                y_known=y_known,
                act_arc = le_img.act,
                bct_arc = le_img.bct,
                r_le_in_arc = le_img.r_le_in,
                r_le_out_arc = le_img.r_le_out

            ))
        return out
    
    def _assemble_surface(self, out):

        def list_convert(attr, out):
            temp_xknow = []
            for ip in out:
                if isinstance(getattr(ip, attr), Iterable):
                    temp_xknow.extend(getattr(ip, attr).flatten())
                else:
                    temp_xknow.append(getattr(ip, attr))
            return temp_xknow

        self.results_simulation.all_x_inter_values = np.array(list_convert('x_inter', out))
        print(self.results_simulation.all_x_inter_values.shape)
        self.results_simulation.all_y_inter_values = np.array(list_convert('y_inter', out))
        self.results_simulation.all_z_inter_values = np.array(list_convert('z_inter', out))
        self.results_simulation.all_surface = np.array(list_convert('surface', out)).reshape(len(self.results_simulation.all_x_inter_values))

        print("Before angle clipping")
        logger.info(f'{max(self.state.surface_known), min(self.state.surface_known)}')

        
        self.results_simulation.x_known = np.array(list_convert('x_known', out))
        self.results_simulation.y_known = np.array(list_convert('y_known', out))
        self.results_simulation.act_all = np.array(list_convert('act_arc', out))
        self.results_simulation.act_all = np.insert(self.results_simulation.act_all, [0], [self.state.act_all])
        self.results_simulation.bct_all = np.array(list_convert('bct_arc', out))
        self.results_simulation.bct_all = np.insert(self.results_simulation.bct_all, [0], [self.state.bct_all])


        self.results_simulation.all_r_le_in_arc = np.array(list_convert('r_le_in_arc', out))
        self.results_simulation.all_r_le_out_arc = np.array(list_convert('r_le_out_arc', out))

        surface_mag_arsec2 = self.results_simulation.all_surface / (0.2**2) / (3e10 / ((self.state.cfg.wavel*1e-4))**2 * (1e-5)) # c/wave^2 * deltawave, e.g 100nm
        
        logger.info(f"FLUX {np.mean(self.results_simulation.all_surface)}")
        mags_plot = -2.5*np.log10(surface_mag_arsec2[surface_mag_arsec2>0])-48.6
        mags_plot = np.nan_to_num(mags_plot, nan=0.0, posinf=0.0, neginf=0.0)
        logger.info(f"mean surface {np.mean(mags_plot)}")
        # logger.info(self.times_lc)
        # logger.info(self.times_lc_use)
        logger.info(f"{self.results_simulation.all_surface.shape}")

        return SurfaceAssembler(self.state.x_bins, self.state.y_bins, self.results_simulation.all_surface, self.state.surface_known)


    def run(self) -> SimulationResult:
        # 1) estimate image grid
        res = self._prepare_image_grid() 
        if res is None:
            return None
        self.state.surface_known = []
        # 2) loop phases -> collect PhaseResult
        results_phases = self._simulate_phase()
        # 3) assemble interpolated surface
        surf_assemble = self._assemble_surface(results_phases)
        inter_surface, grid_x, grid_y = surf_assemble.create_le_surface_interpolate(self.results_simulation.x_known, self.results_simulation.y_known)
        # 4) apply masks
        r_in = np.nanmedian(np.asarray(self.results_simulation.all_r_le_in_arc[-1], dtype=float))
        r_out = np.nanmedian(np.asarray(self.results_simulation.all_r_le_out_arc[0], dtype=float))
        center_x = np.nanmean(abs(self.results_simulation.act_all[0]-self.results_simulation.act_all[-1]))
        center_y = np.nanmean(abs(self.results_simulation.bct_all[0]-self.results_simulation.bct_all[-1]))
        (inter_surface, self.results_simulation.x_img_arc, 
         self.results_simulation.y_img_arc, self.results_simulation.z_img_ly) = surf_assemble.apply_annulus_mask(inter_surface, 
                                                                                                                 grid_x, grid_y, 
                                                                                                                 r_in, r_out,
                                                                                                                 center_x, center_y, self.state.z_all_ly)
        if self.too_big == False:
            self.results_simulation.surface_image = surf_assemble.apply_angle_mask(inter_surface,
                                                                                   x_img_arc=self.results_simulation.x_img_arc, 
                                                                                   y_img_arc=self.results_simulation.y_img_arc, 
                                                                                   ini_angle=self.state.cfg.angles_deg[0], 
                                                                                   end_angle=self.state.cfg.angles_deg[1])
        
        self.results_simulation.surface_image = self.results_simulation.surface_image / (0.2**2) / (3e10 / ((self.state.cfg.wavel*1e-4))**2 * (1e-5)) # c/wave^2 * deltawave, e.g 100nm

        self.results_simulation.metadata = {
            "run_id": self.context.run_id,
            "n_phases": len(results_phases),
            "ct_years": self.config.ct_years,
        }
        # 5) return SimulationResult (no writing here)
        return self.results_simulation

        
                
        
       
        
        ...