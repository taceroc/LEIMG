from dataclasses import dataclass
from typing import List, Literal
DustEnv = Literal["mw", "lmc"]
Composition = Literal["both", "S", "C"]
from utils import fix_constants as fc

@dataclass
class InfPlaneConfig:
    dt0_years: float
    so_d_ly: float
    dz0_ly: float
    ct_years: float
    a: float
    ay: float
    az: float
    z0ly: float     # [a, ay, az, z0_pc]
    angles_deg: List[float]             # [ini, end]
    wavel: float
    dust_env: DustEnv
    composition_s_c: Composition = "both"
    # extra_name: str = ""
    # loc_to_fits_subdir: str = ""

    @classmethod
    def from_yaml_entry(cls, parameters: dict) -> "InfPlaneConfig":

        required = ["dt0", "d", "dz0", "ct", "plane_coefficients", 
                    "angles", "wave", "dust_env"]
        missing = [k for k in required if k not in parameters]
        if missing:
            raise ValueError(f"Missing required keys: {missing}")

        plane = parameters["plane_coefficients"]
        if not (isinstance(plane, (list, tuple)) and len(plane) == 4):
            raise ValueError("plane_coefficients must be [a, ay, az, z0_pc]")

        angles = parameters["angles"]
        if not (isinstance(angles, (list, tuple)) and len(angles) == 2):
            raise ValueError("angles must be [ini_deg, end_deg]")
        
        ini_angle = parameters['angles'][0]
        end_angle = parameters['angles'][1]

        cfg = cls(
        dt0_years = parameters['dt0'] * fc.dtoy, #years
        so_d_ly = parameters['d'] * fc.pctoly,#parameters['d']
        dz0_ly = parameters['dz0'] * fc.pctoly,
        ct_years = parameters['ct'] * fc.dtoy,#in y
        a = parameters['plane_coefficients'][0],
        ay = parameters['plane_coefficients'][1],
        az = parameters['plane_coefficients'][2] ,
        z0ly = parameters['plane_coefficients'][3] * fc.pctoly,
        angles_deg = [ini_angle, end_angle],
        wavel = parameters['wave'], 
        dust_env = parameters['dust_env'],
        composition_s_c= parameters['composition'],
        )
        # self.bool_save = args[0]
        # self.bool_show_plots = args[1]
    
        return cfg
    

    def validate(self) -> None:
        if self.so_d_ly <= 0:
            raise ValueError("d must be > 0")
        if self.dz0_ly <= 0:
            raise ValueError("dz0 must be > 0")
        if self.wavel <= 0:
            raise ValueError("wave must be > 0")
        if self.dust_env not in ("mw", "lmc"):
            raise ValueError("dust_env must be 'mw' or 'lmc'")
        if self.composition_s_c not in ("both", "S", "C"):
            raise ValueError("composition must be 'both', 'S', or 'C'")