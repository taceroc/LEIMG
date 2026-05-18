from dataclasses import dataclass
import numpy as np
from typing import Dict, List, Optional
from .LE_img import LEImageAnalytical

@dataclass
class PhaseResult:
    tt_years: float
    x_inter: np.ndarray
    y_inter: np.ndarray
    z_inter: np.ndarray
    surface: np.ndarray
    r_le_in: np.ndarray
    r_le_out: np.ndarray
    le_img: LEImageAnalytical
    x_known: np.ndarray
    y_known: np.ndarray
    act_arc: float
    bct_arc: float
    r_le_in_arc: np.ndarray
    r_le_out_arc: np.ndarray


# tt_years=float(tt_years),
#                 x_inter=x_inter_values,
#                 y_inter=y_inter_values,
#                 z_inter=z_inter_values,
#                 surface=surface,
#                 r_le_in=LE_plane1source1_tt.r_le_in,
#                 r_le_out=LE_plane1source1_tt.r_le_out,
#                 le_img=le_img,
#                 x_known=x_known,
#                 y_known=y_known,
#                 act_arc = le_img.act,
#                 bct_arc = le_img.bct,
#                 r_le_in_arc = le_img.r_le_in,
#                 r_le_out_arc = le_img.r_le_out

@dataclass
class SimulationResult:
    surface_image: np.ndarray
    x_img_arc: np.ndarray
    y_img_arc: np.ndarray
    z_img_ly: np.ndarray
    all_surface: np.ndarray
    all_x_inter_values: np.ndarray
    all_y_inter_values: np.ndarray
    all_z_inter_values: np.ndarray
    x_known: np.ndarray
    y_known: np.ndarray
    act_all: np.ndarray
    bct_all: np.ndarray
    all_r_le_in_arc: np.ndarray
    all_r_le_out_arc: np.ndarray
    metadata: Dict[str, object]
    outputs: Optional[Dict[str, str]] = None  # populated by writer