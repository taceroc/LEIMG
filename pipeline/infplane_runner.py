from pathlib import Path
from simulations.base import SimulationContext
from simulations.infplane.config import InfPlaneConfig
from simulations.infplane.simulator import InfPlaneSimulator
from io_le import writers
from io_le import create_fits as write_fits
from dataclasses import asdict
import logging
logger = logging.getLogger(__name__)

def run_infplane_entry(cfg: InfPlaneConfig, run_id: int, output_root: Path, save: bool) -> dict:
    logger.info('doing runid: %s', run_id)
    sim = InfPlaneSimulator(cfg, SimulationContext(run_id=run_id, save=save))
    result = sim.run()
    if result is None:
        logger.warning('Moving to next yml instance')
        return {"run_id": run_id, "outputs": {}, "meta": {}}
    result.metadata.update({"params": asdict(cfg)}) 
    outputs = {}
    if save:
        outputs = writers.ResultWriter(f'{output_root}/{run_id}').write_arrays_and_metadata(result, cfg)
        writers.ResultWriter(f'{output_root}/{run_id}').plot_save(result, cfg)

    # optional immediate FITS write from in-memory arrays
    write_fits.main(result.surface_image, 
                    result.x_img_arc, result.y_img_arc, f'{output_root}/{run_id}')
    return {"run_id": run_id, "outputs": outputs, "meta": result.metadata}