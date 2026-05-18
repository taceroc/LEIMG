import pipeline.loader as loader
import pipeline.infplane_runner as infplane_runner
import argparse
from pathlib import Path
import yaml
import logging
import sys

def configure_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(processName)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stdout,
    )
    

def parse_args():
    p = argparse.ArgumentParser("le-runner")
    p.add_argument("funcsim", choices=["SimulateLEInfPlane"])
    p.add_argument("-file_to_parameters", required=True)          # multi-run YAML
    p.add_argument("-outdir", required=True)                      # output root
    p.add_argument("--bool_save", action=argparse.BooleanOptionalAction, default=True)
    return p.parse_args()

def main():
    args = parse_args()
    configure_logging(verbose=True)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = []
    configs, keys = loader.load_configs_from_yaml(args.file_to_parameters)
    
    for run_id, cfg in zip(keys, configs):
        rec = infplane_runner.run_infplane_entry(cfg, run_id, outdir, args.bool_save)
        manifest.append(rec)

    with open(outdir / "manifest.yml", "w") as f:
        yaml.safe_dump(manifest, f, sort_keys=False)

if __name__ == "__main__":
    main()

