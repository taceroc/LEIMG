import yaml
from simulations.infplane.config import InfPlaneConfig

def load_configs_from_yaml(yaml_path: str) -> list[InfPlaneConfig]:
    with open(yaml_path, "r") as f:
        data = yaml.safe_load(f)
    configs = []
    print(data)
    for row, run_id in zip(data, data.keys()):
        if not isinstance(data[row], dict):
            raise ValueError(f"Entry {i} is not a mapping.")
        cfg = InfPlaneConfig.from_yaml_entry(data[row])  # one entry at a time
        cfg.validate()
        configs.append(cfg)
    return configs, data.keys()