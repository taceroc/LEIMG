from abc import ABC, abstractmethod
from dataclasses import dataclass
from .infplane.models import SimulationResult, PhaseResult

@dataclass
class SimulationContext:
    run_id: int
    save: bool

class BaseSimulation(ABC):
    @abstractmethod
    def run(self) -> "SimulationResult":
        pass