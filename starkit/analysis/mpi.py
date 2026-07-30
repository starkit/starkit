from __future__ import annotations

import logging
import tempfile
from dataclasses import dataclass
from pathlib import Path

@dataclass(frozen=True)
class RuntimeContext:
    rank: int
    size: int
    run_dir: Path
    enabled: bool

    @property
    def is_root(self) -> bool:
        return self.rank == 0


def _get_mpi():
    try:
        from mpi4py import MPI
    except ImportError:  # pragma: no cover
        return None
    return MPI


def build_runtime_context(enabled: bool, requested_run_dir: Path | None) -> RuntimeContext:
    mpi = _get_mpi() if enabled else None
    if not enabled or mpi is None:
        run_dir = requested_run_dir or Path(tempfile.mkdtemp(prefix="starkit_analysis_"))
        run_dir.mkdir(parents=True, exist_ok=True)
        return RuntimeContext(rank=0, size=1, run_dir=run_dir, enabled=False)

    comm = mpi.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    if rank == 0:
        if requested_run_dir is None:
            run_dir = Path(tempfile.mkdtemp(prefix="starkit_analysis_"))
        else:
            run_dir = requested_run_dir
            run_dir.mkdir(parents=True, exist_ok=True)
    else:
        run_dir = None
    run_dir = comm.bcast(run_dir, root=0)
    return RuntimeContext(rank=rank, size=size, run_dir=run_dir, enabled=True)


class RootRankFilter(logging.Filter):
    def __init__(self, rank: int) -> None:
        super().__init__()
        self.rank = rank

    def filter(self, record: logging.LogRecord) -> bool:
        return self.rank == 0


def configure_logging(rank: int) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    rank_filter = RootRankFilter(rank)
    root_logger = logging.getLogger()
    for handler in root_logger.handlers:
        handler.addFilter(rank_filter)
