"""MRSIPrep command entry point."""

from __future__ import annotations

import sys

from mrsiprep.cli.parser import parse_args
from mrsiprep.utils.logging import setup_logging
from mrsiprep.workflows.participant import run_participant_workflow


def main(argv: list[str] | None = None) -> int:
    config = parse_args(argv)
    logger = setup_logging(config.verbose)
    if config.analysis_level != "participant":
        logger.error("Only participant analysis level is currently supported.")
        return 2
    statuses = run_participant_workflow(config)
    failed = [status for status in statuses if status.status != "success"]
    succeeded = [status for status in statuses if status.status == "success"]
    logger.info("MRSIPrep finished: %d succeeded, %d failed", len(succeeded), len(failed))
    for status in failed:
        logger.error("FAILED sub-%s%s: %s", status.subject, f" ses-{status.session}" if status.session else "", status.error)
    return 1 if failed and not succeeded else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
