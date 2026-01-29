import io
import logging
import os
import subprocess
import threading
from collections.abc import Mapping
from typing import Optional, Union


class LiveProcessResult(subprocess.CompletedProcess):
    """
    A CompletedProcess-compatible result object with real-time captured output.
    """

    def __init__(self, args, returncode, stdout: str, stderr: str):
        super().__init__(args=args, returncode=returncode, stdout=stdout, stderr=stderr)


def run_command(
    cmd: Union[tuple[str], list[str]],
    verbose: bool = False,
    env: Optional[Mapping[str, str]] = None,
) -> subprocess.CompletedProcess[str]:
    """
    Execute a command with real-time output streaming, while capturing stdout and stderr.

    Parameters:
    - cmd: List or tuple of command and arguments.
    - verbose: If True, prints output in real time and logs errors.
    - env: Optional environment variables to pass to the subprocess.

    Returns:
    - A subprocess.CompletedProcess-compatible object with .args, .returncode, .stdout, .stderr.

    Raises:
    - RuntimeError: if the command fails and verbose is True.
    """
    if verbose:
        logging.info(f"Launching command: {' '.join(cmd)}")

    # Clone and patch environment for macOS if needed
    patched_env = os.environ.copy()
    if env:
        patched_env.update(env)

    stdout_lines: list[str] = []
    stderr_lines: list[str] = []

    def stream_reader(pipe: io.IOBase, collector: list[str], label: str):
        for line in iter(pipe.readline, ""):
            if verbose:
                logging.info(f"[{label}] {line.rstrip()}")
            collector.append(line)
        pipe.close()

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
        env=patched_env,
    )

    t1 = threading.Thread(target=stream_reader, args=(process.stdout, stdout_lines, "STDOUT"))
    t2 = threading.Thread(target=stream_reader, args=(process.stderr, stderr_lines, "STDERR"))
    t1.start()
    t2.start()

    process.wait()
    t1.join()
    t2.join()

    stdout_text = "".join(stdout_lines)
    stderr_text = "".join(stderr_lines)

    if process.returncode != 0 and verbose:
        raise RuntimeError(f"--> Command failed:\n{'-' * 79}\n{stderr_text.strip()}\n{'-' * 79}")

    return LiveProcessResult(
        args=cmd,
        returncode=process.returncode,
        stdout=stdout_text,
        stderr=stderr_text,
    )
