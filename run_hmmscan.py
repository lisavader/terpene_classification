import os
from subprocess import Popen, PIPE, TimeoutExpired
from io import StringIO
from Bio import SearchIO
from typing import Any, Dict, IO, List, Optional, Union

# These functions are adapted from AntiSMASH v.7.0.1 (antismash/antismash/common/subprocessing/base.py, antismash/antismash/common/subprocessing/hmmscan.py)

class RunResult:
    """ A container for simplifying the results of running a command """
    def __init__(self, command: List[str], stdout: bytes, stderr: bytes,
                 return_code: int, piped_out: bool, piped_err: bool) -> None:
        self.command = command
        self.stdout_piped = piped_out
        self.stderr_piped = piped_err
        if piped_out:
            self.stdout = stdout.decode()
        if piped_err:
            self.stderr = stderr.decode()
        self.return_code = return_code

    def __getattribute__(self, attr: str) -> Union[bool, str, int]:
        if attr == 'stdout' and not self.stdout_piped:
            raise ValueError("stdout was redirected to file, unable to access")
        if attr == 'stderr' and not self.stderr_piped:
            raise ValueError("stderr was redirected to file, unable to access")
        return super().__getattribute__(attr)

    def successful(self) -> bool:
        """ Returns True if the command exited with an exit code of zero """
        return not self.return_code

    def get_command_string(self) -> str:
        """ Returns the command that was run to obtain this result """
        return " ".join(self.command)

def execute(commands: List[str], stdin: Optional[str] = None, stdout: Union[int, IO[Any], None] = PIPE,
            stderr: Union[int, IO[Any], None] = PIPE, timeout: int = None,
            environment_overrides: Dict[str, str] = None) -> RunResult:
    """ Executes commands in a system-independent manner via a child process.

        By default, both stderr and stdout will be piped and the outputs
        accessible.

        Arguments:
            commands: a list of arguments to execute
            stdin: None or input to be piped into the child process
            stdout: if a file is provided, stdout from the child process
                    will be piped to that file instead of the parent process
            stderr: if a file is provided, stderr from the child process
                    will be piped to that file instead of the parent process
            timeout: if provided, the child process will be terminated after
                     this many seconds
            environment_overrides: if given, the specified environment variables
                                   will be overriden with the supplied values

        Returns:
            a RunResult object containing any piped output
    """

    if stdin is not None:
        stdin_redir: Optional[int] = PIPE
        input_bytes: Optional[bytes] = stdin.encode("utf-8")
    else:
        stdin_redir = None
        input_bytes = None

    env = os.environ.copy()
    if environment_overrides:
        env.update(environment_overrides)

    with Popen(commands, stdin=stdin_redir, stdout=stdout, stderr=stderr, env=env) as proc:
        try:
            out, err = proc.communicate(input=input_bytes, timeout=timeout)
        except TimeoutExpired:
            proc.kill()
            assert isinstance(timeout, int)
            raise RuntimeError(f"Child process {commands} timed out after {timeout} seconds")

        return RunResult(commands, out, err, proc.returncode, stdout == PIPE,
                         stderr == PIPE)

def _find_error(output: list[str]) -> str:
    """ Returns the most descriptive line in error output from hmmscan """
    # is there a line that explicitly starts with the error logging?
    for i, line in enumerate(output):
        if line.startswith("Error:"):
            if i + 1 < len(output):
                return f"{line.strip()} {output[i + 1].strip()}"
            return line.strip()
    # if not, take the first non-empty line
    for line in output:
        line = line.strip()
        if line:
            return line
    # in the worst case, return a default
    return "unknown error"

def run_hmmscan(target_hmmfile: str, fastafile: str, opts: List[str] = None,
                results_file: str = None) -> SearchIO._model.query.QueryResult:
    """ Runs hmmscan on the inputs and return a list of QueryResults

        Arguments:
            target_hmmfile: the path to a HMM file to use in scanning
            fastafile: fasta file containing input sequences
            opts: a list of extra arguments to pass to hmmscan, or None
            results_file: a path to keep a copy of hmmscan results in, if provided

        Returns:
            a list of QueryResults as parsed from hmmscan output by SearchIO

    """

    command = ["hmmscan"]
    if opts is not None:
        command.extend(opts)
    command.extend([target_hmmfile, fastafile])

    result = execute(command)
    if not result.successful():
        raise RuntimeError("".join([
            f"hmmscan returned {result.return_code}: ",
            f"'{_find_error((result.stderr or result.stdout).splitlines())}'",
        ]))
    if results_file is not None:
        with open(results_file, "w", encoding="utf-8") as handle:
            handle.write(result.stdout)

    return SearchIO.parse(StringIO(result.stdout), 'hmmer3-text')