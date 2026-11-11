"""
parallel.py
-----------
Runs the slow pairwise metrics (alignment and edit distance) across CPU cores.

Biopython's aligner holds Python's GIL, so threads do not help; a process
pool does (about 6-7x faster on a 12-core machine for nucleotide search).
Results are identical to the serial functions in metrics.py.

Functions:
    alignment_similarities(query, seqs, mode, gap_cost) -> list of floats (0.0 to 1.0)
    edit_similarities(query, seqs)                      -> list of floats (0.0 to 1.0)
"""

import os
import sys
import threading
from concurrent.futures import ProcessPoolExecutor

from metrics import get_protein_aligner, dna_aligner, edit_distance_similarity

# Below this many comparisons the pool overhead is not worth it.
MIN_PARALLEL = 8

_pool = None
_pool_lock = threading.Lock()


def _workers():
    # Leave one core for the web server; more than 8 adds little for ~300 sequences.
    return max(1, min(8, (os.cpu_count() or 2) - 1))


def _start_all_workers(pool):
    """
    Starts every worker now, without re-running the main script in each one.

    On Windows a new process normally re-imports the main script (app.py),
    which loads Flask, matplotlib, scipy and pyJASPAR: ~170 MB per worker.
    Workers only need this module, so the main script's path is hidden
    while they start; they then import just parallel.py and metrics.py.
    """
    main = sys.modules.get('__main__')
    saved = {name: main.__dict__[name] for name in ('__file__', '__spec__') if main and name in main.__dict__}
    try:
        for name in saved:
            if name == '__spec__':
                main.__spec__ = None
            else:
                del main.__dict__[name]
        # The pool starts one process per task submitted while none is idle.
        futures = [pool.submit(_dna_score, ("A", "A")) for _ in range(_workers())]
    finally:
        for name, value in saved.items():
            main.__dict__[name] = value
    for f in futures:
        f.result()


def _get_pool():
    """Creates the shared process pool, with all workers started, on first use."""
    global _pool
    with _pool_lock:
        if _pool is None:
            pool = ProcessPoolExecutor(max_workers=_workers())
            _start_all_workers(pool)
            _pool = pool
        return _pool


# Worker functions: module-level so they can be sent to other processes.

def _protein_score(task):
    query, seq, gap_cost = task
    return get_protein_aligner(*gap_cost).score(query, seq)


def _dna_score(task):
    query, seq = task
    return dna_aligner.score(query, seq)


def _edit_similarity(task):
    query, seq = task
    return edit_distance_similarity(query, seq)


def _run(fn, tasks):
    """Maps fn over tasks in the pool, or serially for small or failed runs."""
    global _pool
    if len(tasks) < MIN_PARALLEL:
        return [fn(t) for t in tasks]
    try:
        chunk = max(1, len(tasks) // (_workers() * 4))
        return list(_get_pool().map(fn, tasks, chunksize=chunk))
    except Exception as e:
        # A broken pool (e.g. a worker crashed) should not break the search.
        print(f"  [parallel] Pool failed ({e}); running serially.")
        with _pool_lock:
            _pool = None
        return [fn(t) for t in tasks]


def warm_up():
    """
    Starts the worker processes ahead of time. On Windows each worker takes
    a few seconds to start, which would otherwise delay the first search.
    """
    try:
        _get_pool()
    except Exception as e:
        print(f"  [parallel] Warm-up failed ({e}); searches will start the pool.")


def alignment_similarities(query, seqs, mode='protein', gap_cost=(11, 1)):
    """
    Normalized local alignment scores of query against each sequence:
    score / self-score of the query, clipped at 0. Same values as
    blosum_similarity (protein) and dna_similarity (DNA) in metrics.py,
    but the query's self-score is computed once instead of per sequence.
    """
    if mode == 'protein':
        aligner = get_protein_aligner(*gap_cost)
        fn, tasks = _protein_score, [(query, s, tuple(gap_cost)) for s in seqs]
    else:
        aligner = dna_aligner
        fn, tasks = _dna_score, [(query, s) for s in seqs]

    self_score = aligner.score(query, query)
    if self_score <= 0:
        return [0.0] * len(seqs)
    return [max(0.0, score / self_score) for score in _run(fn, tasks)]


def edit_similarities(query, seqs):
    """Normalized edit distance similarity of query against each sequence."""
    return _run(_edit_similarity, [(query, s) for s in seqs])
