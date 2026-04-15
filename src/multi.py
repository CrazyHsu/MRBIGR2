#!/usr/bin/env python3
"""
MultiProcessMCP - Parallel Processing Module
Pure Python implementation - No R dependencies

Provides utilities for parallel task execution.
"""
import multiprocessing as mp
import subprocess
import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from typing import Callable, List, Any, Optional
import warnings

warnings.filterwarnings("ignore")


def _safe_call(func: Callable, item: Any, error_value: Any = None) -> Any:
    try:
        return func(item)
    except Exception:
        return error_value


# ========== Basic Process Control ==========

def run_cmd(cmd: str) -> int:
    """Run a shell command using subprocess.
    
    Args:
        cmd: Shell command to execute
    
    Returns:
        Return code (0 for success)
    """
    null_f = open(os.devnull, 'w')
    result = subprocess.call(cmd, shell=True, stdout=null_f, stderr=null_f)
    null_f.close()
    return result


def run_cmd_with_output(cmd: str) -> tuple:
    """Run a shell command and capture output.
    
    Args:
        cmd: Shell command to execute
    
    Returns:
        Tuple of (return_code, stdout, stderr)
    """
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.returncode, result.stdout, result.stderr


# ========== Parallel Execution ==========

def parallel_run(func: Callable, values: List[tuple], num_threads: int = 1) -> List[Any]:
    """Parallel execution using multiprocessing Pool.
    
    Args:
        func: Function to execute
        values: List of argument tuples for func
        num_threads: Number of parallel processes
    
    Returns:
        List of results
    """
    try:
        with mp.Pool(processes=num_threads) as p:
            results = [p.apply_async(func, v) for v in values]
            p.close()
            p.join()
            return [r.get() for r in results]
    except KeyboardInterrupt:
        p.terminate()
        p.join()
        raise


def parallel_map(func: Callable, items: List[Any], num_threads: int = 1, 
                 use_threads: bool = False) -> List[Any]:
    """Parallel map using ProcessPoolExecutor or ThreadPoolExecutor.
    
    Args:
        func: Function to apply to each item
        items: List of items to process
        num_threads: Number of parallel workers
        use_threads: Use threads instead of processes
    
    Returns:
        List of results in same order
    """
    executor_class = ThreadPoolExecutor if use_threads else ProcessPoolExecutor
    
    with executor_class(max_workers=num_threads) as executor:
        futures = [executor.submit(func, item) for item in items]
        results = [future.result() for future in as_completed(futures)]
    
    # Restore original order
    item_map = {id(futures[i]): i for i in range(len(futures))}
    ordered_results = [None] * len(items)
    for future in futures:
        idx = item_map[id(future)]
        ordered_results[idx] = future.result()
    
    return ordered_results


def parallel_starmap(func: Callable, args_list: List[tuple], num_threads: int = 1,
                      use_threads: bool = False) -> List[Any]:
    """Parallel starmap - func(*args) for each args in list.
    
    Args:
        func: Function to execute
        args_list: List of argument tuples
        num_threads: Number of parallel workers
        use_threads: Use threads instead of processes
    
    Returns:
        List of results
    """
    executor_class = ThreadPoolExecutor if use_threads else ProcessPoolExecutor
    
    with executor_class(max_workers=num_threads) as executor:
        futures = [executor.submit(func, *args) for args in args_list]
        results = [future.result() for future in futures]
    
    return results


# ========== Batch Processing ==========

def batch_process(items: List[Any], batch_size: int, func: Callable = None, 
                  num_threads: int = 1) -> List[Any]:
    """Process items in batches.
    
    Two modes:
    1. With func: Apply func to each batch
    2. Without func: Just split into batches
    
    Args:
        items: List of items to process
        batch_size: Size of each batch
        func: Optional function to apply to each batch
        num_threads: Number of parallel workers
    
    Returns:
        List of batch results or list of batches
    """
    batches = [items[i:i+batch_size] for i in range(0, len(items), batch_size)]
    
    if func is None:
        return batches
    
    return parallel_map(func, batches, num_threads=num_threads)


def batch_split(items: List[Any], batch_size: int) -> List[List[Any]]:
    """Split items into batches without applying any function.
    
    Args:
        items: List of items to split
        batch_size: Size of each batch
    
    Returns:
        List of batches
    """
    return [items[i:i+batch_size] for i in range(0, len(items), batch_size)]


def chunk_process(items: List[Any], chunk_size: int, func: Callable,
                  num_threads: int = 1) -> List[Any]:
    """Process items in chunks with result aggregation.
    
    Args:
        items: List of items to process
        chunk_size: Size of each chunk
        func: Function to apply to each chunk (should return list)
        num_threads: Number of parallel workers
    
    Returns:
        Combined list of all results
    """
    chunks = [items[i:i+chunk_size] for i in range(0, len(items), chunk_size)]
    chunk_results = parallel_map(func, chunks, num_threads=num_threads)
    
    # Flatten results
    combined = []
    for result in chunk_results:
        if isinstance(result, list):
            combined.extend(result)
        else:
            combined.append(result)
    
    return combined


# ========== Progress Tracking ==========

def parallel_map_with_progress(func: Callable, items: List[Any], 
                               num_threads: int = 1,
                               progress_callback: Optional[Callable] = None) -> List[Any]:
    """Parallel map with progress tracking.
    
    Args:
        func: Function to apply to each item
        items: List of items to process
        num_threads: Number of parallel workers
        progress_callback: Callback function called with (completed, total)
    
    Returns:
        List of results
    """
    results = [None] * len(items)
    completed = 0
    
    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = {executor.submit(func, item): i for i, item in enumerate(items)}
        
        for future in as_completed(futures):
            idx = futures[future]
            results[idx] = future.result()
            completed += 1
            
            if progress_callback:
                progress_callback(completed, len(items))
    
    return results


# ========== Error Handling ==========

def parallel_map_safe(func: Callable, items: List[Any], 
                      num_threads: int = 1,
                      error_value: Any = None) -> List[Any]:
    """Parallel map with error handling - returns error_value on failure.
    
    Args:
        func: Function to apply to each item
        items: List of items to process
        num_threads: Number of parallel workers
        error_value: Value to return on error
    
    Returns:
        List of results (or error_value on failure)
    """
    executor_class = ProcessPoolExecutor

    with executor_class(max_workers=num_threads) as executor:
        futures = [executor.submit(_safe_call, func, item, error_value) for item in items]
        return [future.result() for future in futures]


# ========== Utility Functions ==========

def get_optimal_threads(max_threads: Optional[int] = None) -> int:
    """Get optimal number of threads based on CPU count.
    
    Args:
        max_threads: Maximum threads to use (None = auto)
    
    Returns:
        Number of threads to use
    """
    cpu_count = mp.cpu_count()
    if max_threads is None:
        return cpu_count
    return min(cpu_count, max_threads)


def set_affinity(cpu_ids: List[int]) -> bool:
    """Set CPU affinity for current process (Linux only).
    
    Args:
        cpu_ids: List of CPU IDs to use
    
    Returns:
        True if successful
    """
    try:
        os.sched_setaffinity(0, set(cpu_ids))
        return True
    except (AttributeError, OSError):
        return False


__all__ = [
    'run_cmd', 'run_cmd_with_output',
    'parallel_run', 'parallel_map', 'parallel_starmap',
    'batch_process', 'chunk_process',
    'parallel_map_with_progress',
    'parallel_map_safe',
    'get_optimal_threads', 'set_affinity'
]
