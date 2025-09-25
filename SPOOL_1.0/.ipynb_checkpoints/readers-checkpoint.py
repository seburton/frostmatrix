from obspy import read, UTCDateTime
import numpy as np
from scipy import signal

def reader(file, start_time, end_time):
    """
    Reads SAC file and returns the portion of data within the given time window.
    Handles partial overlaps where the requested range may extend outside the file.

    PARAMS:
        file       (str)        : File path to SAC file
        start_time (np.datetime64): Requested start time
        end_time   (np.datetime64): Requested end time

    RETURNS:
        np.ndarray : Trace data for overlapping window (may be empty if no overlap)
    """
    trace = read(file)[0]
    sr = trace.stats.sampling_rate
    file_start = trace.stats.starttime
    file_end = trace.stats.endtime

    # Convert np.datetime64 inputs to UTCDateTime
    start_time = UTCDateTime(str(start_time))
    end_time = UTCDateTime(str(end_time))

    # If the requested window is completely outside this file, return empty array
    if end_time <= file_start or start_time >= file_end:
        return np.array([])

    # Clamp the requested window to the file's available range
    clipped_start = max(start_time, file_start)
    clipped_end = min(end_time, file_end)

    # Convert to indices relative to the trace start
    start_idx = int(round((clipped_start - file_start) * sr))
    end_idx = int(round((clipped_end - file_start) * sr))

    # Handle edge cases explicitly to avoid slicing errors
    start_idx = max(0, start_idx)
    end_idx = min(len(trace.data), end_idx)

    data_seg = trace.data[start_idx:end_idx]

    # Apply BP filter [1 10] Hz
    nyquist = 0.5 * sr
    low = 1 / nyquist
    high = 10 / nyquist
    b, a = signal.butter(4, [low, high], btype='band')
    data_seg = signal.filtfilt(b, a, data_seg)

    # apply demean
    data_seg = data_seg - np.mean(data_seg)

    # sensitivity correction
    sensitivity = 6.02383e+08  # V/m/s
    data_seg = data_seg / sensitivity
    
    return data_seg
