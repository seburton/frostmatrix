from obspy import read
from glob import glob
import pandas as pd
import os
import numpy as np
import warnings
from datetime import datetime, timedelta, timezone
import scipy.io as io
import warnings
#warnings.filterwarnings("ignore")


def julian_to_datetime(julian_day_number):
    """
    Converts a Julian Day Number to a UTC datetime object.
    """
    # Reference date: November 17, 1858, 00:00:00 UTC (JDN 2400000.5)
    # We use a datetime object for 1858-11-17 00:00:00 UTC
    # and adjust for the 0.5 offset (noon) when calculating the timedelta.
    reference_date = datetime(2008, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
    reference_jdn = 1

    # Calculate the number of days difference
    days_difference = julian_day_number - reference_jdn

    # Add the timedelta to the reference date
    converted_datetime = reference_date + timedelta(days=days_difference)
    return converted_datetime

def make_catalog_from_root(root, to_file=False, filename=None):
    """
    Makes a DataFrame catalog of data from root directory.
    PARAMS:
        root (str): Path to data folder.
        to_file (bool): True means save catalog to file
        filename (str): filename of saved catalog

    RETURNS:
        catalog (DataFrame): catalog of dataset including
            Filename -> full file path
            Network -> seismic network
            Station -> seismic station
            Channel -> trace channel
            Start_Time -> start time from SAC file (datetime)
            End_Time -> end time from SAC file (datetime)

    Filename format should be [network].[station]..[channel].#.[%Y].[%j].[%H][%M][%S].SAC for proper sorting. The sorting assumes that filenames follow this pattern and the catalog/spool workflow may not work if this is not the case.
    There should be no duplicate files or data.
    """
    
    filenames = []
    starttimes = []
    endtimes = []
    networks = []
    stations = []
    channels = []

    files = glob(os.path.join(root, "**", "*.SAC"), recursive=True)
    for file in files:
                try:
                    data = read(file)
                    tr = data[0]  # SAC files contain one trace each
                    filenames.append(file)
                    starttimes.append(tr.stats.starttime.datetime)
                    endtimes.append(tr.stats.endtime.datetime)
                    networks.append(tr.stats.network)
                    stations.append(tr.stats.station)
                    channels.append(tr.stats.channel)
                    print(f"added file {file} to catalog")
                except Exception as e:
                    warnings.warn(f"Warning: Could not read {file} -> {e}")
                    
    if not filenames:
        raise FileNotFoundError(f"No SAC files found.")

    # Build the DataFrame
    catalog = pd.DataFrame({
        "Filename": filenames,
        "Network": networks,
        "Station": stations,
        "Channel": channels,
        "Start_Time": starttimes,
        "End_Time": endtimes,
    })
    catalog = catalog.sort_values(by='Filename')

    if to_file:
        if not filename:
            raise ValueError("Filename must be provided if to_file=True.")
        catalog.to_csv(filename, index=False)

    return catalog


def make_detection_catalog(detects, delta=0.005, length_seconds=20, cutw_seconds=2):
    """
    Build a DataFrame of detection start/end times from a MATLAB detection file.

    Parameters
    ----------
    detects : str
        Path to .mat file created by your detection script.
    delta : float
        Sample spacing in seconds (default 0.005 for 200 Hz).
    length_seconds : float
        Duration of each detection window to record in seconds.
    cutw_seconds : float
        Pre-window to include before detection (in seconds).
    """
    mat_data = io.loadmat(
        detects, squeeze_me=True, struct_as_record=False, spmatrix=False
    )

    ndata = mat_data["ndata"]            

    template_detects = mat_data["template_detects"]
    LOCS   = template_detects.LOCS         # list of arrays of sample indices
    PKS    = template_detects.PKS          # list of arrays of correlation values for picks
    days   = template_detects.day          # julian days
    hours  = template_detects.hour         # hours (0–23)

    cutw = int(round(cutw_seconds / delta))
    lenw = int(round(length_seconds / delta))

    start_times = []
    end_times   = []
    PKS_list    = []
    # loop over each template/hour entry
    for k, detTimes in enumerate(LOCS):
        if detTimes is None:
            continue

        day  = int(days[k])
        hour = int(hours[k])

        try:
            for l, det in enumerate(detTimes):
                det = int(det)                 # MATLAB findpeaks index (1-based)
                idx0 = max(det - cutw, 1)      # keep >=1 for MATLAB indexing
                idx1 = min(det + lenw, ndata)
    
                # convert to Python 0-based for delta multiplication
                start_sec = (idx0 - 1) * delta
                end_sec   = start_sec + length_seconds
    
                start_dt = julian_to_datetime(day) + timedelta(hours=hour, seconds=start_sec)
                end_dt   = start_dt + timedelta(seconds=length_seconds)
    
                start_times.append(np.datetime64(start_dt))
                end_times.append(np.datetime64(end_dt))
                PKS_list.append(PKS[k][l])
        except:
            idx0 = max(det - cutw, 1)      # keep >=1 for MATLAB indexing
            idx1 = min(det + lenw, ndata)
    
            # convert to Python 0-based for delta multiplication
            start_sec = (idx0 - 1) * delta
            end_sec   = start_sec + length_seconds

            start_dt = julian_to_datetime(day) + timedelta(hours=hour, seconds=start_sec)
            end_dt   = start_dt + timedelta(seconds=length_seconds)
    
            start_times.append(np.datetime64(start_dt))
            end_times.append(np.datetime64(end_dt))
            PKS_list.append(PKS[k])
    return pd.DataFrame({"Start_Times": start_times,
                         "End_Times":   end_times,
                        "PKS":PKS_list})
