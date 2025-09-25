## Idea from work with Dr. Ge Jin at Colorado School of Mines Reservoir Characterization Project (RCP)
import pandas as pd
import numpy as np
import os
from obspy import read
import warnings
import scipy.io
import sys
import datetime as datetime


class spool:
    """
    Spool manages organized seismic dataset access with efficient metadata handling.

    A Pandas DataFrame (`df`) and a reader function are required 
    to construct a `spool` object. The `spool` handles caching SAC headers, 
    time-based queries, and waveform extraction.

    Parameters
    ----------
    df : pandas.DataFrame
        The catalog DataFrame describing the dataset. Must contain at least:
        - Filename  -> Full path to SAC file
        - Network   -> Seismic network code
        - Station   -> Seismic station code
        - Channel   -> Trace channel code
        - Start_Time -> Trace start time (datetime or np.datetime64)
        - End_Time   -> Trace end time (datetime or np.datetime64)
    reader : callable
        A function to load waveform data for a given SAC file.  
        Must support partial reads:
            `data = reader(filename, start_time, end_time)`
    sr : float
        Sample rate (Hz). Required for accurate waveform slicing.
    desired_sac_headers : list of str, optional
        A list of SAC header field names to extract into metadata when 
        fetching data. Defaults to:
            [
                "kstnm", "kevnm", "khole", "ko", "ka", "kt0", "kt1", "kt2", "kt3", "kt4",
                "kt5", "kt6", "kt7", "kt8", "kt9", "kf", "kuser0", "kuser1", "kuser2",
                "kcmpnm", "knetwk", "kdatrd", "kinst", "nzyear", "nzjday", "nzhour", "nzmin",
                "nzsec", "nzmsec", "sec", "nvhdr", "npts", "iftype", "idep", "iztype",
                "iqual", "isynth", "leven", "lpspol", "lovrok", "lcalda", "delta", "depmin",
                "depmax", "scale", "odelta", "b", "e", "o", "a", "t0", "t1", "t2", "t3",
                "t4", "t5", "t6", "t7", "t8", "t9", "resp0", "resp1", "resp2", "resp3",
                "resp4", "resp5", "resp6", "resp7", "resp8", "resp9", "stla", "stlo",
                "stel", "stdp", "evla", "evlo", "evel", "evdp", "user0", "user1", "user2",
                "user3", "user4", "user5", "user6", "user7", "user8", "user9", "dist",
                "az", "baz", "gcarc", "depmen", "cmpaz", "cmpinc"
            ]


    Attributes
    ----------
    _df : pandas.DataFrame or None
        Catalog of SAC metadata.
    _reader : callable or None
        User-provided function for loading waveform data.
    _sr : float or None
        Sample rate in Hz.
    _desired_sac_headers : list of str
        Subset of SAC headers to include in metadata output.
    cached_headers : bool
        True if SAC headers have been preloaded into the catalog.
    station_index : dict
        Cached lookup table for station/channel file ranges.
    file_duration : float or None
        Duration (seconds) of each SAC file. Assumes uniform file lengths.

    Notes
    -----
    Dataset Requirements:
        1. Constant sample rate across all files.
        2. Equal file lengths.
        3. Continuous coverage for each station over requested time ranges.
           The spool assumes SAC files are perfectly contiguous.
    """

    def __init__(self, df=None, reader=None, sr=None, desired_sac_headers=None):
        """
        Initialize a spool instance for managing seismic datasets.
    
        Parameters
        ----------
        df : pandas.DataFrame, optional
            Dataset catalog containing file paths, channels, stations, and time ranges.
            Default is None.
        reader : callable, optional
            Function to load waveform data for a file:
                `reader(filename, start_time, end_time)`
            Default is None.
        sr : float, optional
            Sampling rate (Hz). Required to compute the number of samples in the 
            requested time range. Default is None.
        desired_sac_headers : list of str, optional
            List of SAC header fields to include in the output metadata catalog.
            Defaults to [
                "kstnm", "kevnm", "khole", "ko", "ka", "kt0", "kt1", "kt2", "kt3", "kt4",
                "kt5", "kt6", "kt7", "kt8", "kt9", "kf", "kuser0", "kuser1", "kuser2",
                "kcmpnm", "knetwk", "kdatrd", "kinst", "nzyear", "nzjday", "nzhour", "nzmin",
                "nzsec", "nzmsec", "sec", "nvhdr", "npts", "iftype", "idep", "iztype",
                "iqual", "isynth", "leven", "lpspol", "lovrok", "lcalda", "delta", "depmin",
                "depmax", "scale", "odelta", "b", "e", "o", "a", "t0", "t1", "t2", "t3",
                "t4", "t5", "t6", "t7", "t8", "t9", "resp0", "resp1", "resp2", "resp3",
                "resp4", "resp5", "resp6", "resp7", "resp8", "resp9", "stla", "stlo",
                "stel", "stdp", "evla", "evlo", "evel", "evdp", "user0", "user1", "user2",
                "user3", "user4", "user5", "user6", "user7", "user8", "user9", "dist",
                "az", "baz", "gcarc", "depmen", "cmpaz", "cmpinc"
            ]
        """

        # Default SAC headers if not provided
        if desired_sac_headers is None:
            desired_sac_headers = [
                "kstnm", "kevnm", "khole", "ko", "ka", "kt0", "kt1", "kt2", "kt3", "kt4",
                "kt5", "kt6", "kt7", "kt8", "kt9", "kf", "kuser0", "kuser1", "kuser2",
                "kcmpnm", "knetwk", "kdatrd", "kinst", "nzyear", "nzjday", "nzhour", "nzmin",
                "nzsec", "nzmsec", "sec", "nvhdr", "npts", "iftype", "idep", "iztype",
                "iqual", "isynth", "leven", "lpspol", "lovrok", "lcalda", "delta", "depmin",
                "depmax", "scale", "odelta", "b", "e", "o", "a", "t0", "t1", "t2", "t3",
                "t4", "t5", "t6", "t7", "t8", "t9", "resp0", "resp1", "resp2", "resp3",
                "resp4", "resp5", "resp6", "resp7", "resp8", "resp9", "stla", "stlo",
                "stel", "stdp", "evla", "evlo", "evel", "evdp", "user0", "user1", "user2",
                "user3", "user4", "user5", "user6", "user7", "user8", "user9", "dist",
                "az", "baz", "gcarc", "depmen", "cmpaz", "cmpinc"
            ]

        # Store parameters
        self._df = df
        self._reader = reader
        self._sr = sr
        self._desired_sac_headers = desired_sac_headers

        # Warn on missing essentials
        if df is None:
            warnings.warn("No catalog DataFrame (`df`) provided to spool. Set it later with `set_database(df)`.")
        if reader is None:
            warnings.warn("No reader function provided to spool. Set it later with `set_reader(reader)`.")
        if sr is None:
            warnings.warn("No sample rate (`sr`) provided. This may cause issues when fetching waveform data.")

        # Initialize station indices only if df is provided
        self.station_index = {}
        self.file_duration = None
        if self._df is not None:
            self._precompute_station_indices()

        # Track SAC header cache state
        self.cached_headers = False


    def set_database(self, df):
        self._df = df
        self._precompute_station_indices()

    def _preload_sac_headers(self):
        """
        Preload and cache SAC header metadata for all files in the dataset.
    
        This method iterates over all SAC files listed in the catalog, reads their
        headers using ObsPy, and merges the extracted header fields into 
        `self._df`. All available SAC header fields (up to 95) are added as new
        columns in the catalog.
    
        If a file cannot be read or lacks SAC metadata, its header entry will be
        stored as an empty dictionary, and missing values will be represented as NaN.
    
        After successful completion, `self.cached_headers` is set to True.
    
        Notes
        -----
        - This function modifies `self._df` **in place**.
        - SAC headers are read in "head-only" mode for efficiency.
        - Uses `obspy.read(..., headonly=True)` internally.
    
        Raises
        ------
        ValueError
            If the dataset catalog (`self._df`) is missing or empty.
    
        Examples
        --------
        >>> sp._preload_sac_headers()
        >>> sp._df.columns[-5:]  # Inspect last few SAC header columns
        Index(['stla', 'stlo', 'evla', 'evlo', 'mag'], dtype='object')
        """
        header_dicts = []
        for file in self._df["Filename"]:
            try:
                st = read(file, headonly=True)[0]
                headers = dict(st.stats.sac) if hasattr(st.stats, "sac") else {}
            except:
                headers = {}
            header_dicts.append(headers)
        
        # Merge SAC headers into main dataframe
        headers_df = pd.DataFrame(header_dicts)
        self._df = pd.concat([self._df.reset_index(drop=True), headers_df], axis=1)
        self.cached_headers = True
        

    def _precompute_station_indices(self):
        """
        Precompute and cache start/end time indices for each station and channel.
    
        This method groups the dataset catalog (`self._df`) by station and channel,
        sorts files chronologically, and builds a lookup table (`self.station_index`)
        for fast access to time ranges. The resulting structure allows efficient
        binary search when retrieving waveform data.
    
        The function also computes `self.file_duration` using the first file's 
        start and end times, assuming all files have equal lengths and uniform 
        sample rates.
    
        Notes
        -----
        - The computed indices are stored in `self.station_index`, structured as:
          
          ```
          {
              station_code: {
                  channel_code: {
                      "Start_Times": np.ndarray,  # Sorted file start times
                      "End_Times": np.ndarray,    # Sorted file end times
                  },
                  ...
              },
              ...
          }
          ```
        - Assumes **contiguous, non-overlapping** SAC files for each channel.
        - If multiple channels exist for the same station, indices are computed separately.
    
        Raises
        ------
        ValueError
            If the dataset catalog (`self._df`) is missing or empty.
    
        Examples
        --------
        >>> sp._precompute_station_indices()
        >>> sp.station_index["ST01"]["HHZ"]["Start_Times"][:3]
        array(['2023-01-01T00:00:00', '2023-01-01T00:10:00', '2023-01-01T00:20:00'],
              dtype='datetime64[ns]')
        """
    
        self.station_index = {}
        self.file_duration = None
    
        for (sta, chan), group in self._df.groupby(["Station", "Channel"]):
            # Sort files by start time per channel
            chan_df = group.sort_values("Start_Time")
    
            start_times = chan_df["Start_Time"].to_numpy()
            end_times = chan_df["End_Time"].to_numpy()
    
            if sta not in self.station_index:
                self.station_index[sta] = {}
            self.station_index[sta][chan] = {
                "Start_Times": start_times,
                "End_Times": end_times
            }
    
            # Store file duration (constant across all files)
            if self.file_duration is None:
                self.file_duration = (
                    end_times[0] - start_times[0]
                ) / np.timedelta64(1, "s")

    def set_reader(self, reader):
        self._reader = reader

    def print_df(self):
        print(self._df)
        pass
            
    def check_time_in_range(self, station, channel, start_time, end_time):
        """
        Check if a requested time range is fully available for a given station/channel.
    
        Uses the precomputed `self.station_index` (built by `_precompute_station_indices`)
        to efficiently verify that both `start_time` and `end_time` fall within the
        available dataset for the specified station and channel.
    
        Parameters
        ----------
        station : str
            Station code (e.g., "ANMO").
        channel : str
            Channel code (e.g., "BHZ").
        start_time : datetime or np.datetime64
            Requested start time for the waveform data.
        end_time : datetime or np.datetime64
            Requested end time for the waveform data.
    
        Returns
        -------
        bool
            True if the entire time range is contained within the dataset;
            False otherwise.
    
        Raises
        ------
        KeyError
            If the station or channel is not present in the catalog.
        ValueError
            If `start_time` is after `end_time`.
        RuntimeError
            If `_precompute_station_indices` has not been called and
            `self.station_index` is undefined.
    
        Notes
        -----
        - Uses **binary search** via `np.searchsorted` for O(log N) lookups.
        - Assumes SAC files for each station/channel are **contiguous** and sorted.
        - If any portion of the requested range falls outside the dataset,
          the function returns False.
    
        Examples
        --------
        >>> sp.check_time_in_range("ANMO", "BHZ",
        ...                       start_time=np.datetime64("2023-01-01T00:10:00"),
        ...                       end_time=np.datetime64("2023-01-01T00:20:00"))
        True
    
        >>> sp.check_time_in_range("ANMO", "BHZ",
        ...                       start_time=np.datetime64("2025-01-01T00:00:00"),
        ...                       end_time=np.datetime64("2025-01-01T00:10:00"))
        False
        """
        # Check indices are available
        if not self.station_index:
            raise RuntimeError("Station indices are not initialized. "
                               "Call `_precompute_station_indices()` first.")
    
        # Validate station/channel existence
        if station not in self.station_index:
            raise KeyError(f"Station '{station}' not found in catalog.")
        if channel not in self.station_index[station]:
            raise KeyError(f"Channel '{channel}' not found for station '{station}'.")
    
        # Validate time ordering
        if start_time >= end_time:
            raise ValueError("`start_time` must be strictly earlier than `end_time`.")
    
        # Get cached start and end times for this station/channel
        sta_chan_index = self.station_index[station][channel]
        start_times = sta_chan_index["Start_Times"]
        end_times = sta_chan_index["End_Times"]
    
        # If the entire requested range is outside the available data, return False immediately
        if start_time < start_times[0] or end_time > end_times[-1]:
            return False
    
        # Find the first file that ends *after* start_time
        first_idx = np.searchsorted(end_times, start_time, side="right")
        # Find the last file that starts *before* end_time
        last_idx = np.searchsorted(start_times, end_time, side="left") - 1
    
        # Check if we have full coverage (contiguous segments)
        if first_idx <= last_idx:
            return True
        else:
            return False

    def get_data(
        self,
        start_time,
        end_time,
        save_dir=None,
        load_headers=True
    ):
        """
        Fetch waveform data and associated metadata for a given time range.
    
        Efficiently selects only the required SAC files, fetches their waveform
        data via the provided reader function, and optionally saves both the data
        and a metadata catalog to disk.
    
        Parameters
        ----------
        start_time : datetime or np.datetime64
            Start of the requested time window.
        end_time : datetime or np.datetime64
            End of the requested time window.
        save_dir : str, optional
            Directory to save waveform `.mat` file.
            If None, results are not saved.
        load_headers : bool, default=True
            Whether to preload all SAC headers if not already cached.
    
        Returns
        -------
        data_array : np.ndarray
            2D array of waveform data, shape = `(num_traces, num_samples)`.
            Each row corresponds to one station-channel trace.
        metadata_catalog : pandas.DataFrame
            Metadata catalog for the returned traces. Includes:
                - SAC headers from `self._desired_sac_headers`
    
        Raises
        ------
        ValueError
            If the requested time range is outside the dataset bounds.
        RuntimeError
            If the reader function is missing.
        """
        # -----------------------------
        # Define Constants
        # -----------------------------
        LL_lat = -77.5274
        LL_lon = 167.1645
        LL_elev = 3572  # meters
        daystr = str(start_time.to_pydatetime().toordinal() - datetime.datetime(2008,1,1).toordinal() +1)
        # -----------------------------
        # Ensure reader is defined
        # -----------------------------
        if self._reader is None:
            raise RuntimeError("No reader function set. Use `set_reader(reader)` first.")
    
        # -----------------------------
        # Preload SAC headers if requested
        # -----------------------------
        if load_headers and not self.cached_headers:
            self._preload_sac_headers()
            
        # -----------------------------
        # Select files overlapping the requested window
        # -----------------------------
        ind = np.where(
            (self._df["Start_Time"] < end_time) &
            (self._df["End_Time"] > start_time)
        )[0]
        df_subset = self._df.iloc[ind]
    
        if df_subset.empty:
            warnings.warn(f"No data found between {start_time} and {end_time}.")
            return np.array([]), pd.DataFrame()

        # ------------------------------
        # Preallocate data array
        # ------------------------------
        len_trace_seconds = (end_time - start_time) / np.timedelta64(1, "s")
        len_trace_series = int(len_trace_seconds * self._sr)+1
        num_traces = len(df_subset[["Station", "Channel"]].drop_duplicates())
        data_array = np.zeros((num_traces, len_trace_series), dtype=np.float32)
        metadata_rows = []
        trace_idx = 0
        ranges = np.zeros(num_traces, dtype=np.float64)
        # -----------------------------
        # Group by station & channel
        # -----------------------------
        grouped = df_subset.groupby(["Station", "Channel"])

        for (_, chan), chan_df in grouped:
            # Only process certain channels
            if chan not in ['ELZ','ELN','ELE']:
                continue

            # Sort files by start time
            chan_df = chan_df.sort_values("Start_Time").reset_index(drop=True)
            chan_data_list = []
            
            # -----------------------------
            # Loop over SAC files for this channel
            # -----------------------------

            for _, row in chan_df.iterrows():
                file = row["Filename"]
                try:
                    # Read waveform data for this file segment
                    dat = self._reader(file, start_time, end_time)
                    
                    # Skip dead reads
                    if dat.size == 0:
                        continue
                    if np.all(dat == 0):
                        continue
                    
                    chan_data_list.append(dat)
    
                    # Collect SAC headers into metadata
                    sac_headers = {h: row[h] if h in row else -12345 for h in self._desired_sac_headers}
                    
                    # Compute range (distance from lava lake LL)
                    if "stla" in sac_headers and "stlo" in sac_headers and sac_headers["stla"] != -12345 and sac_headers["stlo"] != -12345:
                        # 3D Euclidean distance (lat, lon, elev)
                        R = 6371.0  # Earth radius in km
                        lat1 = np.radians(LL_lat)
                        lon1 = np.radians(LL_lon)
                        lat2 = np.radians(sac_headers["stla"])
                        lon2 = np.radians(sac_headers["stlo"])
                        elev1 = LL_elev / 1000.0  # convert meters to km
                        elev2 = sac_headers.get("stel", 0) / 1000.0  # convert meters to km

                        # Convert spherical to Cartesian coordinates
                        x1 = (R + elev1) * np.cos(lat1) * np.cos(lon1)
                        y1 = (R + elev1) * np.cos(lat1) * np.sin(lon1)
                        z1 = (R + elev1) * np.sin(lat1)

                        x2 = (R + elev2) * np.cos(lat2) * np.cos(lon2)
                        y2 = (R + elev2) * np.cos(lat2) * np.sin(lon2)
                        z2 = (R + elev2) * np.sin(lat2)

                        distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
                        ranges[trace_idx] = distance
                    else:
                        ranges[trace_idx] = -1.0  # Unknown range
                        warnings.warn(f"Unknown range for {file}")

                except Exception as e:
                    warnings.warn(f"Failed to read {file}: {e}")
    
            # -----------------------------
            # Merge all segments for this channel
            # -----------------------------
            if chan_data_list:
                chan_data = np.concatenate(chan_data_list)
                metadata_rows.append(sac_headers)
            else:
                chan_data = np.zeros(len_trace_series, dtype=np.float32)
            
            data_array[trace_idx, :len(chan_data)] = chan_data[:len_trace_series]
            trace_idx += 1
    
        # -----------------------------
        # Build final metadata DataFrame
        # -----------------------------
        metadata_catalog = pd.DataFrame(metadata_rows)

        # -----------------------------
        # Sort data and metadata by range and drop zero traces
        # -----------------------------
        sort_indices = np.argsort(ranges)
        data_array = data_array[sort_indices, :]
        metadata_catalog = metadata_catalog.iloc[sort_indices].reset_index(drop=True)
        ranges = ranges[sort_indices]


        # Drop zero traces
        non_zero_indices = np.where(ranges > 0)[0]
        data_array = data_array[non_zero_indices, :]
        ranges = ranges[non_zero_indices]
        metadata_catalog = metadata_catalog.iloc[non_zero_indices].reset_index(drop=True)
        
        # -----------------------------
        # Optionally save outputs
        # -----------------------------
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            scipy.io.savemat(os.path.join(save_dir, f"template_{daystr}_{start_time.hour}0000_{start_time.minute*60 + start_time.second}_{end_time.minute*60 + end_time.second}.mat"), 
{'template':data_array.T,
                                                                               'template_nchans':num_traces,
                                                                                'template_range':ranges,
                                                                                'template_sachdr':metadata_catalog.to_dict(orient='records'),
                                                                                'x1': (
                                                                                    start_time.minute * 60
                                                                                    + start_time.second
                                                                                ) * self._sr,
                                                                               'x2': (
                                                                                   end_time.minute * 60
                                                                                   + end_time.second
                                                                                ) * self._sr                                                                           
                                                                               })
    
        return data_array

