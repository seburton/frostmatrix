import numpy as np
import pandas as pd
import scipy.io as io
from datetime import datetime, timedelta, timezone

import catalog_management as cm
from Spool import spool
from readers import reader
from obspy import read


# get existing catalogs and make sure start and end times are np.datetime
data_catalog = pd.read_csv('data_catalog.csv')
data_catalog['Start_Time'] = pd.to_datetime(data_catalog['Start_Time']).values.astype('datetime64[ns]')
data_catalog['End_Time'] = pd.to_datetime(data_catalog['End_Time']).values.astype('datetime64[ns]')


# construct a catalog of detections
detections_catalog = cm.make_detection_catalog('/Users/sburton/Library/CloudStorage/OneDrive-Colostate/Work/bulk_data_extraction/erupt_detects_345.mat')
detections_catalog['Start_Times'] = pd.to_datetime(detections_catalog['Start_Times']).values.astype('datetime64[ns]')
detections_catalog['End_Times'] = pd.to_datetime(detections_catalog['End_Times']).values.astype('datetime64[ns]')

threshold = detections_catalog['PKS'].quantile(0.995)
top_tenp_detections = detections_catalog[detections_catalog['PKS'] >= threshold]

# initiate spool using data catalog and reader function (imported from readers)
sp = spool(data_catalog, reader, sr=200.0)

# extract top 0.005% of detections
for i, det in top_tenp_detections.iterrows():
    starttime = det['Start_Times']
    endtime = det['End_Times']

    sp.get_data(starttime, endtime, save_dir='Test/Templates/')