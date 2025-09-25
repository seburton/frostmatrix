import catalog_management as cm
from Spool import spool
from readers import reader
import numpy as np
import pandas as pd

#define root folder
root_folder = '/Volumes/EXTERNAL_1/Data_Tomo_Erebus/Data_345'

# construct catalog using catalog_management module
#catalog = cm.make_catalog_from_root(root_folder, to_file=True, filename='data_catalog_345.csv')
catalog = pd.read_csv('Test/data_catalog_345.csv')
catalog['Start_Time'] = pd.to_datetime(catalog['Start_Time'])
catalog['End_Time'] = pd.to_datetime(catalog['End_Time'])

# initiate spool using data catalog and reader function (imported from readers)
sp = spool(catalog, reader, sr=200.0)

# define start and end time of data to grab
starttime = np.datetime64("2008-12-10T12:00:00") + np.timedelta64(1306, 's')
endtime = np.datetime64("2008-12-10T12:00:00") + np.timedelta64(1326, 's')

merged_data_array = sp.get_data(starttime, endtime, save_dir='Test')