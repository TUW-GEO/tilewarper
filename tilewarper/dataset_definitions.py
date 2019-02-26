from collections import OrderedDict
import pandas as pd
import os
import datetime
from geopathfinder.folder_naming import build_smarttree, SmartPath


class SRTM(object):
    """
    Represents Digital Elevation Model.
    """
    def __init__(self, root_dirpath):
        """
        Initialises folder and filenaming structure
        """
        # target directory path settings
        root_dirpath_ext = os.path.join(root_dirpath, 'DEM', 'SRTM')
        self.root_dirpath = root_dirpath_ext
        # target filenaming settings
        product_name = "CGIAR-CSI_v4.1_90M"
        self.fields_def = OrderedDict([('product_name', {'len': len(product_name)}),
                                      ('grid', {'len': 6, 'delim': True}),
                                      ('tile', {'len': 10})])
        self.fields_fixed = {'product_name': product_name}
        self.out_ndv = -9999
        # source data directory path
        self.ds_src_dirpath = r"/eodc/private/tuwgeo/datapool_raw/srtm/vfp/mosaics/"



class CorineLandCover(object):
    """

    Represents the Corine Land Cover map (see http://land.copernicus.eu/pan-european/corine-land-cover).

    Methods
    -------
    create_filename(self)
        Concatenates the internal filename dictionary to a string.

    """
    def __init__(self):
        """

        Initialises metadata attributes of the external dataset (e.g. the structure of the internal filename).

        """
        self.product_id = 'Land_Cover'
        self.filename_structure = OrderedDict()
        self.filename_structure["release"] = "2012"
        self.filename_structure["product_name"] = 'CLC'
        self.filename_structure["res"] = '100M'
        self.filename_structure["version"] = 'V18-5'
        self.dataset_dl_dir = '/eodc/private/tuwgeo/datapool_raw/external/corine_land_cover'
        self.filename = self.create_filename()

    def create_filename(self):
        """

        Concatenates the internal filename dictionary to a string.

        Returns
        -------
        str
            Internal filename of the external dataset.

        """
        return '_'.join(self.filename_structure.values())


class Landsat(object):
    """

    Represents the Corine Land Cover map (see http://land.copernicus.eu/pan-european/corine-land-cover).

    Methods
    -------
    create_filename(self)
        Concatenates the internal filename dictionary to a string.

    """
    def __init__(self):
        """

        Initialises metadata attributes of the external dataset (e.g. the structure of the internal filename).

        """
        self.product_id = 'Landsat_Overview'
        self.filename_structure = OrderedDict()
        self.filename_structure["product_name"] = "20170125_LS7_Overview"
        self.dataset_dl_dir = '/eodc/private/tuwgeo/users/cnavacch/data/external/Landsat/'
        self.filename = self.create_filename()

    def create_filename(self):
        """

        Concatenates the internal filename dictionary to a string.

        Returns
        -------
        str
            Internal filename of the external dataset.

        """
        return '_'.join(self.filename_structure.values())


# TODO: classes below need to be implemented and tested!
class MultiTemporal:
    def __init__(self):
        self.product_id = 'MuliTemp'
        self.pattern = []
        s_time = datetime.datetime(2014, 1, 1)
        e_time = datetime.datetime(2017, 1, 1)
        self.times = [timestamp.strftime('%Y%m%d') for timestamp in pd.date_range(s_time, e_time)]
        self.filename_structure = OrderedDict()
        self.filename_structure["timestamp"] = '{}'
        self.filename_structure["product_name"] = 'MULTTEMP'
        self.filename = self.create_filename()

    def create_filename(self):
        filenames = []
        for timestamp in self.times:
            self.pattern.append('{0}_.*'.format(timestamp))
            self.filename_structure["timestamp"] = timestamp
            filenames.append('_'.join(self.filename_structure.values()))
        return filenames




# dictionary which links the dataset definition classes with a dataset id/name
DATASET_MAP = {'SRTM': SRTM}


# TODO: also print folder and filenaming structure
def print_ds():
    """

    Prints all external dataset definitions.

    """
    for dataset_key in DATASET_MAP.keys():
        print(dataset_key)


def get_dataset_def(var_name):
    """

    Retrieve external dataset class from a given variable name.

    Parameters
    ----------
    var_name: str
        Internal id of the TileWarper.

    Returns
    -------
    object
        External dataset class.

    """
    if var_name not in DATASET_MAP.keys():
        print('The following external data sets are implemented: ')
        print_ds()
        raise KeyError('Choose one from above or add the unknown data set to dataset_definitions.py.')
    return DATASET_MAP[var_name]