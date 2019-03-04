from collections import OrderedDict
import os


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


class EGM2008(object):
    """
    Represents Digital Elevation Model.
    """
    def __init__(self, root_dirpath):
        """
        Initialises folder and filenaming structure
        """
        # target directory path settings
        root_dirpath_ext = os.path.join(root_dirpath, 'GEOID', 'EGM2008')
        self.root_dirpath = root_dirpath_ext
        # target filenaming settings
        product_name = "EGM2008_2.5M"
        self.fields_def = OrderedDict([('product_name', {'len': len(product_name)}),
                                      ('grid', {'len': 6, 'delim': True}),
                                      ('tile', {'len': 10})])
        self.fields_fixed = {'product_name': product_name}
        self.out_ndv = -9999
        # source data directory path
        self.ds_src_dirpath = r"/eodc/private/tuwgeo/users/cnavacch/data/external/geoid/EGM2008/EGM2008.tif"


class SRTM_EGM2008(object):
    """
    Represents Digital Elevation Model.
    """

    def __init__(self, root_dirpath):
        """
        Initialises folder and filenaming structure
        """
        # target directory path settings
        root_dirpath_ext = os.path.join(root_dirpath, 'SRTM_EGM2008')
        self.root_dirpath = root_dirpath_ext
        # target filenaming settings
        product_name = "ellipsoid_heights_SRTM_EGM2008"
        self.fields_def = OrderedDict([('product_name', {'len': len(product_name)}),
                                       ('grid', {'len': 6, 'delim': True}),
                                       ('tile', {'len': 10})])
        self.fields_fixed = {'product_name': product_name}
        self.out_ndv = -9999
        # source data directory path
        self.ds_src_dirpath = None

# dictionary which links the dataset definition classes with a dataset id/name
DATASET_MAP = {'SRTM': SRTM, 'EGM2008': EGM2008, 'SRTM_EGM2008': SRTM_EGM2008}


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
