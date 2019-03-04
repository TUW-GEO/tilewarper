import glob
import re
import os
from collections import OrderedDict
from tilewarper.dataset_definitions import get_dataset_def
from tilewarper.tilewarper import TileWarper
from pyraster.geotiff import read_tiff


class TileReader(object):
    """

    Reads a .tiff file (or a stack of .tiff files).

    Methods
    -------
    read(self, tiles, pattern=None, resample=False, in_filepath=None, out_filepath=None)
        Reads .tiff files from hard disk and returns them as an array or dictionary of arrays.
    get_filepath(self)
        Builds file path of the external dataset to tile level.
    query_data(filepath, pattern=None):
        Searches for .tif files in a given directory according to the provided pattern.
    """
    def __init__(self, root_dirpath, ds_id, gridname="EQUI7", sub_gridname='EU500M'):
        """

        Parameters
        ----------
        root_dirpath: str
            Filepath/directory root of the internal representation of the external dataset (or to resample to).
        ds_id: str
            Internal id of the external dataset (e.g. 'DEMSLOPE', 'CLC', ...).
        gridname: str (optional)
            Name/tag of the tiled projection system. (default 'EQUI7')
        sub_gridname: str (optional)
            Name/tag of the tiled projection + resolution. (default 'EU500M')

        """
        self.root_dirpath = root_dirpath
        self.gridname = gridname
        self.sub_gridname = sub_gridname
        self.ds_id = ds_id

    def read(self, tilenames, pattern=None, resample=False, root_dirpath=None, gdal_path=None):
        """

        Reads .tiff files from hard disk and returns them as an array or dictionary of arrays.

        Parameters
        ----------
        tilenames: str/list
            Name of the tiles to read (e.g. E048N012T6).
        pattern: str (optional)
            Regular expression to select a unique subset of files to read. (default None)
        resample: bool (optional)
            If True and the tile does not exist in the internal structure, SgrtImport() resamples the tile.
            (default False)
        in_filepath: str (optional)
            File path is necessary to specify a directory where the unchanged external data is located. If None,
            the file path provided to SgrtReader() is used. (default None)
        root_dirpath: str (optional)
            File path is necessary to specify a directory where the resampled external data should be written to.
            If None, the file path provided to SgrtReader() is used. (default None)

        """
        return_dict = False
        if type(tilenames) == str:
            tilenames = [tilenames]
        elif type(tilenames) == list:
            if len(tilenames) > 1:
                return_dict = True
        else:
            raise ValueError('Data type of tile is not supported.')
        data_dict = dict()
        ds_def = get_dataset_def(self.ds_id)
        for tilename in tilenames:
            ftilename = self.sub_gridname + '_' + tilename
            tile_dirpath = os.path.dirname(TileWarper.create_output_filepath(ds_def, self.gridname, self.sub_gridname,
                                                                             ftilename))
            tile_filepaths = self.query_data(tile_dirpath, pattern=pattern)
            if len(tile_filepaths) == 0:
                print('The region {0} has not been resampled yet.'.format(tilename))
                if resample:
                    print('Resampling {0} ...'.format(tilename))
                    if root_dirpath is None:
                        root_dirpath = self.root_dirpath
                    tilewarper = TileWarper(self.ds_id, root_dirpath, gridname=self.gridname,
                                            sub_gridname=self.sub_gridname, gdal_dirpath=gdal_path)
                    tilewarper.run(tilenames=tilename)
                    tile_filepaths = self.query_data(tile_dirpath, pattern=pattern)
                    if len(tile_filepaths) == 0:
                        print('Notice: for the tile {0} the pattern {1} matches no filename.'.format(tilename, pattern))
                else:
                    print('Use TileWarper with tile {0} for resampling first. Array is set to None.'.format(tilename))
                    tile_filepaths = None

            if tile_filepaths is not None:
                data_dict[tilename], _ = read_tiff(tile_filepaths[0])
            else:
                data_dict[tilename] = None

        if not return_dict:
            return list(data_dict.values())[0]
        else:
            return data_dict

    def get_filepaths(self, tilenames, pattern=None):

        if type(tilenames) == str:
            tilenames = [tilenames]

        data_dict = OrderedDict()
        ds_def = get_dataset_def(self.ds_id)(self.root_dirpath)
        for tilename in tilenames:
            ftilename = self.sub_gridname + '_' + tilename
            tile_dirpath = os.path.dirname(TileWarper.create_output_filepath(ds_def, self.gridname, self.sub_gridname,
                                                                             ftilename))
            filepaths = self.query_data(tile_dirpath, pattern=pattern)
            if len(filepaths) == 0:
                print('The region {0} has not been resampled yet.'.format(tilename))
                continue

            if filepaths is not None:
                data_dict[tilename] = filepaths[0]
            else:
                data_dict[tilename] = None

        return data_dict

    @staticmethod
    def query_data(filepath, pattern=None):
        """

        Searches for .tif files in a given directory according to the provided pattern.

        Returns
        -------
        filepath: str
           File path to .tif level (e.g. DEMSlope/SRTM-VFP/EQUI7_EU500M/E048N012T6).
        pattern: str (optional)
            Regular expression to select a subset of files to read. (default None)

        """
        file_ids = glob.glob(os.path.join(filepath, '*.tif'))
        if pattern is not None:
            regex = re.compile(pattern)
            file_ids = [x for x in file_ids if regex.match(x)]
        return file_ids


class CorineLandCoverReader(TileReader):
    """

    Reader for external dataset 'Corine Land Cover'.

    """
    def __init__(self, root_dirpath=None, gridname="EQUI7", sub_gridname='EU500M'):
        """

        Parameters
        ----------
        filepath: str
            Filepath/directory root of the internal representation of the external dataset (or to resample to).
        grid: str (optional)
            Name/tag of the tiled projection system. (default 'EQUI7')
        sub_grid: str (optional)
            Name/tag of the tiled projection + resolution. (default 'EU500M')

        """
        super(CorineLandCoverReader, self).__init__(root_dirpath, 'CLC', gridname=gridname, sub_gridname=sub_gridname)


class SRTMReader(TileReader):
    """

    Reader for external dataset 'DEM'.

    """
    def __init__(self, root_dirpath=None, gridname="EQUI7", sub_gridname='EU500M'):
        """

        Parameters
        ----------
        filepath: str
            Filepath/directory root of the internal representation of the external dataset (or to resample to).
        grid: str (optional)
            Name/tag of the tiled projection system. (default 'EQUI7')
        sub_grid: str (optional)
            Name/tag of the tiled projection + resolution. (default 'EU500M')

        """
        super(SRTMReader, self).__init__(root_dirpath, 'SRTM', gridname=gridname, sub_gridname=sub_gridname)


class SRTM_EGM2008Reader(TileReader):
    """

    Reader for external dataset 'DEM'.

    """
    def __init__(self, root_dirpath=None, gridname="EQUI7", sub_gridname='EU500M'):
        """

        Parameters
        ----------
        filepath: str
            Filepath/directory root of the internal representation of the external dataset (or to resample to).
        grid: str (optional)
            Name/tag of the tiled projection system. (default 'EQUI7')
        sub_grid: str (optional)
            Name/tag of the tiled projection + resolution. (default 'EU500M')

        """
        super(SRTM_EGM2008Reader, self).__init__(root_dirpath, 'SRTM_EGM2008', gridname=gridname,
                                                 sub_gridname=sub_gridname)