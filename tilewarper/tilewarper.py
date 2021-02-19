# Copyright (c) 2014, Vienna University of Technology (TU Wien), Department
# of Geodesy and Geoinformation (GEO).
# All rights reserved.
#
# All information contained herein is, and remains the property of Vienna
# University of Technology (TU Wien), Department of Geodesy and Geoinformation
# (GEO). The intellectual and technical concepts contained herein are
# proprietary to Vienna University of Technology (TU Wien), Department of
# Geodesy and Geoinformation (GEO). Dissemination of this information or
# reproduction of this material is forbidden unless prior written permission
# is obtained from Vienna University of Technology (TU Wien), Department of
# Geodesy and Geoinformation (GEO).

'''
Created in February, 2019

Resamples a specific geospatial dataset (one or multiple GeoTiff files) to an internal structure (provided in
dataset_definitions.py) for a given grid and sub-grid

@author: Claudio Navacchi claudio.navacchi@geo.tuwien.ac.at
'''

from tilewarper.dataset_definitions import get_dataset_def
import re
from osgeo import osr
from pyraster import gdalport
import os
import glob
from collections import OrderedDict
import math
import gdal
from geopathfinder.folder_naming import SmartPath
from geopathfinder.file_naming import SmartFilename
# tiled projection systems imports:
from equi7grid.equi7grid import Equi7Grid
#from latlongrid.latlongrid import LatLonGrid

# TODO check multitemporal datasets and files with multiple bands
# TODO implement query for input data sets (using geopathfinder)

PROJ_TILE_SYS = {'EQUI7': Equi7Grid}  # fill it up with other projections


class TileWarper(object):

    """

    Warps and resamples pre-defined data sets (see dataset_definitions.py).

    Methods
    -------
    start(self, tiles=None)
        Starts the resampling and renaming process.
    _find_overlapping_tiles(self, ftile)
        Finds all source .tiff file/s, which overlap with the given tile.
    _create_folder_struct(self, ftile)
        Creates folder structure for the external dataset for a given tile id.
    query_data_all(filepath)
        Searches for .tif files in all sub-folders.
    resample(self, infiles, outfile, res, resampling_type="near", ndv=None, dst_srs=None, extent=None, gdal_dirpath=None)
        Resamples .tiff files.

    """

    def __init__(self, ds_name, root_dirpath, in_filepath=None, gridname='EQUI7', sub_gridname='EU500M',
                 gdal_dirpath=None):
        """

        Parameters
        ----------
        ds_name: str
            Internal name of the external dataset (e.g. 'SRTM', 'CLC', ...).
        in_filepath: str
            File path to .tif files or tile folders of the external dataset.
        root_dirpath: str
            Output file path where the resampled .tiff files should be written to.
        gridname: str (optional)
            Name/tag of the tiled projection system. (default 'EQUI7')
        sub_gridname: str (optional)
            Name/tag of the tiled projection + resolution. (default 'EU500M')

        """

        self.ds_name = ds_name.upper()
        self.ds_def = get_dataset_def(self.ds_name)(root_dirpath)
        if in_filepath is None:
            self.in_filepath = self.ds_def.ds_src_dirpath
        else:
            self.in_filepath = in_filepath
        self.gdal_dirpath = gdal_dirpath

        # get all files
        self.in_filepaths = self.query_data_all(self.in_filepath)

        # general TiledProjectionSystem settings
        TPS = PROJ_TILE_SYS[gridname]
        self.res = self.decode_res(gridname, sub_gridname)
        self.grid = TPS(self.res)
        self.sub_grid = getattr(self.grid, sub_gridname[:2].upper())

    @staticmethod
    def decode_res(gridname, sub_gridname):
        TPS = PROJ_TILE_SYS[gridname]
        res_string_options = [sub_gridname, sub_gridname[2:-1], sub_gridname[2:]]
        res = None
        for res_string_option in res_string_options:
            try:
                res = TPS.decode_sampling(res_string_option)
                break
            except:
                pass

        if res is None:
            raise Exception('Cannot find spatial resolution in variable "sub_gridname".')
        else:
            return res

    def run(self, tilenames=None):
        """

        Finds common tiles (overlap of the tiled projection and the covered region of the external dataset) and
        resamples the source .tif file/s.

        Parameters
        ----------
        tilenames: list (optional)
            Tile names of the target tiled projection system.
            If None, all tile names of the sub grid are used. (default None)

        """
        sub_gridname_parts = self.sub_grid.name.split('_')
        if tilenames is None:
            ftilenames = self.sub_grid.search_tiles_in_geometry(self.sub_grid.polygon_geog, coverland=False)
        elif type(tilenames) != list:
            ftilename = sub_gridname_parts[-1] + '_' + tilenames
            ftilenames = [ftilename]
        else:
            ftilenames = [tilename if len(tilename) > 11 else sub_gridname_parts[-1] + '_' + tilename
                          for tilename in tilenames]

        for ftilename in ftilenames:
            overlapping_filepaths = self._find_overlapping_filepaths(ftilename)
            if len(overlapping_filepaths) == 0:
                print('No overlap for tile {0}'.format(ftilename))
                continue

            out_filepath = self.create_output_filepath(self.ds_def, sub_gridname_parts[0], sub_gridname_parts[1],
                                                       ftilename)
            tile = self.sub_grid.tilesys.create_tile(name=ftilename)
            tile_extent = tile.limits_m()
            tile_sref = self.sub_grid.core.projection.wkt
            if hasattr(self.ds_def, 'pattern'):
                pattern = self.ds_def.pattern
                pass
            else:
                self.warp(overlapping_filepaths, out_filepath, self.res, dst_srs=tile_sref, extent=tile_extent,
                          gdal_dirpath=self.gdal_dirpath)

    def _find_overlapping_filepaths(self, ftilename):
        """

        Finds all source .tiff file/s, which overlap with the given tile.

        Parameters
        ----------
        ftilename: str
            Tile id consisting of sub grid and tile name (e.g. EU500M_E048N012T6).

        Returns
        -------
        overlapping_filepaths: list
            List containing overlapping tiles of the external dataset.

        """
        overlapping_filepaths = []
        for filepath in self.in_filepaths:
            # get and project extent of external tiff
            ds = gdal.Open(filepath, gdal.GA_ReadOnly)
            ds_img = gdalport.GdalImage(ds, filepath)
            ds_sr = osr.SpatialReference()
            ds_sr.ImportFromWkt(ds_img.projection())
            overlapping_tiles = self.grid.search_tiles_in_roi(extent=ds_img.get_extent(), osr_spref=ds_sr,
                                                              subgrid_ids=[self.sub_grid.core.tag], coverland=False)

            if ftilename in overlapping_tiles:
                overlapping_filepaths.append(filepath)

            # free memory
            ds = None
            ds_img = None

        return overlapping_filepaths

    @staticmethod
    def create_output_filepath(ds_def, gridname, sub_gridname, ftilename):
        """

        Creates folder structure for the external dataset for a given tile id.

        Parameters
        ----------
        ftile: str
            Tile id consisting of sub grid and tile name (e.g. EU500M_E048N012T6).

        Returns
        -------
        tile_path: str
            File path to the tile (e.g. E048N012T6) folder. It has the following structure:
            [product id]/[product name]/sub grid]/[tile name] (e.g. DEMSlope/SRTM-VFP/EQUI7_EU500M/E048N012T6).

        """
        ftilename_parts = ftilename.split('_')
        hierarchy = ['root', 'grid', 'tile']
        levels = {'root': ds_def.root_dirpath, 'grid': gridname + "_" + sub_gridname, 'tile': ftilename_parts[1]}
        smart_path = SmartPath(levels, hierarchy, make_dir=True)
        fields = ds_def.fields_fixed
        fields['grid'] = ftilename_parts[0]
        fields['tile'] = ftilename_parts[1]
        smart_filename = SmartFilename(fields, ds_def.fields_def, ext='.tif')
        out_filepath = os.path.join(smart_path.get_dir(), str(smart_filename))

        return out_filepath

    @staticmethod
    def query_data_all(filepath):
        """

        Searches for .tif files in all sub-folders. If the file path is actually a file,
        it is directly returned as a list.

        Parameters
        ----------
        filepath: str
            Root file path of the .tif sub-folder search.

        Returns
        -------
        filepaths: list
            List of .tif file names.

        """
        filepaths = []
        if not os.path.isfile(filepath):
            for dir_path, _, _ in os.walk(filepath):
                filepaths.extend(glob.glob(os.path.join(dir_path, '*.tif')))
        else:
            filepaths = [filepath]
        return filepaths

    def warp(self, in_filepaths, out_filepath, res, resampling_type="cubic", ndv=None, dst_srs=None, extent=None,
             gdal_dirpath=None, overwrite=True):
        """

        Resamples .tiff files. Includes merging (if multiple files are given) and warping.

        Parameters
        ----------
        in_filepaths: list
            Contains full filepath of input file/s.
        out_filepath: str
            Full filename of output file.
        res: int
            Resolution of output file.
        resampling_type: str (optional)
            Method, which should be used for interpolation (e.g. "near", "bilinear", "cubic", ...). (default "near")
        ndv: int (optional)
            No data value of output file. (default None)
        dst_srs: str
            Representation of spatial reference system (e.g. EPSG code ("EPSG:n"), WKT string, ...). If None,
            the source spatial reference system is used. (default None)
        extent: list (optional)
            Contains lower left and upper right coordinates of the region of interest [xmin, ymin, xmax, ymax].
            (default None)
        gdal_dirpath: str (optional)
            Path to the GDAL executables. (default None)

        """
        if os.path.exists(out_filepath) and overwrite:
            os.remove(out_filepath)

        # compression settings
        comp = ["COMPRESS=LZW"]
        tilesize = 512
        if tilesize:
            tilesize = int(tilesize)
            # make sure the tilesize is exponent of 2
            tilesize = 2 ** int(round(math.log(tilesize, 2)))
            comp.append("TILED=YES")
            comp.append("BLOCKXSIZE={:d}".format(tilesize))
            comp.append("BLOCKYSIZE={:d}".format(tilesize))

        if type(in_filepaths) != list:
            in_filepaths = [in_filepaths]
        # get no data values and colormap
        out_ndv = ndv
        gdal_ds = gdal.Open(in_filepaths[0], gdal.GA_ReadOnly)
        gdal_img = gdalport.GdalImage(gdal_ds, in_filepaths[0])
        in_ndv = gdal_img.get_band_nodata()
        in_ct = gdal_img.colormap()
        if hasattr(self.ds_def, 'out_ndv'):
            out_ndv = self.ds_def.out_ndv

        # merge files
        merge_success = False
        warp_filepath = in_filepaths[0]
        tmp_filepath = ''
        if len(in_filepaths) > 1:  # call gdal merge
            tmp_filename = os.path.splitext(os.path.basename(warp_filepath))[0] + '_temp.tif'
            tmp_filepath = os.path.join(os.path.dirname(out_filepath), tmp_filename)
            warp_filepath = tmp_filepath
            gdal_merge_opt = OrderedDict()
            gdal_merge_opt['-o'] = tmp_filepath
            if in_ct is not None:
                gdal_merge_opt['-pct'] = ''
            if in_ndv is not None:
                gdal_merge_opt['-n'] = str(in_ndv)
                gdal_merge_opt['-a_nodata'] = str(in_ndv)
                gdal_merge_opt['-init'] = str(in_ndv)
            else:
                in_ndv = 0
                if out_ndv is not None:
                    in_ndv = out_ndv
                else:
                    print('Notice: no data value has been set to 0.')
                gdal_merge_opt['-init'] = str(in_ndv)
            gdal_merge_opt['-co'] = comp

            print('Merging {0} -> {1} ...'.format(', '.join([os.path.basename(in_filepath)
                                                             for in_filepath in in_filepaths]),
                                                  os.path.basename(tmp_filepath)))
            merge_success, _ = gdalport.call_gdal_util('gdal_merge.py', src_files=in_filepaths, options=gdal_merge_opt,
                                                       gdal_path=gdal_dirpath)
            if merge_success:
                print(' - Done')
            else:
                print(' - Failed')

        # warping file
        gdal_warp_opt = OrderedDict()
        gdal_warp_opt['-r'] = resampling_type
        gdal_warp_opt['-tr'] = str(res) + ' ' + str(res)
        gdal_warp_opt['-srcnodata'] = str(in_ndv)
        if out_ndv is None:
            out_ndv = in_ndv
        gdal_warp_opt['-dstnodata'] = str(out_ndv)
        if dst_srs is not None:
            gdal_warp_opt['-t_srs'] = dst_srs
        if extent is not None:
            gdal_warp_opt['-te'] = ' '.join([str(coord) for coord in extent])
        gdal_warp_opt['-co'] = comp

        print('Warping {0} -> {1} ...'.format(os.path.basename(warp_filepath), os.path.basename(out_filepath)))
        warp_success, _ = gdalport.call_gdal_util('gdalwarp', src_files=warp_filepath, dst_file=out_filepath,
                                                  options=gdal_warp_opt, gdal_path=gdal_dirpath)
        if warp_success:
            print(' - Done')
        else:
            print(' - Failed')

        # Finally, delete tmp_merged file from disk
        if merge_success:
            os.remove(tmp_filepath)



