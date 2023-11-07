import numpy as np
import os
import pandas as pd

from astropy.io import fits
from copy import copy, deepcopy
from scipy.io import readsav

cobe_idl_linfrq = np.array(
    [
        115.3,
        230.5,
        345.8,
        424.8,
        461.0,
        492.2,
        556.9,
        576.3,
        691.5,
        808.1,
        1113,
        1460,
        2226,
        1901,
        2060,
        2311,
        2459,
        2589,
        921.8,
    ]
)
# cobe_idl_transitions = np.array(['CO 1', 'CO 2', 'CO 3', 'CO 4', 'C 3', 'CO 5',
#                                  'CO 6', 'CO 7 + C 1', 'C+ 1', 'O 1', 'CO 8'])
cobe_idl_transitions = np.array(
    [
        "CO 1",
        "CO 2",
        "CO 3",
        "CO 4",
        "C 1",
        "CO 5",
        "CO 6",
        "CO 7 + C 2",
        "C+ 1",
        "CO 8",
    ]
)
# cobe_idl_indeces = np.array([0, 1, 2, 4, 5, 7, 8, 9, 13, 14, 18])
cobe_idl_indeces = np.array([0, 1, 2, 4, 5, 7, 8, 9, 13, 18])
missions_2d = ["COBE-FIRAS", "COBE-DIRBE", "Planck"]


class Observation(object):
    """
    This is an object to load individual observations. This is merely for
    the convenience of opening all of the observations in a consistent manner.
    There is an optional argument when initialising to set a base directory,
    which makes it easier to load multiple models in succession.
    """

    def __init__(self, base_dir="", regridded_dir="/regridded/temp/"):
        """
        This initialises the object along with the base directory.
        The owned objects of `base_dir` and `files` are created.
        `files` can be modified again when loading a model, but for now it
        has the default filenames created with `kosmatau3d`.

        :param base_dir: the base directory to use when loading observations.
            Default: `''`.
        :param regridded_dir: the directory to use to load regridded
            observations. Default: `'/regridded/temp/'`.


        """

        if len(base_dir):
            if not base_dir[-1] == "/":
                base_dir += "/"
        self.base_dir = base_dir
        if len(regridded_dir):
            if not base_dir[0] == "/":
                base_dir = "/" + base_dir
            if not base_dir[-1] == "/":
                base_dir += "/"
        self.regridded_dir = regridded_dir
        # self.files = {'intensity': 'synthetic_intensity',
        #               'optical_depth': 'synthetic_optical_depth',
        #               'dust_absorption': 'dust_absorption',
        #               'dust_emissivity': 'dust_emissivity',
        #               'species_absorption': 'species_absorption',
        #               'species_emissivity': 'species_emissivity',
        #               'density': 'voxel_density',
        #               'ensemble_dispersion': 'voxel_ensemble_dispersion',
        #               'ensemble_mass': 'voxel_ensemble_mass',
        #               'fuv_absorption': 'voxel_FUVabsorption',
        #               'fuv': 'voxel_fuv',
        #               'position': 'voxel_position',
        #               'velocity': 'voxel_velocity',
        #               'los_count': 'sightlines',
        #               'log': 'log', }

        return

    def reset_attributes(self):
        """
        Reinitialise instance attributes.
        """

        self.obs = []
        self.obs_error = []
        self.obs_error_complete = []
        self.obs_error_conf = []
        self.obs_header = []
        self.obs_data = []
        self.obs_error_data = []
        self.obs_error_complete_data = []
        self.obs_error_conf_data = []
        self.obs_lon = []
        self.obs_lat = []
        self.obs_vel = []
        self.obs_iid = []
        self.obs_i_iid = []
        self.obs_wavelength = []
        self.obs_frequency = []
        self.obs_spectra = None
        self.obs_spectra_resampled = None
        self.obs_spectra_reduced = None

        return

    def load_survey(self, directory="", survey=None, lat=0):
        """
        Load a survey into memory.
        """

        if survey is None:
            print("ERROR: Please specify an observational survey")
            return

        self.survey = survey

        if os.path.exists(self.base_dir + directory + self.survey):
            self.directory = directory
        else:
            print(
                f"ERROR: Survey {self.base_dir + directory + survey} does "
                + "not exist. Ensure that base_dir was initialised "
                + "correctly and that you have the correct survey"
            )
            return

        if os.path.exists(
            self.base_dir + self.directory + self.survey + self.regridded_dir
        ):
            full_path = (
                self.base_dir + self.directory + self.survey + self.regridded_dir
            )
            spectra_path = self.base_dir + self.directory + self.survey + "/spectra/"
        else:
            print(
                "ERROR: Either the survey has not been regridded or the "
                + f"directory differs from {self.regridded_dir}."
            )
            return

        self.files = list(
            f.name
            for f in os.scandir(full_path)
            if (f.is_file() and not ("_error" in f.name or ".fuse_hidden" in f.name))
        )

        self.reset_attributes()

        for f in self.files:
            # Open IDL files if necessary
            if ".idl" in f:
                self.obs.append(readsav(full_path + f))
                self.obs_header.append([None])
                self.obs_iid.append(deepcopy(cobe_idl_transitions))
                self.obs_i_iid.append(deepcopy(cobe_idl_indeces))
                self.obs_frequency.append(deepcopy(cobe_idl_linfrq) * 1e9)
                self.obs_wavelength.append(299792458 / self.obs_frequency[-1])
                self.obs_data.append(
                    self.obs[-1]["amplitude"]
                    * (2.9979**3)
                    / ((self.obs_frequency[-1] / 1e9) ** 3 * 2 * 1.38)
                    * 10**8
                )
                self.obs_error.append(None)
                self.obs_error_data.append(
                    self.obs[-1]["sigma"]
                    * (2.9979**3)
                    / ((self.obs_frequency[-1] / 1e9) ** 3 * 2 * 1.38)
                    * 10**8
                )
                error_conf = np.zeros_like(self.obs_data[-1])
                for _ in range(int(np.floor(self.obs_data[-1].shape[0] / 2))):
                    idx = np.array([_, -1 - _])
                    err = self.obs_data[-1][idx, :]
                    error_conf[idx, :] = np.std(
                        err - np.mean(err, axis=0), axis=0, ddof=1
                    ) / np.sqrt(err.shape[0])
                self.obs_error_conf.append(None)
                self.obs_error_conf_data.append(deepcopy(error_conf))
                self.obs_lon.append(deepcopy(self.obs[-1]["long"]))
                self.obs_lat.append(np.array([0]))
                self.obs_vel.append([None])
            # By default open data stored in a FITS file
            else:  # if '.fits' in f:
                self.obs.append(fits.open(full_path + f))
                self.obs_data.append(self.obs[-1][0].data)
                self.obs_error.append(
                    fits.open(full_path + f.replace(".fits", "_error.fits"))
                )
                self.obs_error_data.append(self.obs_error[-1][0].data)
                if os.path.exists(
                    full_path + f.replace(".fits", "_complete_error.fits")
                ):
                    self.obs_error_complete.append(
                        fits.open(
                            full_path + f.replace(".fits", "_complete_error.fits")
                        )
                    )
                    self.obs_error_complete_data.append(
                        self.obs_error_complete[-1][0].data
                    )
                else:
                    self.obs_error_complete.append([])
                    self.obs_error_complete_data.append([])
                if os.path.exists(full_path + f.replace(".fits", "_error_conf.fits")):
                    self.obs_error_conf.append(
                        fits.open(full_path + f.replace(".fits", "_error_conf.fits"))
                    )
                    self.obs_error_conf_data.append(self.obs_error_conf[-1][0].data)
                else:
                    self.obs_error_conf.append(None)
                    self.obs_error_conf_data.append(0)
                self.obs_header.append(self.obs[-1][0].header)
                self.obs_lon.append(
                    np.linspace(
                        self.obs_header[-1]["CRVAL1"]
                        - self.obs_header[-1]["CDELT1"]
                        * (self.obs_header[-1]["CRPIX1"] - 1),
                        self.obs_header[-1]["CRVAL1"]
                        + self.obs_header[-1]["CDELT1"]
                        * (
                            self.obs_header[-1]["NAXIS1"]
                            - self.obs_header[-1]["CRPIX1"]
                        ),
                        num=self.obs_header[-1]["NAXIS1"],
                    )
                )
                self.obs_lat.append(
                    np.linspace(
                        self.obs_header[-1]["CRVAL2"]
                        - self.obs_header[-1]["CDELT2"]
                        * (self.obs_header[-1]["CRPIX2"] - 1),
                        self.obs_header[-1]["CRVAL2"]
                        + self.obs_header[-1]["CDELT2"]
                        * (
                            self.obs_header[-1]["NAXIS2"]
                            - self.obs_header[-1]["CRPIX2"]
                        ),
                        num=self.obs_header[-1]["NAXIS2"],
                    )
                )
                if self.obs_header[-1]["NAXIS"] == 3 and not self.survey in missions_2d:
                    self.obs_vel.append(
                        np.linspace(
                            (
                                self.obs_header[-1]["CRVAL3"]
                                - self.obs_header[-1]["CDELT3"]
                                * (self.obs_header[-1]["CRPIX3"] - 1)
                            ),
                            (
                                self.obs_header[-1]["CRVAL3"]
                                + self.obs_header[-1]["CDELT3"]
                                * (
                                    self.obs_header[-1]["NAXIS3"]
                                    - self.obs_header[-1]["CRPIX3"]
                                )
                            ),
                            num=self.obs_header[-1]["NAXIS3"],
                        )
                    )
                else:
                    self.obs_vel.append([None])

                self.obs_frequency.append([None])
                self.obs_wavelength.append([None])
                self.obs_iid.append(np.array(self.obs_header[-1]["TRANSL"].split(", ")))
                self.obs_i_iid.append(
                    np.array(self.obs_header[-1]["TRANSI"].split(", ")).astype(int)
                )

        if os.path.exists(spectra_path):
            self.files.append(self.survey + "_resampled_conf.csv")
            self.obs_iid.append(np.array(["C+ 1"]))
            self.obs_i_iid.append(np.array([0]))
            self.obs_spectra = pd.read_csv(spectra_path + self.survey + "_combined.csv")
            self.obs_spectra_resampled = pd.read_csv(
                spectra_path + self.survey + "_resampled.csv"
            )
            self.obs_spectra_mid = pd.read_csv(
                spectra_path + self.survey + "_resampled_conf.csv"
            )
            self.obs_data.append(self.obs_spectra_mid.Tmb.to_numpy())
            self.obs_lon.append(self.obs_spectra_mid.glon.to_numpy())
            self.obs_lat.append(np.zeros(1))
            self.obs_vel.append(self.obs_spectra_mid.Vel.to_numpy())
            self.obs_error_data.append(self.obs_spectra_mid.sigma.to_numpy())
            self.obs_error_conf_data.append(self.obs_spectra_mid.sigma_conf.to_numpy())
            self.obs_spectra_reduced = self.obs_spectra.loc[
                self.obs_spectra.glat == lat
            ]
            self.obs_spectra_resampled_reduced = self.obs_spectra_resampled.loc[
                self.obs_spectra_resampled.glat == lat
            ]

        return

    # Open spectra as a dataframe if included
    def get_obs_extent(self, filename=None, idx=None, kind="extent", verbose=False):
        """
        Return the extent of the observation with usable data (nonzero and
        non-NaN).
        Can return either the dimension values ('extent') or the dimension
        indeces ('index'), specified with `kind`.

        Returns tuple of (lon, lat, vel)
        """

        if filename is None and idx is None:
            print("ERROR: Please specify a filename or index")
            return
        elif not filename is None:
            idx = self.files.index(filename)
        else:
            filename = self.files[idx]

        if verbose:
            print(filename)

        data = self.obs_data[idx]

        if ".idl" in filename:
            lat = copy(self.obs_lat[idx])
            lon = copy(self.obs_lon[idx])
            extent = (lon, lat, None)
            i_extent = (np.arange(lon.size), np.zeros(0))
        elif self.obs_lat[idx].size == 1:
            lon = copy(self.obs_lon[idx])
            lat = copy(self.obs_lat[idx])
            vel = copy(self.obs_vel[idx])
            extent = (lon, lat, vel)
            i_extent = (np.arange(lon.size), np.zeros(1), np.arange(vel.size))
        elif data.ndim == 2:
            i_nan = np.isnan(data)
            lon = self.obs_lon[idx][~i_nan.all(0)]
            lat = self.obs_lat[idx][~i_nan.all(1)]
            extent = (lon, lat, None)
            i_extent = (~i_nan.all(1), ~i_nan.all(0))
        elif data.ndim == 3:
            i_nan = np.isnan(data)
            if self.survey in missions_2d:
                lon = self.obs_lon[idx][~i_nan.all(1).all(0)]
                lat = self.obs_lat[idx][~i_nan.all(2).all(0)]
                extent = (lon, lat, None)
                i_extent = (
                    np.arange(data.shape[0]),
                    ~i_nan.all(2).all(0),
                    ~i_nan.all(1).all(0),
                )
            else:
                vel = self.obs_vel[idx][~i_nan.all(2).all(1)]
                lon = self.obs_lon[idx][~i_nan.all(1).all(0)]
                lat = self.obs_lat[idx][~i_nan.all(2).all(0)]
                extent = (lon, lat, vel)
                i_extent = (
                    ~i_nan.all(2).all(1),
                    ~i_nan.all(2).all(0),
                    ~i_nan.all(1).all(0),
                )
        else:
            print("ERROR: Choose a valid filename.")
            return

        if kind == "extent":
            return extent
        elif kind == "index":
            return i_extent
        else:
            print("ERROR: Choose a valid kind ('extent' or 'index').")
            return

    def get_intensity(
        self,
        filename=None,
        idx=None,
        integrated=False,
        log=False,
        nan=True,
        trimmed=False,
        verbose=False,
    ):
        """
        This method will return the intensity data
        """

        if filename is None and idx is None:
            print("ERROR: Please specify a filename or index")
            return
        elif not filename is None:
            idx = self.files.index(filename)
        else:
            filename = self.files[idx]

        if verbose:
            print(filename)

        if integrated and not self.survey in (*missions_2d, "GOT_C+"):
            data = np.trapz(self.obs_data[idx], self.obs_vel[idx], axis=0)
        else:
            data = self.obs_data[idx]
        i_zero = data <= 0
        i_nan = np.isnan(data)

        data_temp = np.zeros_like(data)

        if log:
            data_temp[i_nan | i_zero] = np.nan
            data_temp[~(i_nan | i_zero)] = np.log10(data[~(i_nan | i_zero)])
        else:
            data_temp[i_nan] = np.nan
            data_temp[~i_nan] = data[~i_nan]

        if nan is False:
            index = np.isnan(data_temp)
            return data_temp[index]
        elif trimmed:
            index = np.ix_(*self.get_obs_extent(idx=idx, kind="index"))
            return data_temp[index]
        else:
            return data_temp
