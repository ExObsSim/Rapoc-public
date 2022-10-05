import astropy.units as u
import numpy as np
import struct
import os
import glob

from .loader import FileLoader


class DACEFileLoader(FileLoader):
    '''
    File loader for DACE files. It works for opacities files in the binary format.
    These files can be downloaded from `DACE website <https://dace.unige.ch/opacityDatabase/?>`_.

    The scripts used to read the DACE files are inspired by `petitRADTRANS <https://petitradtrans.readthedocs.io/en/latest/content/opa_add.html#owtoprt>`_.

    Notes
    -----
    This class is implicitly called by the model when initialised if the matched input is used.

    Examples
    ---------
    Let's build the Rosseland method using an DACE files as input.
    The binary files of the desired molecule should be stored in a `dace_dir`.

    >>> from rapoc import Rosseland
    >>> daceFileList = '/dace_dir'
    >>> model = Rosseland(input_data=daceFileList)

    Now the model is ready to be used

    Warnings
    --------
    This method reads all the binary data in the indicated folder: the user must be cautios to indicate a directory that
    contains only the binary files of the desired molecule.

    '''

    def _open(self):
        # Read the fiducial petitRADTRANS wavelength grid
        path_to_files = self.filename
        filelist = [f for f in os.listdir(path_to_files) if f.endswith('.bin')]
        opened_file = {'p': [], 't': [], 'wn': None, 'opac': None}

        # prepare pressures and temperature grids
        for file in filelist:
            t = file[16:21]
            t = t.lstrip('0')
            opened_file['t'].append(t)
            p = file[22:26]
            opened_file['p'].append(p)
        opened_file['t'] = list(set(opened_file['t']))
        opened_file['p'] = list(set(opened_file['p']))

        # prepare wavelength grids
        opa = _read_bin_single(os.path.join(path_to_files, filelist[0]))
        wl_start = int(filelist[0][4:9])
        wl_end = int(filelist[0][10:15])
        wlen = np.linspace(wl_start, wl_end, len(opa))
        opened_file['wn'] = wlen

        # prepare opac grids
        opened_file['opac'] = np.empty([len(opened_file['p']), len(opened_file['t']), len(opened_file['wn'])])

        # populate opac grid
        print('loading data')
        from tqdm import tqdm
        for i, p in enumerate(tqdm(opened_file['p'])):
            filelist_p = [f for f in filelist if p in f]
            for j, t in enumerate(opened_file['t']):
                filelist_t = [f for f in filelist_p if str(t) in f]
                file = filelist_t[0]
                opa = _read_bin_single(os.path.join(path_to_files, file))
                opened_file['opac'][i, j] = opa

        # convert pressure an temperature grids
        list_p = []
        for p in opened_file['p']:
            if p[0] == 'p':
                sign = 1
            elif p[0] == 'n':
                sign = -1
            pres = float(p[1]) + float(p[2:]) / 100
            list_p.append(sign * pres)
        opened_file['p'] = np.array(list_p)
        opened_file['p'] = np.power(10, opened_file['p'])

        list_t = [float(t) for t in opened_file['t']]
        opened_file['t'] = list_t

        return opened_file

    def _close(self, opened_file):
        pass

    def _read_molecule_name(self, opened_file):
        mol_name = None
        path_to_files = self.filename
        filelist = os.listdir(path_to_files)
        for f in filelist:
            if f.endswith('.ref'):
                import re
                mol_name_list = []
                mol_code = f.split('_')[0]
                mol_elem = mol_code.split('-')
                for elem in mol_elem:
                    tmp = re.findall('\d+|\D+', elem)
                    for i in range(len(tmp)):
                        try:
                            tmp[i] = int(tmp[i])
                        except ValueError:
                            mol_name_list.append(tmp[i])
                            try:
                                mol_name_list.append(tmp[i + 1])
                            except IndexError:
                                continue
                mol_name = ''.join(mol_name_list)
        return mol_name

    def _get_mol_mass(self, mol):
        from molmass import Formula
        f = Formula(mol)
        self.mol_mass = (f.mass * u.u / u.ct).to(u.g / u.ct)
        return self.mol_mass

    def _read_pressure_grid(self, opened_file):
        pressure_grid = np.array(opened_file['p']) * u.bar
        return pressure_grid.si

    def _read_temperature_grid(self, opened_file):
        temperature_grid = np.array(opened_file['t']) * u.K
        return temperature_grid.si

    def _read_wavenumber_grid(self, opened_file):
        return opened_file['wn'] / u.cm

    def _read_opacities(self, opened_file):
        opac = np.array(opened_file['opac']) * u.cm ** 2 / u.g
        return opac.si


def _read_bin_single(filename):
    """ Read a binary opacity world file.
    """
    return np.fromfile(filename, dtype=np.float32)
