"""
This file is part of the pyIMC Monte Carlo transport code (py-IMC)
It processes the output of the SOLEDGE code and generates a plasma profile for GITR
"""
import os
from typing import Tuple, Dict, Any
import netCDF4 as nc
import torch
import numpy as np
import collections

PlasmaParams = Tuple[float, float, float, float, float, float, float, float, float]
PlasmaProfile = np.ndarray

def plasma(model: torch.nn.Module, r: np.ndarray, z: np.ndarray) -> PlasmaProfile:
    """
        Args:
            model (Torch function): neural network model trained of SOLEDGE data
            r (float): radial coordinate
            z (float): vertical coordinate
        Returns:
            b_phi, b_r, b_r: magnetic field components
            T_e, n_e, v_e: electron temperature, density, and velocity
            t_i, n_i, v_i: ion temperature, density, and velocity
        """
    Inputs = collections.namedtuple('SOLEDGE', 'R Z')
    L = np.zeros((len(r),len(z),9))

    for i in range(len(r)):
        for j in range(len(z)):
            [b_phi,b_r,b_z,v_e,n_e,T_e,v_i,n_i,T_i]=[model(Inputs(r[i], z[j]))[1][k] for k in range(9)]
            L[i,j,:] = [b_phi,b_r,b_z,v_e,n_e,T_e,v_i,n_i,T_i]

    return L

def write_plasma_profile(profiles_filename: str, plasma_data: Dict[str, Any]) -> None:
    """Writes plasma profile data to a NetCDF file.
    Args:
        profiles_filename: The name of the NetCDF file to write.
        plasma_data: A dictionary containing the following keys: r, z, ne, te, ve, ni, ti, vi, br, bt, bz.
            The values associated with each key are 2D arrays of plasma parameters, where each row corresponds
            to a specific r value and each column corresponds to a specific z value.
    """
    if os.path.exists(profiles_filename):
        os.remove(profiles_filename)
    #  
    rootgrp = nc.Dataset(profiles_filename, "w", format="NETCDF4")
    nR, nZ = plasma_data['bz'].shape
    nr = rootgrp.createDimension("nR", nR)
    nz = rootgrp.createDimension("nZ", nZ)
    r = rootgrp.createVariable("gridR", "f8", ("nR"))
    z = rootgrp.createVariable("gridZ", "f8", ("nZ"))
    ve = rootgrp.createVariable("ve", "f8", ("nZ", "nR"))
    ne = rootgrp.createVariable("ne", "f8", ("nZ", "nR"))
    te = rootgrp.createVariable("te", "f8", ("nZ", "nR"))
    vi = rootgrp.createVariable("vi", "f8", ("nZ", "nR"))
    ni = rootgrp.createVariable("ni", "f8", ("nZ", "nR"))
    ti = rootgrp.createVariable("ti", "f8", ("nZ", "nR")) 
    vpar = rootgrp.createVariable("Vpara", "f8", ("nZ", "nR"))
    br = rootgrp.createVariable("br", "f8", ("nZ", "nR"))
    bt = rootgrp.createVariable("bt", "f8", ("nZ", "nR"))
    bz = rootgrp.createVariable("bz", "f8", ("nZ", "nR"))
    #
    r[:] = plasma_data['r']
    z[:] = plasma_data['z']
    ve[:] = plasma_data['ve'].T
    ne[:] = plasma_data['ne'].T
    te[:] = plasma_data['te'].T
    vi[:] = plasma_data['vi'].T   
    ni[:] = plasma_data['ni'].T
    ti[:] = plasma_data['ti'].T
    
    vpar[:] = 0.0
    br[:] = plasma_data['br'].T
    bt[:] = 0.0
    bz[:] = plasma_data['bz'].T
    rootgrp.close()

#    # Write data to file
"""
Args:
    model: A PyTorch neural network model trained on SOLEDGE data.
    r: An array of radial coordinates.
    z: An array of vertical coordinates.
    profiles_filename: The name of the NetCDF file to write.
"""
nR,nZ=100,200
r=np.linspace(1.8,3.5,nR)
z=np.linspace(-1,1,nZ)
model = torch.load('data.pt')
profiles_filename = 'profiles.nc'
plasma_data = {}
plasma_values = plasma(model, r, z)
# print(plasma(model, r, z).shape)
plasma_data['r'] = r
plasma_data['z'] = z
plasma_data['ve'] = plasma_values[:, :, 3]
plasma_data['ne'] = plasma_values[:, :, 4]
plasma_data['te'] = plasma_values[:, :, 5]

plasma_data['vi'] = plasma_values[:, :, 6]

plasma_data['ni'] = plasma_values[:, :, 7]
plasma_data['ti'] = plasma_values[:, :, 8]

plasma_data['br'] = plasma_values[:, :, 1]
plasma_data['bt'] = np.zeros_like(plasma_data['br'])
plasma_data['bz'] = plasma_values[:, :, 2]
#
write_plasma_profile(profiles_filename, plasma_data)