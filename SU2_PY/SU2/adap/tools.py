#!/usr/bin/env python

## \file tools.py
#  \brief Useful functions for configuring mesh adaptation parameters
#  \author Victorien Menier, Brian Mungu\'ia
#  \version 8.0.1 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import numpy as np
from itertools import islice
import os

def get_mesh_sizes(config):
    """Get prescribed mesh complexities, i.e. desired mesh sizes"""
    return config['ADAP_SIZES'].strip('()').split(',')

def get_mesh_hmaxs(config):
    """Get prescribed mesh hmax, i.e. desired maximum element sizes"""
    if 'ADAP_HMAXS' in config:
        return config['ADAP_HMAXS'].strip('()').split(',')
    elif 'ADAP_HMAX' in config:
        nExt_iter = len(config['ADAP_SIZES'].strip('()').split(','))
        hmaxs = []
        for i in range(nExt_iter):
            hmaxs.append(config['ADAP_HMAX'])
        return hmaxs
    else:
        raise KeyError('Either ADAP_HMAX or ADAP_HMAXS needs to be specified')

def get_mesh_hmins(config):
    """Get prescribed mesh hmin, i.e. desired minimum element sizes"""
    if 'ADAP_HMINS' in config:
        return config['ADAP_HMINS'].strip('()').split(',')
    elif 'ADAP_HMIN' in config:
        nExt_iter = len(config['ADAP_SIZES'].strip('()').split(','))
        hmins = []
        for i in range(nExt_iter):
            hmins.append(config['ADAP_HMIN'])
        return hmins
    else:
        raise KeyError('Either ADAP_HMIN or ADAP_HMINS needs to be specified')

def get_mesh_armaxs(config):
    """Get prescribed mesh aspect ratio"""
    if 'ADAP_ARMAXS' in config:
        return config['ADAP_ARMAXS'].strip('()').split(',')
    elif 'ADAP_ARMAX' in config:
        nExt_iter = len(config['ADAP_SIZES'].strip('()').split(','))
        armaxs = []
        for i in range(nExt_iter):
            armaxs.append(config['ADAP_ARMAX'])
        return armaxs
    else:
        raise KeyError('Either ADAP_ARMAX or ADAP_ARMAXS needs to be specified')

def get_mesh_norms(config):
    """Get prescribed mesh Lp norms"""
    if 'ADAP_NORMS' in config:
        return config['ADAP_NORMS'].strip('()').split(',')
    elif 'ADAP_NORM' in config:
        nExt_iter = len(config['ADAP_SIZES'].strip('()').split(','))
        norms = []
        for i in range(nExt_iter):
            norms.append(config['ADAP_NORM'])
        return norms
    else:
        raise KeyError('Either ADAP_NORM or ADAP_NORMS needs to be specified')

def get_mesh_size(mesh):
    """Get mesh size info from a python mesh structure"""
    elts = {'xy':'Vertices', 'xyz':'Vertices', 'Tetrahedra':'Tetrahedra', 'Triangles':'Triangles', 'Edges':'Edges'}
    nelts = {}
    for k,v in elts.items():
        if v not in mesh.keys(): continue

        nbe = len(mesh[v])
        if nbe > 0: nelts.update({f'{v}': nbe})

    return nelts


def get_sub_iterations(config):
    """Get number of adaptation iterations for each mesh complexity"""
    return config['ADAP_SUBITER'].strip('()').split(',')

def get_adj_iter(config):
    """Get number of adjoint solver iterations for each mesh complexity"""
    if 'ADAP_ADJ_ITER' in config:
        return config['ADAP_ADJ_ITER'].strip('()').split(',')
    else:
        nExt_iter = len(config['ADAP_SIZES'].strip('()').split(','))
        ext_iter = []
        for i in range(nExt_iter):
            ext_iter.append(config['ITER'])
        return ext_iter

def get_flow_iter(config):
    """Get number of primal solver iterations for each mesh complexity"""
    if 'ADAP_FLOW_ITER' in config:
        return config['ADAP_FLOW_ITER'].strip('()').split(',')
    else:
        nExt_iter = len(config['ADAP_SIZES'].strip('()').split(','))
        flow_iter = []
        for i in range(nExt_iter):
            flow_iter.append(config['ITER'])
        return flow_iter

def get_flow_cfl(config):
    """Get initial CFL number for each mesh complexity"""
    if 'ADAP_FLOW_CFL' in config:
        return config['ADAP_FLOW_CFL'].strip('()').split(',')
    else:
        ncfl = len(config['ADAP_SIZES'].strip('()').split(','))
        cfl = []
        for i in range(ncfl):
            cfl.append(config['CFL_NUMBER'])
        return cfl

def get_adj_cfl(config):
    """Get adjoint CFL number for each mesh complexity"""
    if 'ADAP_ADJ_CFL' in config:
        return config['ADAP_ADJ_CFL'].strip('()').split(',')
    else:
        ncfl = len(config['ADAP_SIZES'].strip('()').split(','))
        cfl = []
        for i in range(ncfl):
            cfl.append(config['CFL_NUMBER'])
        return cfl

def set_cfl(config, cfl_iSiz):
    """Set CFL parameters for current mesh complexity"""
    config.CFL_NUMBER = float(cfl_iSiz)
    if 'CFL_ADAPT' in config:
        if config['CFL_ADAPT'] == 'YES':
            cfl_params = [float(x) for x in config['CFL_ADAPT_PARAM'].strip('()').split(',')]
            cfl_params[2] = cfl_iSiz

            cfl_param_str = '( '
            for i, param in enumerate(cfl_params):
                cfl_param_str = f"{cfl_param_str} {param}"
                if (i == len(cfl_params)-1): cfl_param_str = f"{cfl_param_str} )"
                else: cfl_param_str = f"{cfl_param_str}, "
            cfl_param_str = config['CFL_ADAPT_PARAM']

def get_adap_sensors(config):
    """Get adaptation sensors"""
    return config['ADAP_SENSOR'].replace(' ','').strip('()').split(',')

def get_pyadap_options(config):
    """Get pyadap options"""
    pyadap_dict = {}
    pyadap_dict['ADAP_SIZES'] = get_mesh_sizes(config)
    pyadap_dict['ADAP_SUBITER'] = get_sub_iterations(config)
    pyadap_dict['ADAP_HMAXS'] = get_mesh_hmaxs(config)
    pyadap_dict['ADAP_HMINS'] = get_mesh_hmins(config)
    pyadap_dict['ADAP_NORMS'] = get_mesh_norms(config)
    pyadap_dict['ADAP_ARMAXS'] = get_mesh_armaxs(config)
    pyadap_dict['ADAP_ADJ_ITER'] = get_adj_iter(config)
    pyadap_dict['ADAP_FLOW_ITER'] = get_flow_iter(config)
    pyadap_dict['ADAP_FLOW_CFL'] = get_flow_cfl(config)
    pyadap_dict['ADAP_ADJ_CFL'] = get_adj_cfl(config)
    pyadap_dict['ADAP_SENSOR'] = get_adap_sensors(config)

    return pyadap_dict

def set_flow_config_ini(config, cur_solfil, pyadap_dict):
    """Set primal config for initial solution"""
    config.CONV_FILENAME       = 'history'
    config.RESTART_FILENAME    = cur_solfil
    config.HISTORY_OUTPUT      = ['ITER', 'RMS_RES', 'AERO_COEFF', 'FLOW_COEFF', 'CFL_NUMBER']
    config.MATH_PROBLEM        = 'DIRECT'
    config.WRT_RESTART_COMPACT = 'NO'
    if 'GOAL' in pyadap_dict['ADAP_SENSOR']:
        config.VOLUME_OUTPUT  = 'COORDINATES, SOLUTION, PRIMITIVE, CFL_NUMBER, AUXILIARY, RESIDUAL'
        config.COMPUTE_METRIC = 'NO'
    else:
        config.VOLUME_OUTPUT   = 'COORDINATES, SOLUTION, PRIMITIVE, CFL_NUMBER, AUXILIARY, RESIDUAL, METRIC, GRADIENT_ADAPT'
        config.COMPUTE_METRIC  = 'YES'
        config.ADAP_COMPLEXITY = int(pyadap_dict['ADAP_SIZES'][0])
        config.ADAP_HMAX       = float(pyadap_dict['ADAP_HMAXS'][0])
        config.ADAP_HMIN       = float(pyadap_dict['ADAP_HMINS'][0])
        config.ADAP_NORM       = float(pyadap_dict['ADAP_NORMS'][0])
        config.ADAP_ARMAX      = float(pyadap_dict['ADAP_ARMAXS'][0])

def set_adj_config_ini(config, cur_solfil, cur_solfil_adj, pyadap_dict):
    """Set adjoint config for initial solution"""
    config.CONV_FILENAME        = 'history_adj'
    config.RESTART_ADJ_FILENAME = cur_solfil_adj
    config.SOLUTION_FILENAME    = cur_solfil
    config.RESTART_FILENAME     = cur_solfil
    config.WRT_RESTART_COMPACT  = 'NO'
    config.MATH_PROBLEM         = 'DISCRETE_ADJOINT'
    config.VOLUME_OUTPUT        = 'COORDINATES, SOLUTION, PRIMITIVE, CFL_NUMBER, RESIDUAL, METRIC'
    config.HISTORY_OUTPUT       = ['ITER', 'RMS_RES', 'SENSITIVITY']
    config.COMPUTE_METRIC       = 'YES'
    config.ADAP_COMPLEXITY = int(pyadap_dict['ADAP_SIZES'][0])
    config.ADAP_HMAX       = float(pyadap_dict['ADAP_HMAXS'][0])
    config.ADAP_HMIN       = float(pyadap_dict['ADAP_HMINS'][0])
    config.ADAP_NORM       = float(pyadap_dict['ADAP_NORMS'][0])
    config.ADAP_ARMAX      = float(pyadap_dict['ADAP_ARMAXS'][0])


def update_flow_config(config, cur_meshfil, cur_solfil, cur_solfil_ini, pyadap_dict, iter, subiter):
    """Set primal config for current solution"""
    if subiter == int(pyadap_dict['ADAP_SUBITER'][iter])-1 and iter < len(pyadap_dict['ADAP_SIZES'])-1:
        iter += 1
    else:
        iter = iter
    config.MESH_FILENAME     = cur_meshfil
    config.SOLUTION_FILENAME = cur_solfil_ini
    config.RESTART_FILENAME  = cur_solfil
    config.ITER              = int(pyadap_dict['ADAP_FLOW_ITER'][iter])
    if 'GOAL' not in pyadap_dict['ADAP_SENSOR']:
        config.ADAP_COMPLEXITY = int(pyadap_dict['ADAP_SIZES'][iter])
        config.ADAP_HMAX       = float(pyadap_dict['ADAP_HMAXS'][iter])
        config.ADAP_HMIN       = float(pyadap_dict['ADAP_HMINS'][iter])
        config.ADAP_NORM       = float(pyadap_dict['ADAP_NORMS'][iter])
        config.ADAP_ARMAX      = float(pyadap_dict['ADAP_ARMAXS'][iter])  

    set_cfl(config, float(pyadap_dict['ADAP_FLOW_CFL'][iter]))

def update_adj_config(config, cur_meshfil, cur_solfil, cur_solfil_adj, cur_solfil_adj_ini, pyadap_dict, iter, subiter):
    """Set adjoint config for current solution"""
    if subiter == int(pyadap_dict['ADAP_SUBITER'][iter])-1 and iter < len(pyadap_dict['ADAP_SIZES'])-1:
        iter += 1
    else:
        iter = iter
    config.MESH_FILENAME         = cur_meshfil
    config.RESTART_ADJ_FILENAME  = cur_solfil_adj
    config.SOLUTION_ADJ_FILENAME = cur_solfil_adj_ini
    config.SOLUTION_FILENAME     = cur_solfil
    config.RESTART_FILENAME      = cur_solfil
    config.ITER                  = int(pyadap_dict['ADAP_ADJ_ITER'][iter])
    config.ADAP_COMPLEXITY = int(pyadap_dict['ADAP_SIZES'][iter])
    config.ADAP_HMAX       = float(pyadap_dict['ADAP_HMAXS'][iter])
    config.ADAP_HMIN       = float(pyadap_dict['ADAP_HMINS'][iter])
    config.ADAP_NORM       = float(pyadap_dict['ADAP_NORMS'][iter])
    config.ADAP_ARMAX      = float(pyadap_dict['ADAP_ARMAXS'][iter])  

def set_mmg_config(config_su2, dim):
    """Load parameters from the SU2 config file into the AMG config dict"""
    config_mmg = dict()

    if 'ADAP_HGRAD' in config_su2: config_mmg['hgrad'] = float(config_su2['ADAP_HGRAD'])
    if 'ADAP_ANGLE' in config_su2: config_mmg['ar']    = int(config_su2['ADAP_ANGLE'])

    config_mmg['dim']     = int(dim)
    config_mmg['hmax']    = float(get_mesh_hmaxs(config_su2)[0])
    config_mmg['hmin']    = float(get_mesh_hmins(config_su2)[0])
    config_mmg['mmg_log'] = 'mmg.out'
    config_mmg['mmg_err'] = 'mmg.err'

    if '(' in config_su2['ADAP_HAUSD']:
        config_mmg['hausd'] = {}
        parameters = config_su2['ADAP_HAUSD'].lstrip('(').rstrip(')').split(',')

        if len(parameters) % 2:
            raise KeyError('Missing values in ADAP_HAUSD!')
        else:
            for i_par in range(0, len(parameters), 2):
                config_mmg['hausd'][parameters[i_par].strip(' ')] = float(parameters[i_par+1])
    else:
        config_mmg['hausd'] = float(config_su2['ADAP_HAUSD'])
    
    return config_mmg


def update_mmg_config(config_mmg, pyadap_dict, iter):

    config_mmg['hmax'] = float(pyadap_dict['ADAP_HMAXS'][iter])
    config_mmg['hmin'] = float(pyadap_dict['ADAP_HMINS'][iter])

def print_adap_options(config):
    """Print options used for mesh adaptation"""
    pad = 0
    for key, value in config.items():
        if 'ADAP_' in key:
            pad = max(len(key), pad)

    prt = '\nMesh adaptation options:\n'
    for key, value in config.items():
        if 'ADAP_' in key:
            prt += f'{key:<{pad}} : {value}\n'
    return prt

def get_su2_dim(filename):
    """Get dimension from su2 mesh"""
    meshfile = open(filename,'r')

    def mesh_readlines(n_lines=1):
        fileslice = islice(meshfile,n_lines)
        return list(fileslice)

    dim = -1

    keepon = True
    while keepon:

        line = mesh_readlines()

        if not line:
            keepon = False
            break

        # fix white space
        line = line[0]
        line = line.replace('\t',' ')
        line = line.replace('\n',' ')

        # skip comments
        if line[0] == '%':
            pass

        # number of dimensions
        elif 'NDIME=' in line:
            # save to SU2_MESH data
            dim = int( line.split('=')[1].strip() )
            keepon = False

    return dim

def get_su2_npoin(filename):
    """Get number of points from su2 mesh"""
    meshfile = open(filename,'r')

    def mesh_readlines(n_lines=1):
        fileslice = islice(meshfile,n_lines)
        return list(fileslice)

    npoin = -1

    keepon = True
    while keepon:

        line = mesh_readlines()

        if not line:
            keepon = False
            break

        # fix white space
        line = line[0]
        line = line.replace('\t',' ')
        line = line.replace('\n',' ')

        # skip comments
        if line[0] == '%':
            pass

        # number of dimensions
        elif 'NPOIN=' in line:
            # save to SU2_MESH data
            npoin = int( line.split('=')[1].strip() )
            keepon = False

    return npoin

def merge_sol(mesh0, mesh1):
    """Merge 2 solutions (e.g. primal and adjoint)"""
    mesh0['solution'] = np.hstack((mesh0['solution'], \
                                   mesh1['solution'])).tolist()
    mesh0['solution_tag'] = np.hstack((np.array(mesh0['solution_tag']), \
                                       np.array(mesh1['solution_tag']))).tolist()

def split_adj_sol(mesh):
    """Separate adjoint solution from primal"""
    nsol = len(mesh['solution_tag'])

    adj_sol = dict()

    for i in range(nsol):
        if 'Adjoint' in mesh['solution_tag'][i]:
            iAdj = i

            adj_sol['solution'] = np.delete(np.array(mesh['solution']), np.s_[0:iAdj], axis=1).tolist()
            adj_sol['solution_tag'] = np.delete(np.array(mesh['solution_tag']), np.s_[0:iAdj], axis=0).tolist()

            if 'xyz' in mesh:
                adj_sol['xyz'] = mesh['xyz']
            elif 'xy' in mesh:
                adj_sol['xy'] = mesh['xy']

            adj_sol['dimension'] = mesh['dimension']

            mesh['solution'] = np.delete(np.array(mesh['solution']), np.s_[iAdj:nsol], axis=1).tolist()
            mesh['solution_tag'] = np.delete(np.array(mesh['solution_tag']), np.s_[iAdj:nsol], axis=0).tolist()

            break

    return adj_sol

def create_sensor(solution, sensor_tags):
    """
    Store desired sensor for adaptation

    Returns array of scalars for MACH/PRES, or array of tensors
    for GOAL
    """
    Dim = solution['dimension']
    Sol = np.array(solution['solution'])

    nMet = 3*(Dim-1)
    sensor = Sol[:,-nMet:]
    sensor = np.array(sensor).reshape((len(sensor),nMet))

    sensor_wrap = dict()

    sensor_wrap['solution_tag'] = '-'.join(sensor_tags)
    sensor_wrap['xyz'] = solution['xyz']

    sensor_wrap['dimension'] = solution['dimension']
    sensor_wrap['solution']  = sensor

    return sensor_wrap

def print_adap_table(iter, subiter, pyadap_dict, meshDict):
    """Print adapted mesh sizes to a table"""
    dim = meshDict['Dim']
    sizes = pyadap_dict['ADAP_SIZES']
    nsubiter = int(pyadap_dict['ADAP_SUBITER'][iter])

    #--- Header
    if iter == 0 and subiter == 0:
        print('+=================================================================+')
        if dim == 2:
            print('|   Iter   |   Size   | Sub-iter |   Vert   |   Tria   |   Edge   |')
        else:
            print('|   Iter   |   Size   | Sub-iter |   Vert   |   Tetr   |   Tria   |')
        print('+=================================================================+')

    #--- Data
    line = None
    if subiter == 0:
        size = int(sizes[iter])
        line = f'|    {iter:<2}    | {size:<8} '
    else:
        pad_nul = ' '*10
        line = f'|{pad_nul}|{pad_nul}'

    nvert = len(meshDict['Vertices'])
    if dim == 2:
        nedge = sum(len(faces) for faces in meshDict['Edges'].values())
        ntria = len(meshDict['Triangles'])
        line = f'{line}|    {subiter:<2}    | {nvert:<8} | {ntria:<8} | {nedge:<8} |'
    else:
        ntria = sum(len(faces) for faces in meshDict['Triangles'].values())
        ntetr = len(meshDict['Tetrahedra'])
        line = f'{line}|    {subiter:<2}    | {nvert:<8} | {ntetr:<8} | {ntria:<8} |'
    print(line)

    if subiter == nsubiter-1:
        if iter != len(sizes)-1:
            print('+-----------------------------------------------------------------+')
        else:
            print('+=================================================================+')

def plot_results(history_format, filename, iter, npoin):
    """Write primal results to a Tecplot or CSV file"""
    default_spacing = 16
    indent_spacing  = 0

    #--- Format and file name
    if (history_format == 'TECPLOT'):
        solname  = 'history.dat'
        indent_spacing += 10
    else:
        solname  = 'history.csv'
    indent_spacing = ' '*indent_spacing

    #--- Write header on first adaptive iteration
    if iter == 0:
        #--- Get header from solution history
        header = ''

        if (history_format == 'TECPLOT'):
            header     = 'VARIABLES='
            headerline = 1
        else:
            headerline = 0

        with open(solname, 'rb') as f:
            for i, line in enumerate(f):
                if i == headerline:
                    break

        header = f"{header}\"Adap_Iter\", \"NDOFs\", {line.decode('ascii')}"

        plotfile = open(filename,'w')
        plotfile.write(header)

    # --- Append data on all other iterations
    else:
        plotfile = open(filename,'a')

    #--- Get data from last line of file
    with open(solname, 'rb') as f:
        f.seek(-2, os.SEEK_END)
        while f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        last_line = f.readline().decode('ascii')

    plotfile.write(f'{indent_spacing}{iter}, {npoin}, {last_line}')
    plotfile.close()
