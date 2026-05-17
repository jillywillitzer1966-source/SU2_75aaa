#!/usr/bin/env python

## \file mmg.py
#  \brief python script for running mesh adaptation using the MMG Inria library
#  \author Victorien Menier, Brian Mungu\'ia
#  \version 7.3.0 "Blackbird"
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

import os, shutil, copy, time

from .. import io as su2io
from .. import adap as su2adap
from ..run.interface import CFD as SU2_CFD

def mmg(config, meditformat, runCFD = True):
    """
    Runs the a mesh adaptation loop with the MMG library.

    Inputs:
        config - an SU2 config object
    """

    #--- Check config options related to mesh adaptation

    pyadap_options = [ 'ADAP_SIZES', 'ADAP_SUBITER', 'ADAP_HGRAD', 'ADAP_RESIDUAL_REDUCTION', 
                      'ADAP_FLOW_ITER', 'ADAP_ADJ_ITER', 'ADAP_CFL', 'ADAP_HAUSD', 'ADAP_HMAXS', 
                      'ADAP_HMINS', 'ADAP_NORMS', 'ADAP_ARMAXS' ]
    required_options = [ 'ADAP_SIZES', 'ADAP_SUBITER', 'MESH_FILENAME', 
                        'RESTART_SOL', 'MESH_OUT_FILENAME' ]

    if not all (opt in config for opt in required_options):
        err = '\n\n## ERROR : Missing options: \n'
        for opt in required_options:
            if not opt in config:
                err += opt + '\n'
        raise AttributeError(err)
    
    #--- Medit format
    if meditformat:
        medit_ext = '.meshb'
    else:
        medit_ext = '.mesh'
    
    #--- NEMO solver check
    if 'NEMO' in config.SOLVER:
        nemo = True
    else:
        nemo = False

    #--- Print adap options

    print(su2adap.print_adap_options(config))

    #--- Extracting ADAP options and check for sensors

    pyadap_dict = su2adap.get_pyadap_options(config) 

    if len(pyadap_dict['ADAP_SIZES']) != len(pyadap_dict['ADAP_SUBITER']):
        raise ValueError(f'Inconsistent number of mesh sizes and sub-iterations. {len(pyadap_dict['ADAP_SIZES'])} mesh sizes and {len(pyadap_dict['ADAP_SUBITER'])} sub-iterations provided.')

    sensor_avail = ['GOAL', 'MACH', 'PRESSURE', 'TEMPERATURE', 'TEMPERATURE_VE', 'ENERGY', 'ENERGY_VE', 'DENSITY', 'TOTALPRESSURE']

    for sensor in pyadap_dict['ADAP_SENSOR']:
        if sensor not in sensor_avail:
            raise ValueError(f'Unknown adaptation sensor {sensor}. Available options are {sensor_avail}.')

    gol = 'GOAL' in pyadap_dict['ADAP_SENSOR']

    #--- Change current directory

    warn = True
    base_dir = os.getcwd()
    adap_dir = './adap'
    flow_dir = './Flows'
    surfflow_dir = './Flows_surf'

    if os.path.exists(adap_dir):
        print('./adap exists. Removing old mesh adaptation in 10s.')
        if warn : time.sleep(10)
        shutil.rmtree(adap_dir)
        print(f'The {adap_dir} folder was deleted.')
    
    if os.path.exists(flow_dir):
        print(flow_dir+' exists. Removing old Flows in 10s.')
        shutil.rmtree(flow_dir)
        print(f'The {flow_dir} folder was deleted.')
    
    if os.path.exists(surfflow_dir):
        print(surfflow_dir+' exists. Removing old Surface Flows in 10s.')
        shutil.rmtree(surfflow_dir)
        print(f'The {surfflow_dir} folder was deleted.')


    os.makedirs(flow_dir)
    os.makedirs(surfflow_dir)
    su2adap.run_command(f'ln -s ../adap/ite{0}/'+config.VOLUME_FILENAME+'.vtu ./Flows/'+config.VOLUME_FILENAME+'_'+str(0).zfill(5)+'.vtu')
    su2adap.run_command(f'ln -s ../adap/ite{0}/'+config.SURFACE_FILENAME+'.vtu ./Flows_surf/'+config.SURFACE_FILENAME+'_'+str(0).zfill(5)+'.vtu')
    dir = f'{adap_dir}/ite0'
    os.makedirs(dir)
    os.chdir(dir)
    os.symlink(os.path.join(base_dir, config.MESH_FILENAME), config.MESH_FILENAME)

    meshfil = config['MESH_FILENAME']

    #--- Format of history file

    history_format = config.TABULAR_FORMAT
    if (history_format == 'TECPLOT'):
        history_filename = os.path.join(base_dir, 'history_adap.dat')
    else:
        history_filename = os.path.join(base_dir, 'history_adap.csv')

    #--- Get mesh dimension

    dim = su2adap.get_su2_dim(meshfil)
    if ( dim != 2 and dim != 3 ):
        raise ValueError('Wrong dimension number.')

    #--- MMG parameters

    config_mmg = su2adap.set_mmg_config(config, dim)

    #--- Compute initial solution if needed, else link current files

    config_cfd = copy.deepcopy(config)
    config_cfd_ad = copy.deepcopy(config)
    for opt in pyadap_options:
        config_cfd.pop(opt, None)
        config_cfd_ad.pop(opt, None)

    #--- Check config for filenames if restarting
    restart = config['RESTART_SOL'] == 'YES'
    if restart:
        required_options = ['SOLUTION_FILENAME']
        if gol: required_options.append('SOLUTION_ADJ_FILENAME')
        if not all (opt in config for opt in required_options):
            err = 'RESTART_SOL is set to YES, but the solution is missing:\n'
            for opt in required_options:
                if not opt in config:
                    err += opt + '\n'
            raise ValueError(err)

        os.symlink(os.path.join(base_dir, config.SOLUTION_FILENAME), config.SOLUTION_FILENAME)

        print('\nInitial CFD solution is provided.')

    else:
        print('\nRunning initial CFD solution.')

    #--- Setting up initial flow config
    if not gol:
        if '.csv' in config.RESTART_FILENAME:
            config_cfd.READ_BINARY_RESTART = 'NO'
            sol_ext_cfd = '.csv'
            config_cfd.OUTPUT_FILES = ['RESTART_CSV','PARAVIEW','SURFACE_PARAVIEW']
        else:
            sol_ext_cfd = '.dat'
            config_cfd.OUTPUT_FILES = ['RESTART','PARAVIEW','SURFACE_PARAVIEW']

    if gol:
        if '.csv' in config.RESTART_FILENAME:
            config_cfd.READ_BINARY_RESTART = 'NO'
            sol_ext_cfd = '.csv'
            config_cfd.OUTPUT_FILES = ['RESTART_CSV','PARAVIEW','SURFACE_PARAVIEW']
            config_cfd_ad.READ_BINARY_RESTART = 'NO'
            sol_ext_cfd_ad = '.csv'
            config_cfd_ad.OUTPUT_FILES = ['RESTART_CSV','PARAVIEW','SURFACE_PARAVIEW']   
        else:
            sol_ext_cfd = '.dat'
            config_cfd.OUTPUT_FILES = ['RESTART','PARAVIEW','SURFACE_PARAVIEW']
            sol_ext_cfd_ad = '.dat'
            config_cfd_ad.OUTPUT_FILES = ['RESTART','PARAVIEW','SURFACE_PARAVIEW']

    solfil  = f'restart_flow{sol_ext_cfd}'
    su2adap.set_flow_config_ini(config_cfd, solfil, pyadap_dict)

    try: # run with redirected outputs
        #--- Run a single iteration of the flow if restarting to get history info
        if restart:
            config_cfd.ITER = 1

        with su2io.redirect.output('su2.out'): SU2_CFD(config_cfd)

        #--- Set RESTART_SOL=YES for runs after adaptation
        if not nemo and sol_ext_cfd != '.csv':     # SU2 does not interpolate ASCII files for restarts
            config_cfd.RESTART_SOL = 'YES' 

        if gol:
            adjsolfil = f'restart_adj{sol_ext_cfd_ad}'
            su2adap.set_adj_config_ini(config_cfd_ad, solfil, adjsolfil, pyadap_dict)

            #--- If restarting, check for the existence of an adjoint restart
            if restart:
                adjsolfil_ini = config_cfd_ad.SOLUTION_ADJ_FILENAME
                func_name          = config.OBJECTIVE_FUNCTION
                suffix             = su2io.get_adjointSuffix(func_name)

                #--- Run an adjoint if the solution file doesn't exist
                if not (os.path.exists(os.path.join(base_dir, su2io.add_suffix(adjsolfil_ini, suffix)))):
                    config_cfd_ad.ITER        = config.ITER
                    config_cfd_ad.RESTART_SOL = 'NO'

                    print('Running initial adjoint CFD solution.')

                #--- Otherwise just compute the metric
                else:
                    os.symlink(os.path.join(base_dir, adjsolfil_ini), adjsolfil_ini)
                    config_cfd_ad.ITER = 0

                    print('Initial adjoint CFD solution is provided.')

            else:
                print('Running initial adjoint CFD solution.')

            with su2io.redirect.output('su2.out'): SU2_CFD(config_cfd_ad)

    except:
        raise

    #--- Check existence of initial mesh, solution

    required_files = [meshfil, solfil]

    if not all (os.path.exists(fil) for fil in required_files):
        err = "Can't find the following files:\n"
        for fil in required_files:
            if not os.path.exists(fil):
                err += fil + '\n'
        raise Exception(err)

    #--- Start adaptive loop

    global_iter = 0

    #--- Print convergence history

    npoin = su2adap.get_su2_npoin(meshfil)
    su2adap.plot_results(history_format, history_filename, global_iter, npoin)

    print('\nStarting mesh adaptation process.\n')

    nSiz = len(pyadap_dict['ADAP_SIZES'])
    for iSiz in range(nSiz):
        nSub = int(pyadap_dict['ADAP_SUBITER'][iSiz])
        su2adap.update_mmg_config(config_mmg, pyadap_dict, iSiz)
        for iSub in range(nSub):
            
            global_iter += 1

            os.symlink(f'../adap/ite{global_iter}/'+config_cfd.VOLUME_FILENAME+'.vtu', '../../Flows/'+config_cfd.VOLUME_FILENAME+'_'+str(global_iter).zfill(5)+'.vtu')
            os.symlink(f'../adap/ite{global_iter}/'+config_cfd.SURFACE_FILENAME+'.vtu', '../../Flows_surf/'+config_cfd.SURFACE_FILENAME+'_'+str(global_iter).zfill(5)+'.vtu')

            #--- Instantiating the mesh converter
            fileconverter = su2adap.MeshSolConverter()

            #--- Load and read .su2 mesh and dump .mesh file

            fileconverter.SU2ToMeditMesh(meshfil, meshfil.replace('.su2', medit_ext))

            #--- Load and read .csv solution and dump .sol file

            fileconverter.SU2ToMeditSol(solfil, 
                                        solfil.replace(sol_ext_cfd, '.sol'))

            #--- Writing (if necessary) the parameter file for Hausdorff distances
            fileconverter.WriteParamFile(config_mmg, meshfil.rstrip('.su2'))

            #--- Adapt mesh with MMG
            meshin = config_cfd['MESH_FILENAME'].replace('.su2', medit_ext)
            meshout = config_cfd['MESH_OUT_FILENAME'].replace('.su2', medit_ext)

            solfile = config_cfd['RESTART_FILENAME'].replace(sol_ext_cfd,'.sol')
            su2adap.call_mmg(meshin, meshout, solfile, config_mmg)

            meshnewdict = fileconverter.ReadMeshMedit(meshout)

            #--- Dumping a copy of the adapted mesh 
            fileconverter.WriteMeshSU2(meshout.replace(medit_ext,'.su2'))

            #--- Print mesh sizes
            su2adap.print_adap_table(iSiz, iSub, pyadap_dict, meshnewdict)

            if runCFD:

                dir = f'./ite{global_iter}'
                os.makedirs(os.path.join('..',dir))
                os.chdir(os.path.join('..',dir))

                meshfil = config_cfd['MESH_FILENAME']
                solfil_ini  = f'restart_flow{sol_ext_cfd}'

                fileconverter.WriteMeshSU2(meshfil)

                del meshnewdict

                #--- Run su2

                # link to restart if restart file is binary
                if sol_ext_cfd != '.csv':
                    os.symlink(os.path.join(f'../ite{global_iter-1}', config.RESTART_FILENAME), solfil_ini)

                try: # run with redirected outputs

                    su2adap.update_flow_config(config_cfd, meshfil, solfil, solfil_ini, pyadap_dict, iSiz, iSub)

                    with su2io.redirect.output('su2.out'): SU2_CFD(config_cfd)

                    if not os.path.exists(solfil) :
                        raise RuntimeError('SU2_CFD failed.\n')

                    #--- Print convergence history

                    npoin = su2adap.get_su2_npoin(meshfil)
                    su2adap.plot_results(history_format, history_filename, global_iter, npoin)

                    if gol:

                        su2adap.update_adj_config(config_cfd_ad, meshfil, solfil, adjsolfil,
                                                adjsolfil_ini, pyadap_dict, iSiz, iSub)

                        with su2io.redirect.output('su2.out'): SU2_CFD(config_cfd_ad)

                        if not os.path.exists(su2io.add_suffix(adjsolfil, suffix)) :
                            raise RuntimeError('SU2_CFD_AD failed.\n')

                except:
                    raise

            del fileconverter

    #--- Write final files

    fileconverter = su2adap.MeshSolConverter()
    fileconverter.SU2ToMeditMesh(meshfil, meshfil.replace('.su2', medit_ext))

    os.rename(solfil, os.path.join(base_dir, config.RESTART_FILENAME))
    os.rename(meshfil, os.path.join(base_dir, config.MESH_OUT_FILENAME))

    pad_nul = ' '*15
    print('\nMesh adaptation successfully ended.')
    print(f'Results files: {config.MESH_OUT_FILENAME}\n{pad_nul}{config.RESTART_FILENAME}')
