#!/usr/bin/env python

## \file interface.py
#  \brief Wrapper functions for interfacing with the Inria AMG library
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

import sys, os, shutil
import subprocess
import time
import numpy as np
from struct import unpack, pack, calcsize, iter_unpack
import itertools

# ------------------------------------------------------------
#  Setup
# ------------------------------------------------------------
MMG_RUN = os.environ["MMG_RUN"]
sys.path.append(MMG_RUN)
command_mmg2D = MMG_RUN+'/mmg2d_O3'
command_mmg3D = MMG_RUN+'/mmg3d_O3'

# ------------------------------------------------------------
#  SU2-MMG Interface Functions
# ------------------------------------------------------------

class MeshSolConverter():
    """
    Class to convert .su2 to .mesh files and viceversa.
    Class to convert .csv or .dat SS2 solution files to .sol
    Works only with 2D-triangular and 3D-tetrahedral unstructured meshes.
    """

    def __init__(self, verbose=False):
        self.verbose = verbose
        return
    
    def SetDim(self, dim):
        """
        Setting the number of dimensions.
        """
        self.dim = dim

    def GetDim(self):
        """
        Returning the number of dimensions.
        """
        return self.dim
    
    def SetMeshDict(self, meshDict):
        """
        Setting the mesh dictionary once read either from .su2 or .mesh
        """
        self.meshDict = meshDict

    def GetMeshDict(self):
        """
        Returning the mesh dictionary.
        """
        return self.meshDict
    
    def SetMetricDict(self, metricDict):
        """
        Setting the dictionary for adaptation metric values.
        """
        self.metricDict = metricDict

    def GetMetricDict(self):
        """
        Returning the dictionary for adaptation metric values.
        """
        return self.metricDict
        
    def SetSU2MeditMarkersMap(self, su2MarkersList):
        """
        Constructing a unique between SU2 markers and Medit colors.
        """
        self.markersMap = []
        for meditTag, su2Tag in enumerate(su2MarkersList):
            # medit_tag starts from "1" since the tag "0" is left for the volume domain
            self.markersMap.append([str(meditTag+1), su2Tag])

    def GetSU2MeditMarkersMap(self):
        """
        Returning the map bewteen SU2 markers and Medit colors.
        """
        if self.markersMap:
            return self.markersMap
        else:
            raise ValueError('Su2-Medit markers map not set!')  

    def GetSU2Marker(self, meditTag):
        """
        Returning the SU2 markers correspondent to a Medit color.
        """
        for match in self.markersMap:
            if match[0] == meditTag:
                return match[1]

    def GetMeditMarker(self, su2Tag):
        """
        Returning the Medit color correspondent to a SU2 marker.
        """
        for match in self.markersMap:
            if match[1] == su2Tag:
                return match[0]

    def ReadMeshSU2(self, su2Filename):
        """ 
        Reads a .su2 mesh file and returns node coordinates, elements, and boundary markers in a dictionary data structure. 
        """
        meshDict = read_SU2_mesh_ascii(self, su2Filename)
        return meshDict

    def ReadSolSU2(self, su2Filename):
        """
        Reads a .csv/.dat SU2 solution file to obtain the metric. 
        """
        if '.dat' in su2Filename:
            metricDict = read_SU2_restart_binary(su2Filename)
        elif '.csv' in su2Filename:
            metricDict = read_SU2_restart_ascii(su2Filename) 
        
        self.SetMetricDict(metricDict)

        return metricDict
    
    def ReadMeshMedit(self, meditFilename):
        """ 
        Reads a .mesh/b file and returns node coordinates, elements, and boundary markers  in a dictionary data structure. 
        """
        if meditFilename.endswith('mesh'):
            meshDict = read_medit_mesh_ascii(self, meditFilename)
        elif meditFilename.endswith('meshb'):
            meshDict = read_medit_mesh_binary(self, meditFilename, verbose=self.verbose)
        return meshDict

    def WriteMeshSU2(self, su2Filename):
        """ 
        Writes a .su2 mesh file from given mesh data. 
        """
        write_su2_mesh_ascii(self, su2Filename)
        return
    
    def WriteMeshMedit(self, meditFilename):
        """ 
        Writes a .mesh/b mesh file from given mesh data. 
        """
        if meditFilename.endswith('.mesh'):
            write_medit_mesh_ascii(self, meditFilename)
        elif meditFilename.endswith('.meshb'):
            write_medit_mesh_binary(self, meditFilename, verbose=self.verbose)
        return

    def WriteSolMedit(self, meditFilename):
        """ 
        Writes a .sol Medit file from given metric data. 
        """
        metric_dict = self.GetMetricDict()

        dim = metric_dict['Dim']
        numvert = metric_dict['NumberVertices']

        header = 'MeshVersionFormatted 2\nDimension %i\nSolAtVertices\n%i\n1 3\n' % (dim , numvert)
        footer = '\nEnd\n'

        sol_data = np.empty((numvert,0))

        sol_data = np.hstack((sol_data, metric_dict['Metric_xx'][:,np.newaxis]))
        sol_data = np.hstack((sol_data, metric_dict['Metric_xy'][:,np.newaxis]))
        sol_data = np.hstack((sol_data, metric_dict['Metric_yy'][:,np.newaxis]))

        if dim == 3:
            sol_data = np.hstack((sol_data, metric_dict['Metric_xz'][:,np.newaxis]))
            sol_data = np.hstack((sol_data, metric_dict['Metric_yz'][:,np.newaxis]))
            sol_data = np.hstack((sol_data, metric_dict['Metric_zz'][:,np.newaxis]))

        np.savetxt(meditFilename, sol_data, delimiter=' ', header=header, footer=footer, comments='', fmt='%1.5e')
        
        return
    
    def SU2ToMeditMesh(self, su2Filename, meditFilename):
        """
        Full mesh file conversion (reading-writing) from SU2 to Medit
        """
        self.ReadMeshSU2(su2Filename)
        self.WriteMeshMedit(meditFilename)
        if self.verbose:
            print(f"Converted {su2Filename} to {meditFilename}")

    def SU2ToMeditSol(self, su2Filename, meditFilename):
        """
        Full sol file conversion (reading-writing) from SU2 to Medit
        """
        self.ReadSolSU2(su2Filename)
        self.WriteSolMedit(meditFilename)
        if self.verbose:
            print(f"Converted {su2Filename} to {meditFilename}")

    def MeditToSU2Mesh(self, meditFilename, su2Filename):
        """
        Full mesh file conversion (reading-writing) from Medit to SU2
        """
        self.ReadMeshMedit(meditFilename)
        self.WriteMeshSU2(su2Filename)
        if self.verbose:
            print(f"Converted {meditFilename} to {su2Filename}")

    def WriteParamFile(self, configMmg, meshFilename):
        """
        Writing the .mmg2d/.mmg3d parameter file if required. 
        """
        param_required = isinstance(configMmg['hausd'], dict)
        if param_required:
            meshDict = self.GetMeshDict()
            dim = meshDict['Dim']
            if dim == 2:
                boundaries = meshDict['Edges']
                elem_type = 'Edges'
                mmg_ext = '.mmg2d'
            if dim == 3:
                boundaries = meshDict['Triangles']
                elem_type = 'Triangles'
                mmg_ext = '.mmg3d'

            param_filename = meshFilename + mmg_ext
            with open(param_filename, 'w') as f:
                f.write('Parameters\n')
                f.write(str(len(boundaries.keys()))+'\n')
                f.write('\n')

                if len(boundaries.keys()) != len(configMmg['hausd'].keys()):
                    print('WARNING: Different number of markers between SU2 (%i) mesh and MMG parameters (%i). ' \
                                  'For unspecified markers, HAUSD = 0.01 is assumed.' %
                                   (len(boundaries.keys()), len(configMmg['hausd'].keys())))
                
                for su2_tag in configMmg['hausd'].keys():
                    f.write('%s %s %1.2e %1.2e %1.2e\n' % 
                            (self.GetMeditMarker(su2_tag), 
                            elem_type, 
                            configMmg['hmin'], 
                            configMmg['hmax'], 
                            configMmg['hausd'][su2_tag]))       
        
        else:
            pass

        return      

def call_mmg(meshin, meshout, solfile, config_mmg):
    """Adapt mesh using pyamg module"""

    remesh_options = config_mmg
    remesh_options['meshfile'] = meshin
    remesh_options['meshoutfile'] = meshout
    remesh_options['solfile'] = solfile

    dim = config_mmg['dim']

    if dim == 2:
        the_Command = build_command(command_mmg2D, remesh_options)
    elif dim == 3:
        the_Command = build_command(command_mmg3D, remesh_options)
    else:
        raise KeyError('Wrong number of dimensions!')
    
    try:
        run_command(the_Command)
    except:
        raise RuntimeError("mmg failed.")

    return 

# ------------------------------------------------------------
#  Helper functions
# ------------------------------------------------------------

def build_command(command_mmg, options):
    """builds the mmg command with options"""

    the_Command = command_mmg 
    the_Command += ' -in '    + options['meshfile'] 
    the_Command += ' -sol '   + options['solfile'] 
    the_Command += ' -out '   + options['meshoutfile']
    the_Command += ' -hmin '  + str(options['hmin'])
    the_Command += ' -hmax '  + str(options['hmax'])
    if 'hgrad' in options.keys():
        the_Command += ' -hgrad ' + str(options['hgrad'])
    if 'hausd' in options.keys() and not isinstance(options['hausd'], dict):
        the_Command += ' -hausd ' + str(options['hausd'])
    if 'ar' in options.keys():
        the_Command += ' -ar ' + str(options['ar'])
    the_Command += ' > ' + options['mmg_log'] + ' 2> ' + options['mmg_err']
    
    return the_Command


def run_command(Command):
    """runs os command with subprocess
    checks for errors from command
    """

    sys.stdout.flush()

    proc = subprocess.Popen(
        Command, shell=True, stdout=sys.stdout, stderr=subprocess.PIPE
    )
    return_code = proc.wait()
    message = proc.stderr.read().decode()

    return 

tclock = dict()

def tic(ref=0):
    global tclock
    tclock[ref] = time.time()


def toc(ref=0):
    global tclock
    return format(time.time()-tclock[ref], "0.2f")+"s"

# ------------------------------------------------------------
#  Mesh reading
# ------------------------------------------------------------

def read_SU2_mesh_ascii(mesh, meshFilename):
    """ 
    Reads a .su2 mesh file and returns node coordinates, elements, and boundary markers in a dictionary data structure. 
    """
    with open(meshFilename, "r") as f:
        lines = f.readlines()

    vertices = []
    elements = []
    boundaries = {}

    su2_markers_list = []
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("%"):
            i += 1
            continue

        if line.startswith("NDIME="):
            dim = int(line.split("=")[1].strip())
            mesh.SetDim(dim)

        elif line.startswith("NPOIN="):
            n_vertices = int(line.split("=")[1].strip())
            if dim == 2:
                for j in range(n_vertices):
                    x, y = lines[i + 1 + j].split()[:dim]
                    vertices.append([float(x), float(y)])
            elif dim == 3:
                for j in range(n_vertices):
                    x, y, z = lines[i + 1 + j].split()[:dim]
                    vertices.append([float(x), float(y), float(z)])
            i += n_vertices  # Move index past nodes

        elif line.startswith("NELEM="):
            n_elements = int(line.split("=")[1].strip())
            for j in range(n_elements):
                elem_type, *id_vert = lines[i + 1 + j].split()[:(dim+2)]
                elements.append([int(vert) for vert in id_vert])
            i += n_elements  # Move index past elements

        elif line.startswith("NMARK="):
            n_markers = int(line.split("=")[1].strip())
            for _ in range(n_markers):
                i += 1
                marker_tag = lines[i].split("=")[1].strip()
                su2_markers_list.append(marker_tag)
                i += 1
                n_faces = int(lines[i].split("=")[1].strip())
                boundaries[marker_tag] = []
                for j in range(n_faces):
                    face_type, *id_vert = lines[i + 1 + j].split()
                    boundaries[marker_tag].append([int(vert) for vert in id_vert])
                i += n_faces  # Move index past boundary elements

        i += 1

    if int(elem_type) != 5 and int(elem_type) != 10:
        raise Exception('The .su2 mesh file containes volume elements different from triangles/tetrahedra')

    meshDict = {'Dim': dim, 'Vertices': vertices}
    if int(elem_type) == 5:
        meshDict['Triangles'] = elements
    if int(elem_type) == 10:
        meshDict['Tetrahedra'] = elements

    # reordering markers in alphabetical orders
    #boundaries = {key: value for key, value in sorted(boundaries.items())}
    if int(face_type) == 3:
        meshDict['Edges'] =  boundaries
    if int(face_type) == 5:
        meshDict['Triangles'] =  boundaries

    mesh.SetMeshDict(meshDict)
    mesh.SetSU2MeditMarkersMap(su2_markers_list)

    return meshDict


def read_medit_mesh_ascii(mesh, meshFilename):

    with open(meshFilename, "r") as f:
        lines = f.readlines()

    dim = None
    vertices = []
    elements = []
    boundaries = {}

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("Dimension"):
            dim = int(line.split()[1])
            mesh.SetDim(dim)

        elif line.startswith("Vertices"):
            print(line)
            print(lines[i+1])
            nvert = int(lines[i+1])
            if dim == 2:
                for j in range(nvert):
                    x, y, id_dom = lines[i + 2 + j].split()
                    vertices.append([float(x), float(y)])
            elif dim == 3:
                for j in range(nvert):
                    x, y, z, id_dom = lines[i + 2 + j].split()
                    vertices.append([float(x), float(y), float(z)])
                
            i += nvert  # Move index past nodes

        elif (line.startswith("Triangles") and dim == 2) or (line.startswith("Tetrahedra") and dim == 3):
            nelem = int(lines[i+1])
            elements = [None]*nelem
            for j in range(nelem):
                elem = list(map(int, lines[i + 2 + j].split()))
                elem = [el-1 for el in elem]
                elements[j] = elem[:-1]  # Last column is a region marker
            i += nelem  # Move index past elements             

        elif (line.startswith("Edges") and dim == 2) or (line.startswith("Triangles") and dim == 3):
            # These define boundary markers
            nface = int(lines[i+1])
            for j in range(nface):
                face = list(map(int, lines[i + 2 + j].split()))
                marker = str(face[-1])  # Last column is the boundary marker
                face = [fa - 1 for fa in face[:-1]]
                if marker not in boundaries.keys():
                    boundaries[marker] = []
                boundaries[marker].append(face)  # Store only connectivity
            i += nface  # Move index past boundary elements

        i += 1

    meshDict = {'Dim': dim, 'Vertices': vertices}
    if dim == 2:
        meshDict['Triangles'] = elements
    if dim == 3:
        meshDict['Tetrahedra'] = elements

    # reordering markers in alphabetical orders
    if dim == 2:
        meshDict['Edges'] =  boundaries
    if dim == 3:
        meshDict['Triangles'] =  boundaries

    meshDict = check_surplus_points(meshDict)

    mesh.SetMeshDict(meshDict)

    return meshDict


def read_medit_mesh_binary(mesh, meshFilename, verbose=False):
    """Reload the mesh with meshFile (binary format).

    INPUT:
        meshFile :  the path of the .meshb file (binary format)
    """

    vertices = []
    elements = []
    boundaries = {}

    if not meshFilename.endswith(".meshb"):
        meshFilename += ".meshb"

    tic()
    f = open(meshFilename, "rb")

    code = read_int(f)
    if not code:
        raise Exception("Error in reading the binary file " + meshFilename)
    
    meshVersionFormatted = read_int(f)
    if not meshVersionFormatted in [1, 2]:
        raise Exception("Wrong MeshVersionFormatted. Should be 1 or 2.")
    
    if gmfKwdCod[read_int(f)] != 'GmfDimension':
        raise Exception("Error in reading the binary file " + meshFilename + ". "
                        "GmfDimension expected.")
    
    nextPos = read_int(f)
    if verbose: print(f"Next field position: {nextPos}.")
    dim = read_int(f)
    if verbose: print("Dimension " + str(dim))
    if not dim in [2, 3]:
        raise Exception("Error in reading the binary file " + meshFilename +
                        ". The dimension should be 2 or 3.")
    else:
        mesh.SetDim(dim)

    def readField(f, nfields, title):
        nextPos = read_int(f)
        if verbose: print(f"Next field position: {nextPos}.")
        n = read_int(f)
        if verbose: print(f"{n} {title}.")
        code = "i" * nfields * n
        nbytes = 4 * nfields * n
        field = np.asarray(unpack(code, f.read(nbytes)))
        return (n, field.reshape((n, nfields)))

    while True:
        try:
            kwdCod = read_int(f)
        except:
            if verbose: print("Warning: end of file without the END keyword.")
            break
        if verbose: print("Reading "+gmfKwdCod[kwdCod])

        # Vertices
        if gmfKwdCod[kwdCod] == 'GmfVertices':

            nextPos = read_int(f)
            if verbose: print(f"Next field position: {nextPos}.")

            nvert = read_int(f)
            if verbose: print(f"{nvert} Vertices.")

            if meshVersionFormatted == 1:
                # Float precision
                code = "="+("f"*dim + "i") * nvert
                nbytes = calcsize(code)
            elif meshVersionFormatted == 2:
                # Double precision
                code = "="+("d"*dim + "i") * nvert
                nbytes = calcsize(code)

            vertices = np.asarray(unpack(code, f.read(nbytes)))
            vertices = vertices.reshape((nvert, dim+1)).tolist()
            vertices = [vert[:-1] for vert in vertices]

        # Edges
        elif gmfKwdCod[kwdCod] == 'GmfEdges':
            nedge, edges = readField(f, 3, gmfKwdCod[kwdCod][3:])
            edges.tolist()

            if dim == 2:
                for edg in edges:
                    if str(edg[-1]) not in boundaries.keys():
                        boundaries[str(edg[-1])] = []
                    boundaries[str(edg[-1])].append([ed - 1 for ed in edg[:-1]])

        # Triangles
        elif gmfKwdCod[kwdCod] == 'GmfTriangles':
            ntria, triangles = readField(f, 4, gmfKwdCod[kwdCod][3:])
            triangles.tolist()

            if dim == 2: # triangles as elements
                elements = [[tr-1 for tr in tria[:-1]] for tria in triangles]
            elif dim == 3: # triangles as boundaries
                for tria in triangles:
                    if str(tria[-1]) not in boundaries.keys():
                        boundaries[str(tria[-1])] = []
                    boundaries[str(tria[-1])].append([tr-1 for tr in tria[:-1]])

        # Tetrahedra
        elif gmfKwdCod[kwdCod] == 'GmfTetrahedra':
            ntetra, tetrahedra = readField(f, 5, gmfKwdCod[kwdCod][3:])
            tetrahedra.tolist()
            elements = [[te-1 for te in tet[:-1]] for tet in tetrahedra]

        # Fields that need to be read but not saved
        elif (gmfKwdCod[kwdCod] == 'GmfCorners'          or 
              gmfKwdCod[kwdCod] == 'GmfRequiredVertices' or
              gmfKwdCod[kwdCod] == 'GmfRequiredEdges'    or
              gmfKwdCod[kwdCod] == 'GmfRidges'           or
              gmfKwdCod[kwdCod] == 'GmfRequiredTriangles'):
            _, _ = readField(f, 1, gmfKwdCod[kwdCod][3:])

        elif (gmfKwdCod[kwdCod] == 'GmfNormalAtVertices' or 
              gmfKwdCod[kwdCod] == 'GmfTangentAtVertices'):
            _, _ = readField(f, 2, gmfKwdCod[kwdCod][3:])


        elif (gmfKwdCod[kwdCod] == 'GmfNormals' or
              gmfKwdCod[kwdCod] == 'GmfTangents'):
            nextPos = read_int(f)
            if verbose: print(f"Next field position: {nextPos}.")
            nvn = read_int(f)
            if verbose: print(f"{nvn} Normals.")
            if meshVersionFormatted == 1:
                code = "f"*(3*nvn)
                nbytes = calcsize(code)
            else:
                code = "d"*(3*nvn)
                nbytes = calcsize(code)
            _ = np.asarray(unpack(code, f.read(nbytes)))

        # End
        elif gmfKwdCod[kwdCod] == 'GmfEnd':
            if verbose: print("End of mesh.")
            break

        else:
            raise KeyError("Error, field "+gmfKwdCod[kwdCod]+" not supported.")
        
    f.close()
    if verbose: print("Read " + meshFilename + " in "+toc()+".")

    meshDict = {'Dim': dim, 'Vertices': vertices}
    if dim == 2:
        meshDict['Triangles'] = elements
    if dim == 3:
        meshDict['Tetrahedra'] = elements

    # reordering markers in alphabetical orders
    if dim == 2:
        meshDict['Edges'] =  boundaries
    if dim == 3:
        meshDict['Triangles'] =  boundaries

    meshDict = check_surplus_points(meshDict)

    mesh.SetMeshDict(meshDict)
    
    return meshDict


def check_surplus_points(meshDict):
    """
    Routine to check and eliminate for surplus points in medit mesh generation
    """

    dim = meshDict['Dim']

    vertices = meshDict['Vertices']
    nvert = len(vertices)

    if dim == 2:
        elements   = meshDict['Triangles']
        boundaries = meshDict['Edges']
    elif dim == 3:
        elements   = meshDict['Tetrahedra']
        boundaries = meshDict['Triangles']

    print("Checking for surplus points")

    PointIDs = np.array(range(nvert), dtype=int)
    elements = np.array(elements, dtype=int)
    isIn = np.isin(PointIDs, elements)
    whichAreSurplusPoints = np.where(isIn == False)[0]
    
    if (len(whichAreSurplusPoints) > 0):
        print("There are surplus points in mesh file for unknown reasons.")
        print("Surplus points:", whichAreSurplusPoints)
        print("Deleting it and then fix connectivity...")
        SubtractIDs = np.zeros((nvert, ), dtype=int)
        for iPoint in whichAreSurplusPoints:
            IDsOver = np.where(PointIDs > iPoint)[0]
            SubtractIDs[IDsOver] += 1
        
        # Now I can just fix the connectivity
        elements -= SubtractIDs[elements]
        for marker in boundaries.keys():
            markerData = np.array(boundaries[marker], dtype=int)
            markerData -= SubtractIDs[markerData]
            boundaries[marker] = markerData

        whichAreSurplusPointsReverse = np.sort(whichAreSurplusPoints)[::-1]
        for iPoint in whichAreSurplusPointsReverse:
            # Remove from the list of vertices but in reverse order
            vertices.pop(iPoint)

    meshDict['Vertices'] = vertices
    if dim == 2:
        meshDict['Triangles'] = elements
    if dim == 3:
        meshDict['Tetrahedra'] = elements

    # reordering markers in alphabetical orders
    if dim == 2:
        meshDict['Edges'] = boundaries
    if dim == 3:
        meshDict['Triangles'] =  boundaries

    return meshDict


# ------------------------------------------------------------
#  Mesh writing
# ------------------------------------------------------------

def write_su2_mesh_ascii(mesh, meshFilename):
    """ 
    Writes a .su2 mesh file from given mesh data. 
    """
    meshDict = mesh.GetMeshDict()
    dim = meshDict['Dim']
    vertices = meshDict["Vertices"]
    if dim == 2:
        elements = meshDict['Triangles']
        elem_type = 5
        boundaries = meshDict['Edges']
        face_type = 3
    if dim == 3:
        elements = meshDict['Tetrahedra']
        elem_type = 10
        boundaries = meshDict['Triangles']
        face_type = 5


    with open(meshFilename, "w") as f:
        f.write("NDIME= {}\n".format(dim))

        f.write("NELEM= {}\n".format(len(elements)))
        for i, elem in enumerate(elements):
            f.write("{} {} {}\n".format(elem_type, " ".join(map(str, elem)), i))

        f.write("NPOIN= {}\n".format(len(vertices)))
        for i, node in enumerate(vertices):
            f.write("{} {}\n".format(" ".join(map(str, node)), i))

        f.write("NMARK= {}\n".format(len(boundaries)))
        for medit_tag in boundaries.keys():
            f.write("MARKER_TAG= {}\n".format(mesh.GetSU2Marker(medit_tag)))
            f.write("MARKER_ELEMS= {}\n".format(len(boundaries[medit_tag])))
            for face in boundaries[medit_tag]:
                f.write("{} {}\n".format(face_type, " ".join(map(str, face))))

    return


def write_medit_mesh_ascii(mesh, meshFilename):
    """ 
    Writes a .mesh mesh file from given mesh data. 
    """
    meshDict = mesh.GetMeshDict()
    dim = meshDict['Dim']
    vertices = meshDict['Vertices']
    if dim == 2:
        elements = meshDict['Triangles']
        boundaries = meshDict['Edges']
    if dim == 3:
        elements = meshDict['Tetrahedra']
        boundaries = meshDict['Triangles']

    with open(meshFilename, "w") as f:
        f.write("MeshVersionFormatted 2\n")
        f.write("Dimension {}\n".format(meshDict['Dim']))
        
        # Write nodes
        f.write("\nVertices \n{}\n".format(len(meshDict['Vertices'])))
        for vert in vertices:
            f.write(" ".join(map(str, vert[:dim])) + " 0\n")  # 0 is the default region ID
        
        # Write elements (assume triangles for 2D, tetrahedra for 3D)
        if dim == 2:
            f.write("\nTriangles \n{}\n".format(len(elements)))
        elif dim == 3:
            f.write("\nTetrahedra \n{}\n".format(len(elements)))
        
        for elem in elements:
            elem = [el + 1 for el in elem]
            f.write(" ".join(map(str, elem)) + " 0\n")  # Last value is a region ID
        
        # Write boundary elements correctly
        if boundaries:
            if dim == 2:
                f.write("\nEdges \n{}\n".format(sum(len(faces) for faces in boundaries.values())))
            elif dim == 3:
                f.write("\nTriangles \n{}\n".format(sum(len(faces) for faces in boundaries.values())))
            
            for su2_tag in boundaries.keys():
                for face in boundaries[su2_tag]:
                    face = [fa + 1 for fa in face]
                    f.write(" ".join(map(str, face)) + " {}\n".format(mesh.GetMeditMarker(su2_tag)))

        f.write("\nEnd\n")
    
    return
    

def write_medit_mesh_binary(mesh, meshFilename, verbose=False):
    """Save a mesh in the INRIA binary file format"""
    
    if not meshFilename.endswith(".meshb"):
        meshFilename += ".meshb"

    meshDict = mesh.GetMeshDict()
    dim = meshDict['Dim']
    vertices = meshDict['Vertices']
    if dim == 2:
        elements = meshDict['Triangles']
        boundaries = meshDict['Edges']
    if dim == 3:
        elements = meshDict['Tetrahedra']
        boundaries = meshDict['Triangles']

    tic()
    f = open(meshFilename, "wb")
    f.write(pack("i", 1))  # Write code
    f.write(pack("i", 2))  # Write MeshVersionFormatted 2

    f.write(pack("i", indicesGmf['GmfDimension']))
    f.write(pack("i", 20))  # NextPos
    if not dim in [2, 3]:
        raise Exception("Error, the mesh dimension is not 2 or 3")
    f.write(pack("i", dim))

    # Vertices
    if verbose: print(f"Write {len(vertices)} Vertices.")
    f.write(pack("i", indicesGmf['GmfVertices']))
    code = "=" + ("d" * dim + "i") * len(vertices)
    nextPos = f.tell() + calcsize("ii") + calcsize(code)
    if verbose: print(f"Next position: {nextPos}")
    f.write(pack("i", nextPos))
    vertices = [vert+[0] for vert in vertices] # appending the region ID
    f.write(pack("i", len(vertices)))
    data = list(itertools.chain.from_iterable(vertices))
    f.write(pack(code, *data))

    # Edges
    if dim == 2:
        nedge = sum(len(faces) for faces in boundaries.values())
        if verbose: print(f"Write {nedge} "+"Edges.")
        f.write(pack("i", indicesGmf['GmfEdges']))
        code = "i" * 3 * nedge

        nextPos = f.tell()+calcsize(code)+calcsize("ii")
        if verbose: print(f"Next position: {nextPos}")
        f.write(pack("i", nextPos))  # NulPos
        f.write(pack("i", nedge))
        edges = []
        for su2_tag in boundaries.keys():
            for edg in boundaries[su2_tag]:
                edges.append([ed+1 for ed in edg] + [int(mesh.GetMeditMarker(su2_tag))])
        data = list(itertools.chain.from_iterable(edges))
        f.write(pack(code, *data))

    # Triangles
    if dim == 2: # triangles as elements
        ntria = len(elements)
    elif dim == 3: # triangles as boundaries
        ntria = sum(len(faces) for faces in boundaries.values())

    if verbose: print(f"Write {ntria} "+"Triangles.")
    f.write(pack("i", indicesGmf['GmfTriangles']))
    code = "i" * 4 * ntria

    nextPos = f.tell()+calcsize(code)+calcsize("ii")
    if verbose: print(f"Next position: {nextPos}")
    f.write(pack("i", nextPos))  # NulPos
    f.write(pack("i", ntria))
    triangles = []
    if dim == 2:
        for elem in elements:
            triangles.append([el + 1 for el in elem] + [0])
    if dim == 3:
        for su2_tag in boundaries.keys():
            for tria in boundaries[su2_tag]:
                triangles.append([tr + 1 for tr in tria] +[int(mesh.GetMeditMarker(su2_tag))])   

    data = list(itertools.chain.from_iterable(triangles))
    f.write(pack(code, *data))

    # Tetrahedra
    if dim == 3:
        ntetra = len(elements)

        if verbose: print(f"Write {ntetra} "+"Tetrahedra.")
        f.write(pack("i", indicesGmf['GmfTetrahedra']))
        code = "i" * 5 * ntetra

        nextPos = f.tell()+calcsize(code)+calcsize("ii")
        if verbose: print(f"Next position: {nextPos}")
        f.write(pack("i", nextPos))  # NulPos
        f.write(pack("i", ntetra))

        tetrahedra = []
        for elem in elements:
            tetrahedra.append([el + 1 for el in elem] + [0])

        data = list(itertools.chain.from_iterable(tetrahedra))
        f.write(pack(code, *data))

    # End
    f.write(pack("i", indicesGmf['GmfEnd']))
    nextPos = f.tell()+calcsize("i")
    # Final size
    f.write(pack("i",nextPos))
    f.close()

    if verbose: print("Wrote "+meshFilename+" in "+toc()+".")

    return


# ------------------------------------------------------------
#  Restart reading
# ------------------------------------------------------------

CGNS_STRING_SIZE = 33  # Fixed string size per CGNS standard

def read_SU2_restart_binary(filename):
    """
    Read SU2 binary restart file and return fields and data array.

    Returns:
        fields (List[str]): Field names including "Point_ID".
        data (np.ndarray): Data array of shape (nPoints, nFields-1).

    Note that the Point_ID column is implicit in the ordering
    """

    restartFields = []  

    with open(filename, 'rb') as f:
        # Read 5 integers (magic number + metadata)
        header = np.fromfile(f, dtype=np.int32, count=5)
        if header.size != 5:
            raise RuntimeError("Error reading header from restart file.")
        
        magic_number, nFields, nPoints, _, _ = header

        # Check the magic number
        if magic_number != 535532:
            raise RuntimeError(f"{filename} is not a binary SU2 restart file.")

        # Read field names (each is CGNS_STRING_SIZE characters)
        for _ in range(nFields):
            name_bytes = f.read(CGNS_STRING_SIZE)
            name_str = name_bytes.decode('utf-8').strip('\x00').strip()
            restartFields.append(name_str)

        # Read restart data as a flat array of doubles
        data = np.fromfile(f, dtype=np.float64, count=nFields * nPoints)

        if data.size != nFields * nPoints:
            raise RuntimeError("Error reading restart data.")

        # Reshape to 2D: each row is a point, each column is a field
        data = data.reshape((nPoints, nFields))

    metric_dict = {'NumberVertices': nPoints}

    if 'z' in restartFields:
        metric_dict['Dim'] = 3
        fieldsToRead = ['Metric_xx', 'Metric_xy', 'Metric_yy', 'Metric_xz', 'Metric_yz', 'Metric_zz']
    else:
        metric_dict['Dim'] = 2  
        fieldsToRead = ['Metric_xx', 'Metric_xy', 'Metric_yy']

    for field in fieldsToRead:

        try:
            ind_metric = restartFields.index(field)
            metric_dict[field] = data[:, ind_metric]
        except ValueError:
            print('The metric field %s is missing!' % field)
            exit()

    return metric_dict


def read_SU2_restart_ascii(filename):
    """
    Read SU2 ASCII restart file and return fields and data array.

    Returns:
        fields (List[str]): Field names.
        data (np.ndarray): Data array of shape (nPoints, nFields).
    """

    # reading the first line to get the fields name
    try:
        with open(filename, "r") as f:
            line = f.readline()
    except:
        raise("The solution file must be in ASCII format!")

    restartFields = line.lstrip('"').rstrip('"\n')
    restartFields = restartFields.split('","')

    data = np.genfromtxt(filename, delimiter=',', skip_header=1, dtype=np.float64)
    nPoints = data.shape[0]

    metric_dict = {'NumberVertices': nPoints}

    if 'z' in restartFields:
        metric_dict['Dim'] = 3
        fieldsToRead = ['Metric_xx', 'Metric_xy', 'Metric_yy', 'Metric_xz', 'Metric_yz', 'Metric_zz']
    else:
        metric_dict['Dim'] = 2  
        fieldsToRead = ['Metric_xx', 'Metric_xy', 'Metric_yy']

    for field in fieldsToRead:

        try:
            ind_metric = restartFields.index(field)
            metric_dict[field] = data[:, ind_metric]
        except ValueError:
            print('The metric field %s is missing!' % field)
            exit()

    return metric_dict


# ------------------------------------------------------------
#  Other utilities for I/O
# ------------------------------------------------------------
def read_int(f):
    return unpack("i", f.read(4))[0]


def next_line(f):
    while True:
        line = f.readline().strip()
        if line:
            break
    return line


def readField(f, nfields, title, verbose):
    nextPos = read_int(f)
    print(f"Next field position: {nextPos}.")
    n = read_int(f)
    print(f"{n} {title}.")
    code = "i"*nfields*n
    nbytes = 4*nfields*n
    field = np.asarray(unpack(code, f.read(nbytes)))
    return (n, field.reshape((n, nfields)))


gmfKwdCod = ['GmfReserved1',
            'GmfVersionFormatted',
            'GmfReserved2',
            'GmfDimension',
            'GmfVertices',
            'GmfEdges',
            'GmfTriangles',
            'GmfQuadrilaterals',
            'GmfTetrahedra',
            'GmfPentahedra',
            'GmfHexahedra',
            'GmfReserved3',
            'GmfReserved4',
            'GmfCorners',
            'GmfRidges',
            'GmfRequiredVertices',
            'GmfRequiredEdges',
            'GmfRequiredTriangles',
            'GmfRequiredQuadrilaterals',
            'GmfTangentAtEdgeVertices',
            'GmfNormalAtVertices',
            'GmfNormalAtTriangleVertices',
            'GmfNormalAtQuadrilateralVertices',
            'GmfAngleOfCornerBound',
            'GmfReserved5',
            'GmfReserved6',
            'GmfReserved7',
            'GmfReserved8',
            'GmfReserved9',
            'GmfReserved10',
            'GmfReserved11',
            'GmfReserved12',
            'GmfReserved13',
            'GmfReserved14',
            'GmfReserved15',
            'GmfReserved16',
            'GmfReserved17',
            'GmfReserved18',
            'GmfReserved19',
            'GmfReserved20',
            'GmfReserved21',
            'GmfReserved22',
            'GmfReserved23',
            'GmfReserved24',
            'GmfReserved25',
            'GmfReserved26',
            'GmfReserved27',
            'GmfReserved28',
            'GmfReserved29',
            'GmfReserved30',
            'GmfBoundingBox',
            'GmfReserved31',
            'GmfReserved32',
            'GmfReserved33',
            'GmfEnd',
            'GmfReserved34',
            'GmfReserved35',
            'GmfReserved36',
            'GmfReserved37',
            'GmfTangents',
            'GmfNormals',
            'GmfTangentAtVertices',
            'GmfSolAtVertices',
            'GmfSolAtEdges',
            'GmfSolAtTriangles',
            'GmfSolAtQuadrilaterals',
            'GmfSolAtTetrahedra',
            'GmfSolAtPentahedra',
            'GmfSolAtHexahedra',
            'GmfDSolAtVertices',
            'GmfISolAtVertices',
            'GmfISolAtEdges',
            'GmfISolAtTriangles',
            'GmfISolAtQuadrilaterals',
            'GmfISolAtTetrahedra',
            'GmfISolAtPentahedra',
            'GmfISolAtHexahedra',
            'GmfIterations',
            'GmfTime',
            'GmfReserved38']

indicesGmf = dict([(value, i) for i, value in enumerate(gmfKwdCod)])
