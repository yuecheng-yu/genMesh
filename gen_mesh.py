import meshio 
import numpy as np
import vtk
import pyvista as pv
import os

from collections import defaultdict
from vtk.util.numpy_support import numpy_to_vtk as n2v
from vtk.util.numpy_support import vtk_to_numpy as v2n

from vtk_functions import (
    read_geo,
    write_geo,
    extract_surface,
    threshold,
    clean,
)

def main():

    path = "/Users/yuechengyu/Work/Cardiac/svFSIplus/example/mesh/unitCube"
    # path = "./mesh"

    # create cube vol mesh
    # Define points for a unit cube
    points = np.array([
        [0.0, 0.0, 0.0],  # Bottom surface points
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],  # Top surface points
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0]
    ])
    points_data = {
        "GlobalNodeID" : np.arange(len(points)) + 1,
        "Z0" : np.where(points[:, 2] == 0, 1, 0),
        "Z1" : np.where(points[:, 2] == 1, 1, 0),
    }

    # Define cells using the VTK format for hexahedrons
    cells = [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))]
    cells_data = {
        "GlobalElementID": np.expand_dims(np.arange(len(cells)) + 1, axis=1)
    }

    # Create the volume mesh
    
    volume_mesh = meshio.Mesh(points, cells, point_data=points_data, cell_data=cells_data)

    # # Define the bottom face of the cube
    # bot_face = np.where(points[:, 2] == 0, 1, 0)
    # print(bot_face)
    # # Define the top face of the cube
    # top_face = np.where(points[:, 2] == 1, 1, 0)
    # print(top_face)

    # # Add to the point data
    # volume_mesh.point_data["bot_face"] = bot_face
    # volume_mesh.point_data["top_face"] = top_face

    # Save the volume mesh
    volume_mesh_path = path+"/unit_cube_volume.vtu"
    meshio.write(volume_mesh_path, volume_mesh, file_format="vtu")

    # # check arrays of vol mesh 
    # reader = vtk.vtkXMLUnstructuredGridReader()
    # reader.SetFileName(volume_mesh_path)
    # reader.Update()
    # vol = reader.GetOutput()
    # for i in range(vol.GetPointData().GetNumberOfArrays()):
    #     print(vol.GetPointData().GetArrayName(i))
    #     array = v2n(vol.GetPointData().GetArray(i))
    #     print(array)


    # read volume mesh in vtk
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(volume_mesh_path)
    reader.Update()
    vol = reader.GetOutput()

    # Is it necessary for unit cube??
    # reset global ids
    n_array = n2v(np.arange(vol.GetNumberOfPoints()).astype(np.int32) + 1)
    e_array = n2v(np.arange(vol.GetNumberOfCells()).astype(np.int32) + 1)
    n_array.SetName("GlobalNodeID")
    e_array.SetName("GlobalElementID")
    vol.GetPointData().AddArray(n_array)
    vol.GetCellData().AddArray(e_array)

    # map point data to cell data
    p2c = vtk.vtkPointDataToCellData()
    p2c.SetInputData(vol)
    p2c.PassPointDataOn()
    p2c.Update()
    vol = p2c.GetOutput()

    # extract surfaces
    extract = vtk.vtkGeometryFilter()
    extract.SetInputData(vol)
    # extract.SetNonlinearSubdivisionLevel(0)
    extract.Update()
    surfaces = extract.GetOutput()

    # threshold surfaces
    for name in ["Z0", "Z1"]:

        # select only current surface
        thresh = vtk.vtkThreshold()
        thresh.SetInputData(surfaces)
        thresh.SetInputArrayToProcess(0, 0, 0, 0, name)
        thresh.SetUpperThreshold(1)
        thresh.SetLowerThreshold(1)
        thresh.Update()
        surf = thresh.GetOutput()
        if surf.GetNumberOfPoints() > 0:
            surf = clean(extract_surface(surf))

        write_geo(path+"/mesh_surfaces/"+name+'.vtp', extract_surface(surf))

        # # get new GlobalNodeIDs of surface points
        # surf_ids[f + "_" + name] = v2n(
        #     surf.GetPointData().GetArray("GlobalNodeID")
        # ).tolist()

        



    # def write_geo(fname, input):
    #     """
    #     Write geometry to file
    #     Args:
    #         fname: file name
    #     """
    #     _, ext = os.path.splitext(fname)
    #     if ext == '.vtp':
    #         writer = vtk.vtkXMLPolyDataWriter()
    #     elif ext == '.vtu':
    #         writer = vtk.vtkXMLUnstructuredGridWriter()
    #     else:
    #         raise ValueError('File extension ' + ext + ' unknown.')
    #     writer.SetFileName(fname)
    #     writer.SetInputData(input)
    #     writer.Update()
    #     writer.Write()

    # # Define the eight vertices of a unit cube
    # points = vtk.vtkPoints()
    # points.InsertNextPoint(0.0, 0.0, 0.0)
    # points.InsertNextPoint(1.0, 0.0, 0.0)
    # points.InsertNextPoint(1.0, 1.0, 0.0)
    # points.InsertNextPoint(0.0, 1.0, 0.0)
    # points.InsertNextPoint(0.0, 0.0, 1.0)
    # points.InsertNextPoint(1.0, 0.0, 1.0)
    # points.InsertNextPoint(1.0, 1.0, 1.0)
    # points.InsertNextPoint(0.0, 1.0, 1.0)

    # for i in range(1):
    #     # Create a quad on the four points
    #     quad = vtk.vtkQuad()
    #     quad.GetPointIds().SetId(i*4+0, i*4+0)
    #     quad.GetPointIds().SetId(i*4+1, i*4+1)
    #     quad.GetPointIds().SetId(i*4+2, i*4+2)
    #     quad.GetPointIds().SetId(i*4+3, i*4+3)

    #     # Create a cell array to store the quad in
    #     quads = vtk.vtkCellArray()
    #     quads.InsertNextCell(quad)

    #     # Create a polydata to store everything in
    #     polydata = vtk.vtkPolyData()

    #     # Add the points and quads to the dataset
    #     polydata.SetPoints(points)
    #     polydata.SetPolys(quads)

    #     # Write to VTP file
    #     writer = vtk.vtkXMLPolyDataWriter()
    #     writer.SetFileName(f'Z{i}.vtp')
    #     writer.SetInputData(polydata)
    #     writer.Update()
    #     writer.Write()


    # # Define the six faces of the cube using the vertex points
    # # faces = vtk.vtkCellArray()
    # face_ids = [
    #     # [0, 3, 2, 1],  # bottom face
    #     [0, 1, 2, 3],  # bottom face
    #     [4, 5, 6, 7],  # top face
    #     [0, 1, 5, 4],  # front face
    #     [2, 3, 7, 6],  # back face
    #     [0, 4, 7, 3],  # left face
    #     [1, 2, 6, 5]   # right face
    # ]
    # # for id_list in face_ids:
    # #     faces.InsertNextCell(4, id_list)

    # for i in range(2):
    #     face = vtk.vtkCellArray()
    #     face.InsertNextCell(4, face_ids[i])

    #     # Create a polydata object for the surface mesh
    #     surface_mesh = vtk.vtkPolyData()
    #     surface_mesh.SetPoints(points)
    #     surface_mesh.SetPolys(face)

    #     # Write the surface mesh as a VTP file
    #     writer_surface = vtk.vtkXMLPolyDataWriter()
    #     writer_surface.SetFileName(f'unit_cube_Z{i}.vtp')
    #     writer_surface.SetInputData(surface_mesh)
    #     # writer_surface.SetDataModeToBinary()
    #     # writer_surface.SetEncodeAppendedData(False)
    #     writer_surface.Update()
    #     writer_surface.Write()

    # # Define the hexahedral cell for the volume mesh
    # hexahedron = vtk.vtkHexahedron()
    # for i in range(8):
    #     hexahedron.GetPointIds().SetId(i, i)

    # # Create an unstructured grid object for the volume mesh
    # volume_mesh = vtk.vtkUnstructuredGrid()
    # volume_mesh.SetPoints(points)
    # volume_mesh.InsertNextCell(hexahedron.GetCellType(), hexahedron.GetPointIds())

    # # Write the volume mesh as a VTU file
    # writer_volume = vtk.vtkXMLUnstructuredGridWriter()
    # writer_volume.SetFileName('unit_cube_volume.vtu')
    # writer_volume.SetInputData(volume_mesh)
    # writer_volume.Write()











    # # *** Doesn't work with svFSI
    # plane_z0 = pv.Plane((0.5, 0.5, 0.0), i_resolution=1, j_resolution=1)
    # plane_z1 = pv.Plane((0.5, 0.5, 1.0), i_resolution=1, j_resolution=1)

    # plane_z0.save('/Users/yuechengyu/Work/Cardiac/gmsh/unit_cube_Z0.vtp')
    # plane_z1.save('/Users/yuechengyu/Work/Cardiac/gmsh/unit_cube_Z1.vtp')



    # # meshio style:      
    # # *** 'vtu' works
    # # Points of the unit cube
    # points = np.array([
    #     [0.0, 0.0, 0.0],  # Point 0
    #     [1.0, 0.0, 0.0],  # Point 1
    #     [1.0, 1.0, 0.0],  # Point 2
    #     [0.0, 1.0, 0.0],  # Point 3
    #     [0.0, 0.0, 1.0],  # Point 4
    #     [1.0, 0.0, 1.0],  # Point 5
    #     [1.0, 1.0, 1.0],  # Point 6
    #     [0.0, 1.0, 1.0]   # Point 7
    # ])

    # # One hexahedral element connectivity
    # cells = [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))]

    # # Create the volume mesh
    # volume_mesh = meshio.Mesh(points, cells)
    # # Save the volume mesh
    # volume_mesh_path = "unit_cube_volume.vtu"
    # meshio.write(volume_mesh_path, volume_mesh, file_format="vtu")
    # # meshio.write(volume_mesh_path, volume_mesh)

    # # *** Doesn't work with svFSI:
    # # Define the six faces of the cube for surface meshes
    # faces = {
    #     "Z0": ("triangle", np.array([[0, 1, 3], [1, 2, 3]])),
    #     "Z1": ("triangle", np.array([[4, 5, 7], [5, 6, 7]])),
    #     "Y0": ("triangle", np.array([[0, 1, 4], [1, 5, 4]])),
    #     "Y1": ("triangle", np.array([[3, 2, 7], [2, 6, 7]])),
    #     "X0": ("triangle", np.array([[0, 3, 4], [3, 7, 4]])),
    #     "X1": ("triangle", np.array([[1, 2, 5], [2, 6, 5]])),
    # }
    # # faces = {
    # #     "bottom": ("quad", np.array([[0, 1, 2, 3]])),
    # #     "top": ("quad", np.array([[4, 5, 6, 7]])),
    # #     "front": ("quad", np.array([[0, 1, 5, 4]])),
    # #     "back": ("quad", np.array([[3, 2, 6, 7]])),
    # #     "left": ("quad", np.array([[0, 3, 7, 4]])),
    # #     "right": ("quad", np.array([[1, 2, 6, 5]])),
    # # }

    # # Create and save the six surface meshes
    # surface_mesh_paths = {}
    # reader = vtk.vtkSTLReader()
    # writer = vtk.vtkXMLPolyDataWriter()
    # for face_name, (cell_type, cell_data) in faces.items():
    #     cells = [meshio.CellBlock(cell_type, cell_data)]
    #     surface_mesh = meshio.Mesh(points, cells)
    #     # surface_mesh = meshio.Mesh(points, [cell_type, cell_data])
    #     stl_file = f"unit_cube_{face_name}.stl"
    #     # breakpoint()
    #     meshio.write(stl_file, surface_mesh, file_format="stl")
    #     # meshio.write(file_path, surface_mesh)
    #     surface_mesh_paths[face_name] = stl_file

    #     vtp_file = f"unit_cube_{face_name}.vtp"

    #     # # Convert STL to VTP by using 'VTK'
    #     # reader.SetFileName(stl_file)
    #     # writer.SetFileName(vtp_file)
    #     # writer.SetInputConnection(reader.GetOutputPort())
    #     # writer.Write()

    #     # Convert STL to VTP by using 'pyvista'
    #     stl_mesh = pv.read(stl_file)
    #     stl_mesh.save(vtp_file)

    


if __name__ == "__main__":
    main()



