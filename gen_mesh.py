import meshio 
import numpy as np
import vtk
import pyvista as pv
import os

from collections import defaultdict
from vtk.util.numpy_support import numpy_to_vtk as n2v

def threshold(inp, t, name):
    """
    Threshold according to cell array
    Args:
        inp: InputConnection
        t: BC_FaceID
        name: name in cell data used for thresholding
    Returns:
        reader, point data
    """
    thresh = vtk.vtkThreshold()
    thresh.SetInputData(inp)
    thresh.SetInputArrayToProcess(0, 0, 0, 1, name)
    thresh.ThresholdBetween(t, t)
    thresh.Update()
    return thresh

def main():

    # read volume mesh in vtk
    fname = 'unit_cube_volume.vtu'
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fname)
    reader.Update()
    vol = reader.GetOutput()

    surf_dict = defaultdict(list)
    surf_ids = {}
    points_inlet = []
    for f in ["solid", "fluid"]:
        # select sub-mesh
        vol_f = threshold(vol, 1, "ids_" + f).GetOutput()

        # reset global ids
        n_array = n2v(np.arange(vol_f.GetNumberOfPoints()).astype(np.int32) + 1)
        e_array = n2v(np.arange(vol_f.GetNumberOfCells()).astype(np.int32) + 1)
        n_array.SetName("GlobalNodeID")
        e_array.SetName("GlobalElementID")
        vol_f.GetPointData().AddArray(n_array)
        vol_f.GetCellData().AddArray(e_array)

        # map point data to cell data
        p2c = vtk.vtkPointDataToCellData()
        p2c.SetInputData(vol_f)
        p2c.PassPointDataOn()
        p2c.Update()
        vol_f = p2c.GetOutput()

        # extract surfaces
        extract = vtk.vtkGeometryFilter()
        extract.SetInputData(vol_f)
        # extract.SetNonlinearSubdivisionLevel(0)
        extract.Update()
        surfaces = extract.GetOutput()

        # threshold surfaces
        for name in surf_dict.keys():
            # interior quad elements
            if self.p["n_seg"] == 1 and "_zero" in name:
                # threshold circle segments first
                thresh = vtk.vtkThreshold()
                thresh.SetInputData(vol_f)
                thresh.SetInputArrayToProcess(0, 0, 0, 1, "ids_" + name[0] + "_seg")
                thresh.SetUpperThreshold(1)
                thresh.SetLowerThreshold(1)
                thresh.Update()

                # extract surfaces
                extract = vtk.vtkGeometryFilter()
                extract.SetInputData(thresh.GetOutput())
                extract.SetNonlinearSubdivisionLevel(0)
                extract.Update()
                inp = extract.GetOutput()
            else:
                inp = surfaces

            # select only current surface
            thresh = vtk.vtkThreshold()
            thresh.SetInputData(inp)
            thresh.SetInputArrayToProcess(0, 0, 0, 0, "ids_" + name)
            thresh.SetUpperThreshold(1)
            thresh.SetLowerThreshold(1)
            thresh.Update()
            surf = thresh.GetOutput()
            if surf.GetNumberOfPoints() > 0:
                surf = clean(extract_surface(surf))

            fout = os.path.join(self.p["f_out"], f, "mesh-surfaces", name + ".vtp")
            write_geo(fout, extract_surface(surf))

            # get new GlobalNodeIDs of surface points
            surf_ids[f + "_" + name] = v2n(
                surf.GetPointData().GetArray("GlobalNodeID")
            ).tolist()

            # store inlet points (to calculate flow profile later)
            if f == "fluid" and name == "start":
                points_inlet = v2n(surf.GetPoints().GetData())

        # export volume mesh
        write_geo(os.path.join(self.p["f_out"], f, "mesh-complete.mesh.vtu"), vol_f)



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



