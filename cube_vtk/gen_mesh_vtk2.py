import meshio 
import numpy as np
import vtk
import pyvista as pv
import os

from collections import defaultdict
from vtk.util.numpy_support import numpy_to_vtk as n2v

def main():

    # Define the eight vertices of a unit cube
    points = vtk.vtkPoints()
    points.InsertNextPoint(0.0, 0.0, 0.0)
    points.InsertNextPoint(1.0, 0.0, 0.0)
    points.InsertNextPoint(1.0, 1.0, 0.0)
    points.InsertNextPoint(0.0, 1.0, 0.0)
    points.InsertNextPoint(0.0, 0.0, 1.0)
    points.InsertNextPoint(1.0, 0.0, 1.0)
    points.InsertNextPoint(1.0, 1.0, 1.0)
    points.InsertNextPoint(0.0, 1.0, 1.0)

    for i in range(1):
        # Create a quad on the four points
        quad = vtk.vtkQuad()
        quad.GetPointIds().SetId(i*4+0, i*4+0)
        quad.GetPointIds().SetId(i*4+1, i*4+1)
        quad.GetPointIds().SetId(i*4+2, i*4+2)
        quad.GetPointIds().SetId(i*4+3, i*4+3)

        # Create a cell array to store the quad in
        quads = vtk.vtkCellArray()
        quads.InsertNextCell(quad)

        # Create a polydata to store everything in
        polydata = vtk.vtkPolyData()

        # Add the points and quads to the dataset
        polydata.SetPoints(points)
        polydata.SetPolys(quads)

        # Write to VTP file
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(f'Z{i}.vtp')
        writer.SetInputData(polydata)
        writer.Update()
        writer.Write()


    # Define the six faces of the cube using the vertex points
    # faces = vtk.vtkCellArray()
    face_ids = [
        # [0, 3, 2, 1],  # bottom face
        [0, 1, 2, 3],  # bottom face
        [4, 5, 6, 7],  # top face
        [0, 1, 5, 4],  # front face
        [2, 3, 7, 6],  # back face
        [0, 4, 7, 3],  # left face
        [1, 2, 6, 5]   # right face
    ]
    # for id_list in face_ids:
    #     faces.InsertNextCell(4, id_list)

    for i in range(2):
        face = vtk.vtkCellArray()
        face.InsertNextCell(4, face_ids[i])

        # Create a polydata object for the surface mesh
        surface_mesh = vtk.vtkPolyData()
        surface_mesh.SetPoints(points)
        surface_mesh.SetPolys(face)

        # Write the surface mesh as a VTP file
        writer_surface = vtk.vtkXMLPolyDataWriter()
        writer_surface.SetFileName(f'unit_cube_Z{i}.vtp')
        writer_surface.SetInputData(surface_mesh)
        # writer_surface.SetDataModeToBinary()
        # writer_surface.SetEncodeAppendedData(False)
        writer_surface.Update()
        writer_surface.Write()

    # Define the hexahedral cell for the volume mesh
    hexahedron = vtk.vtkHexahedron()
    for i in range(8):
        hexahedron.GetPointIds().SetId(i, i)

    # Create an unstructured grid object for the volume mesh
    volume_mesh = vtk.vtkUnstructuredGrid()
    volume_mesh.SetPoints(points)
    volume_mesh.InsertNextCell(hexahedron.GetCellType(), hexahedron.GetPointIds())

    # Write the volume mesh as a VTU file
    writer_volume = vtk.vtkXMLUnstructuredGridWriter()
    writer_volume.SetFileName('unit_cube_volume.vtu')
    writer_volume.SetInputData(volume_mesh)
    writer_volume.Write()











    


if __name__ == "__main__":
    main()



