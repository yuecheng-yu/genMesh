import meshio 
import numpy as np
import vtk
import pyvista as pv
import os

def main():

    # Create points for the cube
    points = vtk.vtkPoints()
    points.InsertNextPoint(0, 0, 0)
    points.InsertNextPoint(1, 0, 0)
    points.InsertNextPoint(1, 1, 0)
    points.InsertNextPoint(0, 1, 0)
    points.InsertNextPoint(0, 0, 1)
    points.InsertNextPoint(1, 0, 1)
    points.InsertNextPoint(1, 1, 1)
    points.InsertNextPoint(0, 1, 1)

    # Create a hexahedron from the points
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



    # Shrink the cells of the grid slightly to separate them
    shrink = vtk.vtkShrinkFilter()
    shrink.SetInputData(volume_mesh)
    shrink.SetShrinkFactor(0.99)
    shrink.Update()

    # Use vtkDataSetSurfaceFilter to extract the surface from the shrunken unstructured grid
    surface_filter = vtk.vtkDataSetSurfaceFilter()
    surface_filter.SetInputConnection(shrink.GetOutputPort())
    surface_filter.Update()

    # The result is a vtkPolyData object representing the surface mesh
    surface_mesh = surface_filter.GetOutput()

    # Create cell arrays to hold the top and bottom faces
    top_cells = vtk.vtkCellArray()
    bottom_cells = vtk.vtkCellArray()

    # Iterate over each polygon in the surface mesh
    for i in range(surface_mesh.GetNumberOfCells()):
        cell = surface_mesh.GetCell(i)
        cell_points = cell.GetPoints()
        cell_point_ids = cell.GetPointIds()
        z_values = [cell_points.GetPoint(j)[2] for j in range(cell.GetNumberOfPoints())]

        # If all z-values are close to 1, it's the top cell; if close to 0, it's the bottom cell
        if all(z > 0.9 for z in z_values):
            point_ids = [cell_point_ids.GetId(j) for j in range(cell_point_ids.GetNumberOfIds())]
            top_cells.InsertNextCell(len(point_ids), point_ids)
        elif all(z < 0.1 for z in z_values):
            point_ids = [cell_point_ids.GetId(j) for j in range(cell_point_ids.GetNumberOfIds())]
            bottom_cells.InsertNextCell(len(point_ids), point_ids)

    # Create polydata for the top and bottom surfaces
    top_surface = vtk.vtkPolyData()
    bottom_surface = vtk.vtkPolyData()

    top_surface.SetPoints(surface_mesh.GetPoints())
    bottom_surface.SetPoints(surface_mesh.GetPoints())

    top_surface.SetPolys(top_cells)
    bottom_surface.SetPolys(bottom_cells)

    # Write the top surface to a file
    writer_top = vtk.vtkXMLPolyDataWriter()
    writer_top.SetFileName('cube_top_surface.vtp')
    writer_top.SetInputData(top_surface)
    writer_top.Write()

    # Write the bottom surface to a file
    writer_bottom = vtk.vtkXMLPolyDataWriter()
    writer_bottom.SetFileName('cube_bottom_surface.vtp')
    writer_bottom.SetInputData(bottom_surface)
    writer_bottom.Write()




    # # Create a cell array to store the hexahedron
    # cells = vtk.vtkCellArray()
    # cells.InsertNextCell(hexahedron)

    # # Create an unstructured grid to represent the cube
    # cube = vtk.vtkUnstructuredGrid()
    # cube.SetPoints(points)
    # cube.InsertNextCell(hexahedron.GetCellType(), hexahedron.GetPointIds())

    # # Use vtkDataSetSurfaceFilter to extract the surface from the unstructured grid
    # surface_filter = vtk.vtkDataSetSurfaceFilter()
    # surface_filter.SetInputData(cube)
    # surface_filter.Update()

    # # The result is a vtkPolyData object representing the surface mesh
    # surface_mesh = surface_filter.GetOutput()

    # # Extract the top and bottom surfaces
    # top_cells = vtk.vtkCellArray()
    # bottom_cells = vtk.vtkCellArray()

    # # Iterate over each cell (face) in the surface mesh
    # for i in range(surface_mesh.GetNumberOfCells()):
    #     breakpoint()
    #     cell = surface_mesh.GetCell(i)
    #     cell_points = cell.GetPoints()
    #     cell_ids = cell.GetPointIds()
        
    #     # Assume all z-values are the same for each cell's points in a perfect cube
    #     z_value = cell_points.GetPoint(0)[2]
        
    #     # Check if the cell is a top face or bottom face based on the z-value
    #     if z_value == 1.0:
    #         # Create a new cell for the top face
    #         top_cell = vtk.vtkPolygon()
    #         top_cell.GetPointIds().SetNumberOfIds(cell_ids.GetNumberOfIds())
    #         for j in range(cell_ids.GetNumberOfIds()):
    #             top_cell.GetPointIds().SetId(j, cell_ids.GetId(j))
    #         top_cells.InsertNextCell(top_cell)
    #     elif z_value == 0.0:
    #         # Create a new cell for the bottom face
    #         bottom_cell = vtk.vtkPolygon()
    #         bottom_cell.GetPointIds().SetNumberOfIds(cell_ids.GetNumberOfIds())
    #         for j in range(cell_ids.GetNumberOfIds()):
    #             bottom_cell.GetPointIds().SetId(j, cell_ids.GetId(j))
    #         bottom_cells.InsertNextCell(bottom_cell)

    # # Create polydata for the top and bottom surfaces
    # top_surface = vtk.vtkPolyData()
    # bottom_surface = vtk.vtkPolyData()

    # top_surface.SetPoints(points)
    # bottom_surface.SetPoints(points)

    # top_surface.SetPolys(top_cells)
    # bottom_surface.SetPolys(bottom_cells)

    # # Write the top surface to a file
    # writer_top = vtk.vtkXMLPolyDataWriter()
    # writer_top.SetFileName('cube_top_surface.vtp')
    # writer_top.SetInputData(top_surface)
    # writer_top.Write()

    # # Write the bottom surface to a file
    # writer_bottom = vtk.vtkXMLPolyDataWriter()
    # writer_bottom.SetFileName('cube_bottom_surface.vtp')
    # writer_bottom.SetInputData(bottom_surface)
    # writer_bottom.Write()


    

if __name__ == "__main__":
    main()



