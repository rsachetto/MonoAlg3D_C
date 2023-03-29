import numpy as np
from paraview.util.vtkAlgorithm import (
    VTKPythonAlgorithmBase,
    smdomain,
    smhint,
    smproperty,
    smproxy,
)
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
import os

paraview_plugin_version = '0.1'
alg_input_filetypes = ['alg', 'acm']
alg_extensions = ['alg', 'acm']

def LoadAlgMesh(filePath, output):

    algFile = open(filePath, 'r')

    p_dict = {}

    p_id = 0
    points = []
    cells = []
    cell_arrays = None
    extra_arrays_len = 0
    cell_arrays_inited = False

    _, extension = os.path.splitext(filePath)

    acm_file = False

    if('acm' in extension):
        acm_file = True
        next(algFile) #skip the first line as it is not used in acm files

    for line in algFile:

        spl_line = line.split(',')

        if not cell_arrays_inited:
            if acm_file:
                extra_arrays_len = 1
            else:
                extra_arrays_len = len(spl_line) - 6

            cell_arrays = [[] for i in range(extra_arrays_len)]
            cell_arrays_inited = True

        center_x = float(spl_line[0])
        center_y = float(spl_line[1])
        center_z = float(spl_line[2])

        if extra_arrays_len >=0:
            half_face_x = float(spl_line[3])
            half_face_y = float(spl_line[4])
            half_face_z = float(spl_line[5])

            if(acm_file):
                if bool(spl_line[6]) == False: continue
                cell_arrays[0].append(float(spl_line[8].split()[1]))
            else:
                for i in range(extra_arrays_len):
                    cell_arrays[i].append(float(spl_line[i+6]))

        else:
            half_face_x = float(spl_line[3])
            half_face_y = float(spl_line[3])
            half_face_z = float(spl_line[3])

        center_x_plus = center_x + half_face_x
        center_x_minus = center_x - half_face_x
        center_y_plus = center_y + half_face_y
        center_y_minus = center_y - half_face_y

        center_z_plus = center_z + half_face_z
        center_z_minus = center_z - half_face_z

        point1 = (center_x_minus, center_y_minus, center_z_minus)
        point2 = (center_x_plus, center_y_minus, center_z_minus)
        point3 = (center_x_plus, center_y_plus, center_z_minus)
        point4 = (center_x_minus, center_y_plus, center_z_minus)
        point5 = (center_x_minus, center_y_minus, center_z_plus)
        point6 = (center_x_plus, center_y_minus, center_z_plus)
        point7 = (center_x_plus, center_y_plus, center_z_plus)
        point8 = (center_x_minus, center_y_plus, center_z_plus)

        point1_idx = p_dict.get(point1, -1)
        point2_idx = p_dict.get(point2, -1)
        point3_idx = p_dict.get(point3, -1)
        point4_idx = p_dict.get(point4, -1)
        point5_idx = p_dict.get(point5, -1)
        point6_idx = p_dict.get(point6, -1)
        point7_idx = p_dict.get(point7, -1)
        point8_idx = p_dict.get(point8, -1)


        if point1_idx == -1:
            points.append(list(point1))
            p_dict[point1] = p_id
            point1_idx = p_id
            p_id += 1

        if point2_idx == -1:
            points.append(list(point2))
            p_dict[point2] = p_id
            point2_idx = p_id
            p_id += 1

        if point3_idx == -1:
            points.append(list(point3))
            p_dict[point3] = p_id
            point3_idx = p_id
            p_id += 1

        if point4_idx == -1:
            points.append(list(point4))
            p_dict[point4] = p_id
            point4_idx = p_id
            p_id += 1

        if point5_idx == -1:
            points.append(list(point5))
            p_dict[point5] = p_id
            point5_idx = p_id
            p_id += 1

        if point6_idx == -1:
            points.append(list(point6))
            p_dict[point6] = p_id
            point6_idx = p_id
            p_id += 1

        if point7_idx == -1:
            points.append(list(point7))
            p_dict[point7] = p_id
            point7_idx = p_id
            p_id += 1

        if point8_idx == -1:
            points.append(list(point8))
            p_dict[point8] = p_id
            point8_idx = p_id
            p_id += 1

        cells.append([point1_idx, point2_idx, point3_idx, point4_idx, point5_idx, point6_idx, point7_idx, point8_idx])

    ncells = len(cells)
    npoints = 8

    output.SetPoints(points)

    cell_types = np.array([], dtype=np.ubyte)
    cell_offsets = np.array([], dtype=int)
    cell_conn = np.array([], dtype=int)

    cell_types = np.hstack([cell_types, np.full(ncells, 12, dtype=np.ubyte)])

    offsets = len(cell_conn) + (1 + npoints) * np.arange(ncells, dtype=int)
    cell_offsets = np.hstack([cell_offsets, offsets])

    conn = np.hstack([npoints * np.ones((ncells, 1), dtype=int), cells]).flatten()
    cell_conn = np.hstack([cell_conn, conn])
    output.SetCells(cell_types, cell_offsets, cell_conn)

    # Cell data
    if cell_arrays:
        count = 1
        for data in cell_arrays:
            if count == 1 and acm_file:
                name = "Activaton data"
            else:
                name = "Cell data " + str(count)
            output.CellData.append(np.array(data), name)
            count += 1
    algFile.close()

@smproxy.reader(
    name="alg reader",
    extensions=alg_extensions,
    file_description="monoalg3d-supported files",
    support_reload=False,
)
class AlgReader(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1, outputType="vtkUnstructuredGrid"
        )
        self._filename = None
        self._file_format = None

    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(
        extensions=alg_extensions, file_description="monoalg3d-supported alg file"
    )
    def SetFileName(self, filename):
        if self._filename != filename:
            self._filename = filename
            self.Modified()

    @smproperty.stringvector(name="StringInfo", information_only="1")
    def GetStrings(self):
        return alg_input_filetypes

    @smproperty.stringvector(name="FileFormat", number_of_elements="1")
    @smdomain.xml(
        """
        <StringListDomain name="list">
            <RequiredProperties>
                <Property name="StringInfo" function="StringInfo"/>
            </RequiredProperties>
        </StringListDomain>
        """
    )
    def SetFileFormat(self, file_format):
        if self._file_format != file_format:
            self._file_format = file_format
            self.Modified()

    def RequestData(self, request, inInfoVec, outInfoVec):
        output = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(outInfoVec))
        LoadAlgMesh(self._filename, output)
        return 1

'''
        # Point data
        for name, array in mesh.point_data.items():
            output.PointData.append(array, name)

        # Field data
        for name, array in mesh.field_data.items():
            output.FieldData.append(array, name)

@smproxy.writer(
    name="meshio Writer",
    extensions=meshio_extensions,
    file_description="meshio-supported files",
    support_reload=False,
)
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)
class MeshioWriter(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=1, nOutputPorts=0, inputType="vtkUnstructuredGrid"
        )
        self._filename = None

    @smproperty.stringvector(name="FileName", panel_visibility="never")
    @smdomain.filelist()
    def SetFileName(self, filename):
        if self._filename != filename:
            self._filename = filename
            self.Modified()

    def RequestData(self, request, inInfoVec, outInfoVec):
        mesh = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(inInfoVec[0]))

        # Read points
        points = np.asarray(mesh.GetPoints())

        # Read cells
        # Adapted from test/legacy_reader.py
        cell_conn = mesh.GetCells()
        cell_offsets = mesh.GetCellLocations()
        cell_types = mesh.GetCellTypes()
        cells_dict = {}
        for vtk_cell_type in np.unique(cell_types):
            offsets = cell_offsets[cell_types == vtk_cell_type]
            ncells = len(offsets)
            npoints = cell_conn[offsets[0]]
            array = np.empty((ncells, npoints), dtype=int)
            for i in range(npoints):
                array[:, i] = cell_conn[offsets + i + 1]
            cells_dict[vtk_to_meshio_type[vtk_cell_type]] = array
        cells = [meshio.CellBlock(key, cells_dict[key]) for key in cells_dict]

        # Read point and field data
        # Adapted from test/legacy_reader.py
        def _read_data(data):
            out = {}
            for i in range(data.VTKObject.GetNumberOfArrays()):
                name = data.VTKObject.GetArrayName(i)
                array = np.asarray(data.GetArray(i))
                out[name] = array
            return out

        point_data = _read_data(mesh.GetPointData())
        field_data = _read_data(mesh.GetFieldData())

        # Read cell data
        cell_data_flattened = _read_data(mesh.GetCellData())
        cell_data = {}
        for name, array in cell_data_flattened.items():
            cell_data[name] = []
            for cell_type in cells_dict:
                vtk_cell_type = meshio_to_vtk_type[cell_type]
                mask_cell_type = cell_types == vtk_cell_type
                cell_data[name].append(array[mask_cell_type])

        # Use meshio to write mesh
        meshio.write_points_cells(
            self._filename,
            points,
            cells,
            point_data=point_data,
            cell_data=cell_data,
            field_data=field_data,
        )
        return 1

    def Write(self):
        self.Modified()
        self.Update()
'''
