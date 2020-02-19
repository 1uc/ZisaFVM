import os
import h5py

import xml.etree.ElementTree as ET
import xml.dom.minidom


def xdmf_type_info(type):
    if type == "int":
        number_type = "Int"
        precision = "4"
    elif type == "double":
        number_type = "Float"
        precision = "8"
    else:
        raise Exception(f"Unknown `type`. [{type}]")

    return number_type, precision


def data_item(parrent, dims, type, format, body):

    number_type, precision = xdmf_type_info(type)

    di = ET.SubElement(
        parrent,
        "DataItem",
        Dimensions=" ".join(dims),
        NumberType=number_type,
        Precision=precision,
        Format=format,
    )
    di.text = body

    return di


def component(parrent, n_cells, key, body):
    attr = ET.SubElement(
        parrent, "Attribute", Name=key, AttributeType="Scalar", Center="Cell"
    )
    data_item(attr, (n_cells,), "double", "HDF", body)


def extract_grid_info(grid_file):
    with h5py.File(grid_file, "r") as h5_grid:
        n_cells = str(int(h5_grid["n_cells"][()]))
        n_vertices = str(int(h5_grid["n_vertices"][()]))
        max_neighbours = str(int(h5_grid["max_neighbours"][()]))

    if max_neighbours == "4":
        element = "Tetrahedron"
    elif max_neighbours == "3":
        element = "Triangle"
    else:
        raise Exception(f"Unknown number of neighbours. [{max_neighbours}]")

    return n_cells, n_vertices, max_neighbours, element


def extract_time(data_file):
    with h5py.File(data_file, "r") as h5:
        return str(float(h5["time"][()]))


def generate_xdmf(grid_file, data_files, components):
    n_cells, n_vertices, max_neighbours, element = extract_grid_info(grid_file)

    doc = ET.Element("Xdmf", Version="2.0")
    domain = ET.SubElement(doc, "Domain")
    top = ET.SubElement(
        domain, "Topology", TopologyType=element, NumberOfElements=n_cells
    )

    vertex_indices = f"{os.path.basename(grid_file)}:/vertex_indices"
    data_item(top, (n_cells, max_neighbours), "int", "HDF", vertex_indices)

    geo = ET.SubElement(domain, "Geometry", GeometryType="XYZ")
    vertices = f"{os.path.basename(grid_file)}:/vertices"
    data_item(geo, (n_vertices, "3"), "double", "HDF", vertices)

    temporal = ET.SubElement(
        domain, "Grid", GridType="Collection", CollectionType="Temporal"
    )

    for df in data_files:
        snapshot = ET.SubElement(temporal, "Grid", GridType="Uniform")
        ET.SubElement(snapshot, "Topology", Reference="/Xdmf/Domain/Topology[1]")
        ET.SubElement(snapshot, "Geometry", Reference="/Xdmf/Domain/Geometry[1]")
        ET.SubElement(snapshot, "Time", Value=extract_time(df))

        for c in components:
            component(snapshot, n_cells, c, f"{os.path.basename(df)}:/{c}")

    return doc


def xml_to_string(xml_doc):
    return xml.dom.minidom.parseString(ET.tostring(xml_doc)).toprettyxml(indent="  ")
