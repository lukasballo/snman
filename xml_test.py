import xml.etree.ElementTree as ET

# Create the overall structure
network = ET.Element('network')
attributes = ET.SubElement(network, 'attributes')
ET.SubElement(
    attributes, 'attribute', attrib={'name': 'coordinateReferenceSystem', 'class': 'java.lang.string'
    }).text = '2056'
nodes = ET.SubElement(network, 'nodes')
links = ET.SubElement(network, 'links')

# Add Nodes


# Save into xml file
ET.indent(network)
output = ET.tostring(network)
with open("GFG.xml", "wb") as f:
    f.write(output)