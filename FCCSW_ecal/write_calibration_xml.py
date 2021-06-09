import os, sys
from xml.dom import minidom
# python write_calibration_xml.py ../../Detector/DetFCCeeECalInclined/compact/FCCee_ECalBarrel.xml
# careful: if you use the FCC_DETECTOR environment variable, you have to recompile after modifying the xml so that they go in the right place

input_xml_path = sys.argv[1]
output_xml_path_sf = input_xml_path.replace(".xml", "_calibration.xml")
output_xml_path_upstream = input_xml_path.replace(".xml", "_upstream.xml")

list_of_pair_layerThickness_numberOfLayer = []

input_xml = minidom.parse(input_xml_path)
#print input_xml.toprettyxml()
numberOfLayer = 0
#print input_xml.getElementsByTagName('lccdd')
for nodeList in input_xml.getElementsByTagName('lccdd'):
    for node in nodeList.childNodes:
        if node.localName == 'detectors':
            for subnode in node.childNodes:
                if subnode.localName == 'detector':
                    print(subnode.localName)
                    for subsubnode in subnode.childNodes:
                        if subsubnode.localName == 'calorimeter':
                            print("    ", subsubnode.localName)
                            for subsubsubnode in subsubnode.childNodes:
                                if subsubsubnode.localName == 'readout':
                                    print("        ", subsubsubnode.localName)
                                    #print subsubsubnode.getAttribute('sensitive')
                                    subsubsubnode.setAttribute('sensitive', 'true') # here we change the readout as sensitive
                                    #print subsubsubnode.getAttribute('sensitive')
                                if subsubsubnode.localName == 'passive':
                                    print("        ", subsubsubnode.localName)
                                    for subsubsubsubnode in subsubsubnode.childNodes:
                                        if subsubsubsubnode.localName in ['inner', 'glue', 'outer']:
                                            print("            ", subsubsubsubnode.localName)
                                            #print subsubsubsubnode.getAttribute('sensitive')
                                            subsubsubsubnode.setAttribute('sensitive', 'true') # here we change the absorber into sensitive material
                                            #print subsubsubsubnode.getAttribute('sensitive')
                                if subsubsubnode.localName == 'layers':
                                    for subsubsubsubnode in subsubsubnode.childNodes:
                                        if subsubsubsubnode.localName == 'layer':
                                            numberOfLayer += int(subsubsubsubnode.getAttribute('repeat'))
                                            list_of_pair_layerThickness_numberOfLayer.append([subsubsubsubnode.getAttribute('thickness').split('*')[0], subsubsubsubnode.getAttribute('repeat')])
with open(output_xml_path_sf, "w") as f:
    input_xml.writexml(f)
print(output_xml_path_sf, " written.") 
print("Number of layers: %d"%numberOfLayer)
print("Layer layout {depth : number}: ", list_of_pair_layerThickness_numberOfLayer) 

# modify the number of layer in sampling fraction config
os.system("sed -i 's/numLayers.*,/numLayers = %d,/' fcc_ee_samplingFraction_inclinedEcal.py"%numberOfLayer)
print("fcc_ee_samplingFraction_inclinedEcal.py modified")

# modify the layer layout in plot_sampling_fraction script
os.system("sed -i 's/default = \[1\] \*.*,/default = \[1\] \* %d,/' FCC_calo_analysis_cpp/plot_samplingFraction.py"%numberOfLayer)
os.system("sed -i 's/totalNumLayers\", default = .*,/totalNumLayers\", default = %d,/' FCC_calo_analysis_cpp/plot_samplingFraction.py"%numberOfLayer)
string_for_layerWidth = ""
for pair_layerThickness_numberOfLayer in list_of_pair_layerThickness_numberOfLayer:
    string_for_layerWidth += "[%f] * %d + "%(float(pair_layerThickness_numberOfLayer[0]), int(pair_layerThickness_numberOfLayer[1])) 
string_for_layerWidth = string_for_layerWidth[0:-2]
os.system("sed -i 's/layerWidth\", default = .*,/layerWidth\", default = %s,/' FCC_calo_analysis_cpp/plot_samplingFraction.py"%string_for_layerWidth)
print("FCC_calo_analysis_cpp/plot_samplingFraction.py modified")

# modify the number of layer in fcc_ee_upstreamMaterial_inclinedEcal.py
os.system("sed -i 's/numLayers.*,/numLayers = %d,/' fcc_ee_upstreamMaterial_inclinedEcal.py"%numberOfLayer)
print("fcc_ee_upstreamMaterial_inclinedEcal.py modified")

# Write upstream correction xml
# Re-make absorber and readout not sensitive
for nodeList in input_xml.getElementsByTagName('lccdd'):
    for node in nodeList.childNodes:
        if node.localName == 'detectors':
            for subnode in node.childNodes:
                if subnode.localName == 'detector':
                    for subsubnode in subnode.childNodes:
                        if subsubnode.localName == 'calorimeter':
                            for subsubsubnode in subsubnode.childNodes:
                                if subsubsubnode.localName == 'readout':
                                    subsubsubnode.setAttribute('sensitive', 'false')
                                if subsubsubnode.localName == 'passive':
                                    for subsubsubsubnode in subsubsubnode.childNodes:
                                        if subsubsubsubnode.localName in ['inner', 'glue', 'outer']:
                                            subsubsubsubnode.setAttribute('sensitive', 'false') # here we change the abnsorber into sensitive material!
                        if subsubnode.localName == 'cryostat':
                            for subsubsubnode in subsubnode.childNodes:
                                if subsubsubnode.localName == 'front':
                                    subsubsubnode.setAttribute('sensitive', 'true')

with open(output_xml_path_upstream, "w") as f:
    input_xml.writexml(f)
print(output_xml_path_upstream, " written.") 


