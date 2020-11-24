import os, sys
from xml.dom import minidom
# python write_calibration_xml.py ../../Detector/DetFCCeeECalInclined/compact/FCCee_ECalBarrel.xml

input_xml_path = sys.argv[1]
output_xml_path = input_xml_path.replace(".xml", "_calibration.xml")

input_xml = minidom.parse(input_xml_path)
#print input_xml.toprettyxml()
numberOfLayer = 0
#print input_xml.getElementsByTagName('lccdd')
for nodeList in input_xml.getElementsByTagName('lccdd'):
    for node in nodeList.childNodes:
        if node.localName == 'detectors':
            for subnode in node.childNodes:
                if subnode.localName == 'detector':
                    print subnode.localName
                    for subsubnode in subnode.childNodes:
                        if subsubnode.localName == 'calorimeter':
                            print "    ", subsubnode.localName
                            for subsubsubnode in subsubnode.childNodes:
                                if subsubsubnode.localName == 'readout':
                                    print "        ", subsubsubnode.localName
                                    #print subsubsubnode.getAttribute('sensitive')
                                    subsubsubnode.setAttribute('sensitive', 'true')
                                    #print subsubsubnode.getAttribute('sensitive')
                                if subsubsubnode.localName == 'passive':
                                    print "        ", subsubsubnode.localName
                                    for subsubsubsubnode in subsubsubnode.childNodes:
                                        if subsubsubsubnode.localName in ['inner', 'glue', 'outer']:
                                            print "            ", subsubsubsubnode.localName
                                            #print subsubsubsubnode.getAttribute('sensitive')
                                            subsubsubsubnode.setAttribute('sensitive', 'true')
                                            #print subsubsubsubnode.getAttribute('sensitive')
                                if subsubsubnode.localName == 'layers':
                                    for subsubsubsubnode in subsubsubnode.childNodes:
                                        if subsubsubsubnode.localName == 'layer':
                                            numberOfLayer += int(subsubsubsubnode.getAttribute('repeat'))
with open(output_xml_path, "w") as f:
    input_xml.writexml(f)
print output_xml_path, " written." 
print "Number of layers: %d"%numberOfLayer

print "sed -i 's/numLayers.*,/numLayers = %d,/' fcc_ee_samplingFraction_inclinedEcal.py"%numberOfLayer
os.system("sed -i 's/numLayers.*,/numLayers = %d,/' fcc_ee_samplingFraction_inclinedEcal.py"%numberOfLayer)
sys.exit(numberOfLayer)
