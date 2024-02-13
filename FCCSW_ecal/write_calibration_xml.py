import os
import sys
from xml.dom import minidom
# python write_calibration_xml.py ../../k4geo/FCCee/ALLEGRO/compact/ALLEGRO_o1_v02/ECalBarrel_thetamodulemerged.xml

input_xml_path = sys.argv[1]
output_xml_path_sf = input_xml_path.replace(".xml", "_calibration.xml")
output_xml_path_upstream = input_xml_path.replace(".xml", "_upstream.xml")
detDim_xml_path = os.path.join(os.path.dirname(input_xml_path), "DectDimensions.xml")

list_of_pair_layerThickness_numberOfLayer = []

input_xml = minidom.parse(input_xml_path)
# print input_xml.toprettyxml()
numberOfLayer = 0
n_phi_bins = 0
n_modules = 0
eta_bin_size = 0
eta_bin_size_eval = 0
theta_bin_size = 0
theta_bin_size_eval = 0
original_cryo_back_size = 0
# print input_xml.getElementsByTagName('lccdd')

# write xml file for calculation of sampling fractions, turning passive material into active
for nodeList in input_xml.getElementsByTagName('lccdd'):
    for node in nodeList.childNodes:
        # get the cryostat back size
        if node.localName == 'define':
            for subnode in node.childNodes:
                if subnode.localName == 'constant' and subnode.getAttribute('name') == 'CryoBarrelBackCold':
                    original_cryo_back_size_str = subnode.getAttribute('value')
                    if not original_cryo_back_size_str.split("*")[1] == 'mm':
                        print("Error in orignal xml file, cryo thickness expected in mm, exiting...")
                        sys.exit(1)
                    original_cryo_back_size = int(original_cryo_back_size_str.split("*")[0])
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
                                    # print subsubsubnode.getAttribute('sensitive')
                                    subsubsubnode.setAttribute('sensitive', 'true')  # here we change the readout as sensitive
                                    # print subsubsubnode.getAttribute('sensitive')
                                if subsubsubnode.localName == 'passive':
                                    print("        ", subsubsubnode.localName)
                                    for subsubsubsubnode in subsubsubnode.childNodes:
                                        if subsubsubsubnode.localName in ['inner', 'innerMax', 'glue', 'outer']:
                                            print("            ", subsubsubsubnode.localName)
                                            # print subsubsubsubnode.getAttribute('sensitive')
                                            subsubsubsubnode.setAttribute('sensitive', 'true')  # here we change the absorber into sensitive material
                                            # print subsubsubsubnode.getAttribute('sensitive')
                                if subsubsubnode.localName == 'layers':
                                    for subsubsubsubnode in subsubsubnode.childNodes:
                                        if subsubsubsubnode.localName == 'layer':
                                            numberOfLayer += int(subsubsubsubnode.getAttribute('repeat'))
                                            list_of_pair_layerThickness_numberOfLayer.append([subsubsubsubnode.getAttribute('thickness').split('*')[0], subsubsubsubnode.getAttribute('repeat')])
        if node.localName == 'readouts':
            for subnode in node.childNodes:
                if subnode.localName == 'readout' and subnode.getAttribute('name') == 'ECalBarrelPhiEta':
                    for subsubnode in subnode.childNodes:
                        if subsubnode.localName == 'segmentation':
                            n_phi_bins = int(subsubnode.getAttribute('phi_bins'))
                            eta_bin_size = subsubnode.getAttribute('grid_size_eta')
                            eta_bin_size_eval = eval(eta_bin_size)
                if subnode.localName == 'readout' and subnode.getAttribute('name') == 'ECalBarrelModuleThetaMerged':
                    for subsubnode in subnode.childNodes:
                        if subsubnode.localName == 'segmentation':
                            n_modules = int(subsubnode.getAttribute('nModules'))
                            theta_bin_size = subsubnode.getAttribute('grid_size_theta')
                            theta_bin_size_eval = eval(theta_bin_size)

with open(output_xml_path_sf, "w") as f:
    input_xml.writexml(f)
print(output_xml_path_sf, " written.")


# while parsing the input xml we have found out number and depth of layers, and theta/eta grid size and number of phi bins/modules
# so we update the .py files accordingly
print("Number of layers: %d" % numberOfLayer)
print("Layer layout {depth : number}: ", list_of_pair_layerThickness_numberOfLayer)

# modify the number of layer in sampling fraction and upstream config files
os.system("sed -i 's/numLayers.*,/numLayers=%d,/' fcc_ee_samplingFraction_inclinedEcal.py fcc_ee_upstream_inclinedEcal.py" % numberOfLayer)
print("number of layers updated in fcc_ee_samplingFraction_inclinedEcal.py and fcc_ee_upstream_inclinedEcal.py")

# modify the layer layout in plot_sampling_fraction script
os.system("sed -i 's/default = \[1\] \*.*,/default = \[1\] \* %d,/' FCC_calo_analysis_cpp/plot_samplingFraction.py" % numberOfLayer)
os.system("sed -i 's/totalNumLayers\", default = .*,/totalNumLayers\", default = %d,/' FCC_calo_analysis_cpp/plot_samplingFraction.py" % numberOfLayer)
string_for_layerWidth = ""
for pair_layerThickness_numberOfLayer in list_of_pair_layerThickness_numberOfLayer:
    string_for_layerWidth += "[%f] * %d + " % (float(pair_layerThickness_numberOfLayer[0]), int(pair_layerThickness_numberOfLayer[1]))
string_for_layerWidth = string_for_layerWidth[0:-2]
os.system("sed -i 's/layerWidth\", default = .*,/layerWidth\", default = %s,/' FCC_calo_analysis_cpp/plot_samplingFraction.py" % string_for_layerWidth)
print("number and thickness of layers updated in FCC_calo_analysis_cpp/plot_samplingFraction.py")

# modify the number of layers in runUpstreamSlidingWindowAndCaloSim
os.system("sed -i 's/numLayers.*,/numLayers=[%d],/' run*.py" % numberOfLayer)
os.system("sed -i 's/lastLayerIDs.*,/lastLayerIDs=[%d],/' run*.py" % (numberOfLayer - 1))
os.system("sed -i 's/activeVolumesNumber.*,/activeVolumesNumber=%d,/' run*.py" % numberOfLayer)
print("number of layers updated in run*.py")

# modify the number of layers in noise_map.py
os.system("sed -i 's/numRadialLayers.*,/numRadialLayers = %d,/' noise_map.py" % numberOfLayer)
os.system("sed -i 's/activeVolumesNumbers.*,/activeVolumesNumbers = [%d],/' noise_map.py neighbours.py" % numberOfLayer)
print("number of layers updated in noise_map.py and neighbours.py")

# modify the tower definition in clustering algorithms
if eta_bin_size_eval > 0 and n_phi_bins > 0:
    os.system("sed -i 's/deltaEtaTower.*$/deltaEtaTower=%s, deltaPhiTower=2 * _pi \/ %d.,/'  run*.py" % (eta_bin_size, n_phi_bins))
    print("eta-phi tower size updated in run*.py")
else:
    print("eta_bin_size or n_phi_bins are 0 => not using an eta-phi grid readout => not updating eta-phi tower size in run*.py")

# the *4 for theta is because the grid size reflects the width of the strip cells, but the
# tower sizes are defined in units of big cells (corresponding to 4 cells)
if theta_bin_size_eval > 0 and n_modules > 0:
    os.system("sed -i 's#deltaThetaTower.*$#deltaThetaTower=4 * %s, deltaPhiTower=2 * 2 * pi \/ %d.,#'  run_thetamodulemerged.py" % (theta_bin_size, n_modules))
    print("theta-phi tower size updated in run_thetamodulemerged.py")


# Write upstream correction xml
# Re-make absorber and readout not sensitive, make cryostat sensitive and enlarge it to 1.1 m
# First find the detector dimension to enlarge it for the new big cryostat
detDim_xml = minidom.parse(detDim_xml_path)
BarECal_rmax = '0'
for nodeList in detDim_xml.getElementsByTagName('lccdd'):
    for node in nodeList.childNodes:
        if node.localName == 'define':
            for subnode in node.childNodes:
                if subnode.localName == 'constant' and subnode.getAttribute('name') == 'BarECal_rmax':
                    BarECal_rmax = subnode.getAttribute('value')


new_cryo_back_size = 1100  # mm
BarECal_rmax_int = int(BarECal_rmax.split('*')[0])
if not BarECal_rmax.split('*')[1] == 'mm':
    print("Error: dimensions from DectDimensions.xml should be in mm, exiting...")
    sys.exit(1)
new_BarECal_rmax_int = BarECal_rmax_int + (new_cryo_back_size - original_cryo_back_size)

added_rmax = False
for nodeList in input_xml.getElementsByTagName('lccdd'):
    for node in nodeList.childNodes:
        if node.localName == 'define':
            for subnode in node.childNodes:
                if not added_rmax:
                    rMaxNode = input_xml.createElement('constant')
                    rMaxNode.setAttribute('name', 'BarECal_rmax')
                    rMaxNode.setAttribute('value', '%d*mm' % new_BarECal_rmax_int)
                    node.insertBefore(rMaxNode, subnode)
                    added_rmax = True
                if subnode.localName == 'constant' and subnode.getAttribute('name') == "CryoBarrelBackCold":
                    subnode.setAttribute('value', "%d*mm" % new_cryo_back_size)
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
                                        if subsubsubsubnode.localName in ['inner', 'innerMax', 'glue', 'outer']:
                                            subsubsubsubnode.setAttribute('sensitive', 'false')  # here we change the absorber into sensitive material!
                        if subsubnode.localName == 'cryostat':
                            for subsubsubnode in subsubnode.childNodes:
                                if subsubsubnode.localName in ['front', 'side', 'back']:
                                    subsubsubnode.setAttribute('sensitive', 'true')

with open(output_xml_path_upstream, "w") as f:
    input_xml.writexml(f)
print(output_xml_path_upstream, " written.")
