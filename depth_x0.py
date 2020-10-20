from math import *

#FIXME consider presampler, refine cell footprint on circle

# X0 units are in cm from PDG, I use 1/X0 provided by pdg, assuming it mention what thickness of material corresponds to 1 X0 i
# (units 1/cm)
lar_x0 = 1/14.0

lead_x0 = 1/0.5612
W_x0 = 1/0.3504
#absorber_x0 = lead_x0
absorber_x0 = W_x0
glue_x0 = 1/20.209 #LArCaloGlue from geant4 xml dd4hep
steel_x0 = 1/1.775 # lArCaloSteel from geant4 xml dd4hep
cryostat_x0 = 1/8.897 # Aluminiium 
pcb_x0 = 1/15.92 # FR4 (not sure) 

lar_thickness = 1.5 # mm
pcb_thickness = 1.2 # mm
# absorbers
glue_thickness = 0.2 # mm consider immediately the two sides
steel_thickness = 0.4 # mm consider immediately the two sides
absorber_thickness = 1.4 # mm there is only one (sandwiched)

cryostat_thickness_in = 50 # mm
cryostat_thickness_out = 100 # mm

inner_margin_for_services = 10 # mm
outer_margin_for_services = 50 # mm

inclination = 50 #degrees w.r.t. the radial direction

pre_sampler_thickness = 26.9 # mm. without taking into account inclination, number is its bare length

x0_toReach = 22
x0 = 0

cell_size = (lar_thickness * 2 + absorber_thickness + glue_thickness + steel_thickness + pcb_thickness)

# compute the X0 from services, cryostat etc
x0_services = ((cryostat_thickness_in + cryostat_thickness_out) * cryostat_x0 + (inner_margin_for_services + outer_margin_for_services) * lar_x0) / 10.0
x0 += x0_services

# Compute X0 when crossing one cell(abs gap pcb gap), taking into account the inclination w.r.t. radial angle
factor_length_traversed_particle_perpendicular = 1/cos(radians(90 - inclination))
x0_per_normal_cell = (lar_thickness * lar_x0 * 2 + absorber_thickness * absorber_x0 + glue_thickness * glue_x0 + steel_thickness * steel_x0 + pcb_thickness * pcb_x0) * factor_length_traversed_particle_perpendicular / 10.0 # / 10 because thickness provided in mm while X0 in 1/cm 
print "X0 per normal cell: ", x0_per_normal_cell

n_layer = 0
while(x0 < x0_toReach):
   n_layer += 1
   x0 += x0_per_normal_cell

print "After %d layers, we reach %f X0 (at eta = 0)"%(n_layer, x0)

total_calo_thickness = cell_size * n_layer * factor_length_traversed_particle_perpendicular # mm

print "It corresponds to a calo thickness of %f cm"%(total_calo_thickness / 10.0)

total_calo_thickness_withServices = total_calo_thickness + inner_margin_for_services + outer_margin_for_services + cryostat_thickness_in + cryostat_thickness_out # mm

print "Including cryostat and services: it corresponds to a calo thickness of %f cm"%(total_calo_thickness_withServices / 10)

# Compute the total number of layers in phi
cryostat_starting_radius = 2100 # mm

cells_starting_radius = cryostat_starting_radius + cryostat_thickness_in + inner_margin_for_services
print "Readout cells start at %f mm"%cells_starting_radius
circumferance = 2 * pi * cells_starting_radius
cellSize_circle_projection = cell_size * cos(radians(inclination))
print cellSize_circle_projection

n_cells = circumferance / cellSize_circle_projection

print "Need %f cells to cover the 2*pi*r in phi"%n_cells
print "Phi granularity if we read two cells together: %f"%(2*pi/(n_cells/2))
