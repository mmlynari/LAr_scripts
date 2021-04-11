lar_inner = 1.24 # mm
lar_outer = 2.4 # mm
height = 565 # mm
n_plates = 1536
detector_width = 2200 * 2 # mm 
volume_signal_gap = (lar_inner + lar_outer) * height * detector_width / 2.0 # cubic milimeter
volume_signal_gap_cubic_meter = volume_signal_gap * 1e-9
barrel_sensitive_volume = volume_signal_gap_cubic_meter * 2 * n_plates # *2 because two gap for 1 absorber
print "Sensitive LAr volume in cubic meter: %f"%barrel_sensitive_volume
