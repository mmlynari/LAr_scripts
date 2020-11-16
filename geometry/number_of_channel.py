n_cells = 1536
phi_granularity_degradation_factor = 2
n_eta_segment = 2 * 81 
n_trace_per_eta_layer = 13
n_channels_total = (n_cells * n_eta_segment * n_trace_per_eta_layer)/phi_granularity_degradation_factor
print "Number of channels for the barrel = %d"%n_channels_total
