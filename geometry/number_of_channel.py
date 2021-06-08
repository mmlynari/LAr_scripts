import math

n_cells = 1536
phi_granularity_degradation_factor = 2
n_eta_segment = 2 * 81 
n_trace_per_eta_layer = 15
n_channels_total = (n_cells * n_eta_segment * n_trace_per_eta_layer)/phi_granularity_degradation_factor
n_channel_per_feedthrough = 20000.0
n_feedthrough = n_channels_total/n_channel_per_feedthrough
feedthrough_diameter = 60.0 # cm
outerCryoBoundaryRadius = 270 # cm feedtrhouh are placed on the outer side, dans une excroissance par rapport au calo lui meme
outerCryoCircumference = 2 * math.pi * outerCryoBoundaryRadius
outerCryoMaxFeedthrough = 2 * outerCryoCircumference / feedthrough_diameter # 2 because we can put them + and - z position

innerCryoBoundaryRadius = 210 # cm 
innerCryoCircumference = 2 * math.pi * innerCryoBoundaryRadius
innerCryoMaxFeedthrough = 2 * innerCryoCircumference / feedthrough_diameter # 2 because we can put them + and - z position

print "Number of channels for the barrel = %d"%n_channels_total
print "Number of feedthrough for the barrel = %d"%n_feedthrough
print "Maximum number of feedthrough on the outer side if they touch each other: 2 * %f/%f = %f"%(outerCryoCircumference, feedthrough_diameter, outerCryoMaxFeedthrough)
print "Maximum number of feedthrough on the inner side if they touch each other: 2 * %f/%f = %f"%(innerCryoCircumference, feedthrough_diameter, innerCryoMaxFeedthrough)
