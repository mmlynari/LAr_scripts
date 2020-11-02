n_layer = 12
n_cell_strip = 4
layers_extracted_front = 4

n_trace_front = layers_extracted_front - 2 + 4
n_trace_back = n_layer - layers_extracted_front - 1

print "Extracting the first %d layers from front."%layers_extracted_front
print "Number of trace in the front %d, ignoring the fact that we have to extract signal from first layer as well"%n_trace_front
print "Number of trace in the back %d, ignoring the fact that we have to extract signal from last layer as well"%n_trace_back


