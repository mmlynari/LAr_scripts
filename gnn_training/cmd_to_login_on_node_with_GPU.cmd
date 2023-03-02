executable              = dummy.sh
arguments               = $(ClusterId)$(ProcId)
output                  = traiining.$(ClusterId).$(ProcId).out
error                   = training.$(ClusterId).$(ProcId).err
log                     = training.$(ClusterId).log
should_transfer_files   = YES
#transfer_input_files    = matrix.py
when_to_transfer_output = ON_EXIT
request_GPUs = 1
request_CPUs = 1
+JobFlavour = "longlunch"
queue
