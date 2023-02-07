# Scripts for electronics studies

- `analyse_scope_output.py` script to read the scope data, plot them, apply shaper on top of it and derive peak-to-peak cross talk
- `analyse_scope_output_and_compare.py` same as above, but it additionally compares the cross-talk values with the simulated ones
- `apply_s_parameters.py` old script to apply S parameters on the triangular signal to get various output signals
- `illuminated_detector_fraction.py` computes the probablility to have hit overlapping in consecutive events
- `pileup_at_FCC_ee.py` study how out of time pile up affects us
- `prepare_csv_for_ANSYS_from_scope.py` script used once to write the input signal to feed ANSYS, directly taken from the scope output
- `signal_csv.py` script to write the triangular signal in csv format (useful to feed various tools)
- `analyse_ansys_data.py` take the output from ANSYS (plot exported data), apply shaping, derive cross-talk, ...
- `yparam_crosstalk.py` tentative to derive the cross-talk from Y parameters
