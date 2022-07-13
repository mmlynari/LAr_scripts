import ROOT

generator = ROOT.TRandom3()

ill_frac = 0.0017
ill_frac = 0.4
n_trial = 100000
n_overlap = 0
for trial in range(n_trial):
    a = generator.Uniform(0 + ill_frac/2., 1 - ill_frac/2.)
    b = generator.Uniform(0 + ill_frac/2., 1 - ill_frac/2.)
    if(abs(a-b) < ill_frac):
        n_overlap += 1
print("Overlap probability = ", n_overlap / float(n_trial))


