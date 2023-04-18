ddsim --compactFile $PWD/CLD_LAr/FCCee_o1_v05.xml \
      --enableGun \
      --gun.distribution uniform \
      --gun.energy "10*GeV" \
      --gun.particle e- \
      --numberOfEvents 10 \
      --outputFile Step1_edm4hep.root
