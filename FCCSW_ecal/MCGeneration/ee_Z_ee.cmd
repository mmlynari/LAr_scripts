! main03.cmnd.
! This file contains commands to be read in for a Pythia8 run.
! Lines not beginning with a letter or digit are comments.
! Names are case-insensitive  -  but spellings-sensitive!
! The settings here are illustrative, not always physics-motivated.

! 1) Settings used in the main program.
Random:setSeed = on  ! set to off to have real random number each time, set to on with fixed see for reproducibility
Random:seed = 54
Main:numberOfEvents = 1000        ! number of events to generate
Main:timesAllowErrors = 5          ! how many aborts before run stops

! 2) Settings related to output in init(), next() and stat().
Init:showChangedSettings = on      ! list changed settings
Init:showChangedParticleData = off ! list changed particle data
Next:numberCount = 10             ! print message every n events
Next:numberShowInfo = 1            ! print event information n times
Next:numberShowProcess = 1         ! print process record n times
Next:numberShowEvent = 0           ! print event record n times

! 3) Beam parameter settings. Values below agree with default ones.
Beams:idA = 11                   ! first beam, e = 2212, pbar = -2212
Beams:idB = -11                   ! second beam, e = 2212, pbar = -2212

! 4) Hard process : Z->qqbar at Ecm=91 GeV
Beams:eCM = 91  ! CM energy of collision

WeakSingleBoson:ffbar2ffbar(s:gmZ) = on 
WeakZ0:gmZmode = 2 ! means no gamma*, only pure Z0
23:onMode = off ! switch off all decay mode from the Z (pdg Id = 23)
23:onIfMatch = 11 -11 ! switch on the decay of Z to e+ e-
ParticleData:modeBreitWigner = 0 ! means we do not have breit wigner, Z0 fixed mass at 91.2 GeV, avoid comvolution of BW with experimental resolution

! PhaseSpace:pTHatMin = 30 ! try to restrict the phase space to the barrel (could not find theta cut --> ask a minimum P_T, since particel energy is ~fixed and they are back to back, it should lead to some constrain on theta. P_T = p * sin(theta) --> 91*sin(20) = ~30 --- theta max = 20 degrees (90 degrees corresponds to eta = 0). Not sure why but this lead to cut at 40 degree (I leave it like that because it is probably good to avoid being on the calo edge)
