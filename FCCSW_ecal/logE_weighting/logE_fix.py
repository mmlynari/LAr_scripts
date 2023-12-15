import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib as mpl
import math as ma

plt.style.use('bmh')

mpl.rc('font',family='A-OTF Shin Go Pro')

do_resolution_plots = True

samples = {
    '2': ('output_Scurve_prod_evts_5000_pdg_11_2_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 5000),
    '7': ('output_Scurve_prod_evts_5000_pdg_11_7_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 5000),
    '12': ('output_Scurve_prod_evts_5000_pdg_11_12_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 5000),
    '17': ('output_Scurve_prod_evts_5000_pdg_11_17_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 5000),
    '22': ('output_Scurve_prod_evts_5000_pdg_11_22_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 5000),
    '25': ('output_Scurve_prod_evts_1000_pdg_11_25_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '32': ('output_Scurve_prod_evts_5000_pdg_11_32_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 5000),
    '35': ('output_Scurve_prod_evts_1000_pdg_11_35_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '37': ('output_Scurve_prod_evts_5000_pdg_11_37_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 5000),
    '42': ('output_Scurve_prod_evts_1000_pdg_11_42_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '45': ('output_Scurve_prod_evts_1000_pdg_11_45_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '47': ('output_Scurve_prod_evts_1000_pdg_11_47_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '52': ('output_Scurve_prod_evts_1000_pdg_11_52_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '55': ('output_Scurve_prod_evts_1000_pdg_11_55_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '57': ('output_Scurve_prod_evts_1000_pdg_11_57_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '62': ('output_Scurve_prod_evts_1000_pdg_11_62_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '65': ('output_Scurve_prod_evts_1000_pdg_11_65_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '67': ('output_Scurve_prod_evts_1000_pdg_11_67_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '72': ('output_Scurve_prod_evts_1000_pdg_11_72_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '75': ('output_Scurve_prod_evts_1000_pdg_11_75_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '77': ('output_Scurve_prod_evts_1000_pdg_11_77_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '82': ('output_Scurve_prod_evts_1000_pdg_11_82_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '85': ('output_Scurve_prod_evts_1000_pdg_11_85_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '87': ('output_Scurve_prod_evts_1000_pdg_11_87_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '92': ('output_Scurve_prod_evts_1000_pdg_11_92_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '95': ('output_Scurve_prod_evts_1000_pdg_11_95_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '97': ('output_Scurve_prod_evts_1000_pdg_11_97_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 1000),
    '0.1': ('output_Scurve_prod_evts_5000_pdg_11_0.1_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 5000),
    '0.5': ('output_Scurve_prod_evts_5000_pdg_11_0.5_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root', 5000),
}



inDir = "/home/almaloiz/mnt/cernbox/thesis/calo/root_Tong/scan_E/"

# w0_list = np.linspace(4,10,9)
w0_list = np.linspace(2,10, 33)

nocorr_resolution = {} # Resolution obtained with linear E weighting
min_w0 = [] # w0 value for which the resolution is minimum for X GeV
energies = []
resolutions = {}
for sample_energy in samples:
    file_name, Nevts = samples[sample_energy]   

# file_path = "./root_merge/output_evts_"+str(Nevts)+"_pdg_"+str(pdgCode)+"_"+str(momentum)+"_GeV"+"_ThetaMinMax_"+str(thetaMin)+"_"+str(thetaMax)+"_PhiMinMax_"+str(phiMin)+"_"+str(phiMax)+"_MagneticField_"+str(magneticField)+"_Noise"+str(addNoise)+".root"

# output_Scurve_prod_evts_40000_pdg_11_50_GeV_ThetaMinMax_40_140_PhiMinMax_-0.2_0.2.root

    print("Doing "+sample_energy+" GeV sample")
    print("\topening file")
    file = uproot.open(inDir+file_name)
    print("\tgetting events...")
    events = file.get('events')

    print("\t\t Particle x,y,z and reponses")
    part_x = events.get('genParticles.momentum.x').array()
    part_y = events.get('genParticles.momentum.y').array()
    part_z = events.get('genParticles.momentum.z').array()

    part_x = np.array([part_x[i][0] for i in range(0, Nevts)])
    part_y = np.array([part_y[i][0] for i in range(0, Nevts)])
    part_z = np.array([part_z[i][0] for i in range(0, Nevts)])
    part_response = np.arctan2(np.sqrt(part_x*part_x + part_y*part_y), part_z)

    # topo_x = np.array([events.get('CorrectedCaloTopoClusters.position.x').array()[i][0] for i in range(0, Nevts)])
    # topo_y = np.array([events.get('CorrectedCaloTopoClusters.position.y').array()[i][0] for i in range(0, Nevts)])
    # topo_z = np.array([events.get('CorrectedCaloTopoClusters.position.z').array()[i][0] for i in range(0, Nevts)])
    # topo_response = np.arctan2(np.sqrt(topo_x*topo_x + topo_y*topo_y), topo_z)

    # resol = np.subtract(part_response, topo_response)

    print("\t\t EcalBarrel cells x,y,z and energy")
    cells_x = events.get('ECalBarrelPositionedCells.position.x').array()
    cells_y = events.get('ECalBarrelPositionedCells.position.y').array()
    cells_z = events.get('ECalBarrelPositionedCells.position.z').array()
    cells_energy = events.get('ECalBarrelPositionedCells.energy').array()
    E_cluster = events.get('CorrectedCaloTopoClusters.energy').array()
    E_cluster = np.array([np.max(E_cluster[i]) for i in range(0, Nevts)])

    #Calculating resolution and response without correction and with simple log(E) linear correction
    #No correction
    weight = cells_energy/E_cluster
        
    posX = np.sum(cells_x * weight, axis = 1)
    posY = np.sum(cells_y * weight, axis = 1)
    posZ = np.sum(cells_z * weight, axis = 1)
    
    weight = np.sum(weight, axis = 1)

    posX = posX / weight
    posY = posY / weight
    posZ = posZ / weight
    
    topo_response = np.arctan2(np.sqrt(posX*posX + posY*posY), posZ)
    nocorr_response = part_response-topo_response
    nocorr_resolution[sample_energy] = np.std(nocorr_response)*1000
    
    #linear log(E) correction
    weight = np.log(cells_energy)
        
    posX = np.sum(cells_x * weight, axis = 1)
    posY = np.sum(cells_y * weight, axis = 1)
    posZ = np.sum(cells_z * weight, axis = 1)
    
    weight = np.sum(weight, axis = 1)

    posX = posX / weight
    posY = posY / weight
    posZ = posZ / weight
    
    topo_response = np.arctan2(np.sqrt(posX*posX + posY*posY), posZ)
    linlog_response = part_response-topo_response
    linlog_resolution = np.std(linlog_response)*1000
    
    print(f"Linear log correction result : theta_resolution = {linlog_resolution}")
    
    def plot_response(resp, energy, w0):
        plt.figure()
        # plt.grid()
        plt.suptitle(f"theta response for {energy} GeV samples",
                    fontsize=15,
                    fontweight='bold')
        plt.title(f"W = max(0, {w0} + ln(Ecell/Ecluster))",
                fontsize=9,
                loc='left',
                style='italic')
        plt.hist(resp, 100, label=f"w0={w0}", alpha=0.7)
        plt.hist(nocorr_response, 100, alpha=0.3, label="no corr.")
        plt.xlabel('theta_{e}-theta_{cluster}')
        plt.legend()
        plt.show()
        plt.savefig(f"plots/responses/response_{sample_energy}GeV_w0_{w0}.png")
        plt.close()
    

    
    w0_responses = []
    w0_resolutions = []
    min_reso = 9999 # will store the w0 value for which the resol is minimal
    _w0 = 9999
    for w0 in w0_list:
        print(f"\tDoing w0 = {w0}")
        
        weight =  np.maximum(0, w0 + np.log(cells_energy/E_cluster))
        
        posX = np.sum(cells_x * weight, axis = 1)
        posY = np.sum(cells_y * weight, axis = 1)
        posZ = np.sum(cells_z * weight, axis = 1)
        
        weight = np.sum(weight, axis = 1)
 
        posX = posX / weight
        posY = posY / weight
        posZ = posZ / weight
        
        topo_response = np.arctan2(np.sqrt(posX*posX + posY*posY), posZ)
        w0_resp = part_response - topo_response
        w0_responses.append(w0_resp)
        
        plot_response(w0_resp, sample_energy, w0)
        reso = np.std(w0_resp)*1000
        w0_resolutions.append(reso)
        if reso < min_reso and not ma.isnan(reso):
            min_reso = reso
            _w0 = w0
        
    # w0_resolutions = [np.std(resp)*1000 for resp in w0_responses]
    print("Resolutions for " + sample_energy +" GeV are : ")
    print(w0_resolutions) 
    min_w0.append(_w0)
    energies.append(float(sample_energy))
    resolutions[sample_energy] = w0_resolutions




Nsamples = len(samples)
Nplots = Nsamples//3 + Nsamples%3
Nline = Nplots//3 + 1 if Nplots%3 != 0 else 0

fig, axes = plt.subplots(Nline, 3, tight_layout=True, figsize=(18,16), sharey='row')
# fig.subplots_adjust(bottom=0.2)
fig.suptitle("theta resolution with weighting",
    fontsize=13,
    fontweight='bold',
    )

i_plot = 0
colors = ["#348ABD", "#A60628", "#7A68A6"]
# colors = ["#8dd3c7", "#feffb3", "#bfbbd9"]
for s in samples:
    plt.subplot(Nline,3,i_plot//3+1)
    
    plt.plot([],[],' ', label=f"{s} GeV") #hack to have a beautiful legend
    plt.plot(w0_list, resolutions[s], label=f"corr.", color=colors[i_plot%3]) # reminder that resolutions are multiplied by 1000 to ease readability on the axis
    plt.plot([w0_list[0],w0_list[-1]], [nocorr_resolution[s], nocorr_resolution[s]], '--', label='no corr.', alpha = 0.5, color=colors[i_plot%3])
    plt.subplots_adjust(bottom=0.4)
    if i_plot%3 == 2:
        plt.title(f"theta res. for "+ f"{energies[i_plot-2]}, {energies[i_plot-1]}, {energies[i_plot]} GeV", style='italic', fontsize = 9, loc='left')
        plt.legend(ncol=3, fontsize=11, labelspacing=0.05, loc='upper right')
        plt.ylabel('theta resolution*1000')
        plt.xlabel('w0', loc='right')
    if i_plot == Nsamples-1:
        plt.title(f"theta res. for " + (f"{energies[i_plot-1]}, {energies[i_plot]} GeV" if i_plot%3 == 1 else f"{energies[i_plot]} GeV"), style='italic', fontsize = 9, loc='left')
        plt.legend(ncol=3-i_plot%3, fontsize=11, labelspacing=0.05)
        plt.ylabel('theta resolution*1000')
        plt.xlabel('w0', loc='right')
    i_plot+=1
    
plt.show()  
plt.savefig('plots/resolutions_w0.png')
plt.close()

plt.figure()
plt.suptitle("w0_min value with the energy",
    fontsize=13,
    fontweight='bold',
    )
plt.title("w0_min -> w0 such as theta_res(E) is the smallest", style='italic', fontsize = 9, loc='left')
plt.plot(energies, min_w0, 'x')
plt.xlabel('E part. (GeV)', loc='right')
plt.ylabel('min w')
plt.show()
plt.savefig('plots/w0_min_energy')

print("\n min w0 : ")
print(min_w0)
print("Energies : ")
print(energies)


#events->Draw("atan2(sqrt(CorrectedCaloTopoClusters.position.x*CorrectedCaloTopoClusters.position.x+CorrectedCaloTopoClusters.position.y*CorrectedCaloTopoClusters.position.y),CorrectedCaloTopoClusters.position.z)-atan2(sqrt(genParticles.momentum.x*genParticles.momentum.x+genParticles.momentum.y*genParticles.momentum.y),genParticles.momentum.z)>>h2")

