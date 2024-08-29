import ROOT
import os, sys, glob, math

# Sampling fraction  
def get_sampling_fraction(th1, prefix, plot_dir_name): 
    canvas_resol = ROOT.TCanvas(prefix + "_resolution", prefix + "_resolution")
    fit_range_min = th1.GetXaxis().GetBinCenter(th1.GetMaximumBin()) - 2 * th1.GetRMS()
    fit_range_max = th1.GetXaxis().GetBinCenter(th1.GetMaximumBin()) + 2 * th1.GetRMS()
    fit_result = th1.Fit("gaus", "SQ", "", fit_range_min, fit_range_max)
    th1.Draw()
    th1.GetXaxis().SetTitle("E_{rec}/E_{tot}")
    canvas_resol.Print(os.path.join(plot_dir_name, prefix + ".png"))
    th1.Write()
    return fit_result

def function(E_Barrel,total_energy_ECal,total_energy_HCal,energy_HCal_first,energy_ECal_last,a,b,c,d):
    E0_Barrel = total_energy_ECal*a + total_energy_HCal*b + c*math.sqrt(abs(energy_ECal_last*a*energy_HCal_first*b)) + d*pow(total_energy_ECal*a,2)
    fitvalue = pow((E_Barrel-E0_Barrel),2)/E_Barrel
    return fitvalue 