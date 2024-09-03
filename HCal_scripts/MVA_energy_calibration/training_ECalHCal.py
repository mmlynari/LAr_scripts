#!/usr/bin/env python

import argparse
import ROOT
import numpy as np
from datetime import datetime


def main():
    """Read command line and trigger processing"""

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputDir", default="FCCSW_ecal/baseline_LAr_testResol_1/",
                        help="Directory containing the input files to process", type=str)
    parser.add_argument("--writeFeatures", default="", help="File containing the values of the input features", type=str)
    parser.add_argument("--writeTarget", default="", help="File containing the values of the output target", type=str)
    parser.add_argument("--training", default=True, action=argparse.BooleanOptionalAction, help='if no-training is set, the training is skipped')
    parser.add_argument("-o", "--outFile", default="out.json",
                        help="Output XGBoost training file", type=str)
    parser.add_argument("--numLayersECalHCal", default=24, type=int)
    parser.add_argument("--numLayersECal", default=11, type=int)
    parser.add_argument("clusters", help="Cluster collection to use", type=str)
    parser.add_argument("--useShapeParameters", default=False, action=argparse.BooleanOptionalAction, help='calculate energies per layers from cells or read them from shapeParameters')
    args = parser.parse_args()
    run(args.inputDir, args.clusters, args.outFile, args.training, args.writeFeatures, args.writeTarget, args.numLayersECalHCal, args.numLayersECal, args.useShapeParameters)


def run(in_directory, clusters, out_file, runTraining, writeFeatures, writeTarget, numLayersECalHCal, numLayersECal, useShapeParameters):
    """Actual processing"""
    init_stuff(useShapeParameters)

    df = ROOT.ROOT.RDataFrame("events", in_directory + "/*.root")
    num_init = df.Count()
    cells = clusters[:-1] + "Cells"
    print("mici mici")
    df = (
        df
        .Define("E_truth_v", "sqrt(genParticles.momentum.y*genParticles.momentum.y+genParticles.momentum.x*genParticles.momentum.x+genParticles.momentum.z*genParticles.momentum.z)")
        .Define("Truth_E", "E_truth_v[0]")
        )
    if useShapeParameters:
        df = df.Alias("clusters_energy", f"Augmented{clusters}.energy")
    else:
        df = (
            df
            .Define(f"{clusters}_EnergyInLayersBoth", f"getCaloCluster_energyInLayersBoth({clusters}, {cells}, {numLayersECalHCal}, {numLayersECal})")
            .Alias("clusters_energy", f"{clusters}.energy")
            )
    df = (
        df
        .Define("lc_idx", "ArgMax(clusters_energy)")
        .Define("Cluster_E", "clusters_energy[lc_idx]")
        )
    
    if useShapeParameters:
        for i in range(numLayersECalHCal):
             df = df.Define(f"Cluster_E{i}", f"_Augmented{clusters}_shapeParameters[lc_idx*33+{i*3}]")
    else:
        for i in range(numLayersECalHCal):
            df = df.Define(f"{clusters}_E{i}", f"getFloatAt({i})({clusters}_EnergyInLayersBoth)")
            df = df.Define(f"Cluster_E{i}", f"{clusters}_E{i}[lc_idx]")

    cols_to_use = ["Truth_E", "Cluster_E"]
    cols_to_use += [f"Cluster_E{i}" for i in range(numLayersECalHCal)]
    v_cols_to_use = ROOT.std.vector('string')(cols_to_use)

    # Filter to remove weird events and get a proper tree
    # d = df.Filter("Cluster_E5!=0 && Cluster_E!=0")
    d = df.Filter("Cluster_E!=0")
    print("We have run on", num_init.GetValue(), "events")

    # Study weird events where clusters don't have cells properly attached
    # df.Filter("Cluster_E5==0").Snapshot("events", "problems.root")

    # Training is so fast it can be done online
    cols = d.AsNumpy(v_cols_to_use)
    
    if useShapeParameters:
        cluster_E = cols["Cluster_E"]
        normalized_layers = np.array([cols[f"Cluster_E{i}"] for i in range(numLayersECalHCal)])
    else:
        layers = np.array([cols[f"Cluster_E{i}"] for i in range(numLayersECalHCal)])
        ## mm: per event sum energy along the rows to get the total cluster ene
        cluster_E = layers.sum(axis=0)
        normalized_layers = np.divide(layers, cluster_E)

    label = cols["Truth_E"] / cluster_E
    ## two alternatives for weights to be tested 
    ## weights = -np.log10(cols["Truth_E"]) + 3
    weights = 1/label**.5

    data = np.vstack([normalized_layers, cluster_E])

    if writeFeatures != "":
        np.save(writeFeatures + "-" + clusters, data)
        np.savetxt(writeFeatures + "-" + clusters, data.transpose())
    if writeTarget != "":
        np.save(writeTarget + "-" + clusters, label)
        np.savetxt(writeTarget + "-" + clusters, label.transpose())

    if runTraining:
        # Somehow importing those libraries earlier in the script does not work. Bad interaction with
        # ROOT libraries in loadGeometry...
        import xgboost as xgb

        # default loss function is reg:squarederror
        # reg = xgb.XGBRegressor(tree_method="approx", objective='reg:squaredlogerror')
        # for squaredlogerror, the corresponding evaluation metric is rmsle: root mean square log error.
        # This metric reduces errors generated by outliers in dataset but
        # since log is used, all input labels are required to be greater than -1.
        reg = xgb.XGBRegressor(tree_method="hist")

        # Even do fancy Hyperparameter optimization for the lulz
        from sklearn.model_selection import RandomizedSearchCV
        import scipy.stats as stats
        param_dist_XGB = {'max_depth': stats.randint(3, 12),  # default 6
                           'n_estimators': stats.randint(300, 800),  # default 100
                           'learning_rate': stats.uniform(0.1, 0.5),  # def 0.3
                           'subsample': stats.uniform(0.5, 1)}
        
        gsearch = RandomizedSearchCV(estimator=reg,
                                      param_distributions=param_dist_XGB,
                                      n_iter=10,
                                      cv=3)
        gsearch.fit(data.T, label, sample_weight=weights)
        print("Best parameters : ", gsearch.best_params_)
        print("Best score (on train dataset CV) : ", gsearch.best_score_)

        ''' 
        # Did HPO once. Now use best fit values for almost instant training
        # reg.objective = 'reg:squaredlogerror'
        # {'learning_rate': 0.11716471738296005, 'max_depth': 6, 'n_estimators': 627, 'subsample': 0.5345527598172761}
        reg.subsample = 0.547
        reg.max_depth = 7
        reg.n_estimators = 773
        reg.learning_rate = 0.29
        # reg.subsample = 0.8
        # reg.max_depth = 4
        # reg.n_estimators = 700
        # reg.learning_rate = 0.32
        print("Start training")
        beg = datetime.now()
        reg.fit(data.T, label, sample_weight=weights, verbose=True)
        # reg.fit(data.T, label)
        end = datetime.now()
        print("Done in ", end - beg)
        reg.save_model(out_file)

        # res = reg.predict(data.T)
        # a = np.array([label, res])
        # print(a)
        ''' 


def init_stuff(useShapeParameters):
    readoutECal = "ECalBarrelModuleThetaMerged"
    readoutHCal = "HCalBarrelReadout"
    readoutName = "HCalBarrelReadout"
    geometryFile = "../../../k4geo/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml"
    ROOT.gROOT.SetBatch(True)
    if not useShapeParameters:
        ROOT.gSystem.Load("libFCCAnalyses")
        _fcc = ROOT.dummyLoader
        ROOT.gInterpreter.Declare("using namespace FCCAnalyses;")
        ROOT.gInterpreter.Declare("using namespace FCCAnalyses::CaloNtupleizer;")
        ROOT.CaloNtupleizer.loadGeometryBoth(geometryFile, readoutName, readoutECal, readoutHCal)
    # ROOT.ROOT.EnableImplicitMT(32)
    ROOT.ROOT.EnableImplicitMT(96)


if __name__ == "__main__":
    main()

    # ROOT.gInterpreter.Declare("""
    # float generatedE(const ROOT::RDF::RSampleInfo &id) {
    # TString s(id.AsString());
    # TPRegexp tpr("_(\\\\d+)\\\\.root");
    # TString ss = s(tpr);
    # TString e_string = ss(1, ss.Index('.')-1);
    # return e_string.Atof() / 1000.;
    # }
    # """)
