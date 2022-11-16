#!/usr/bin/env python

import numpy as np
import json

class UpDownStreamCorrector:
    def __init__(self, corr_file):
        with open(corr_file, 'r') as jsonfile:
            source = json.load(jsonfile)
            dict_up = {e['name'] : e['value'] for e in source['corr_params'] if e['type']=='upstream'}
            dict_do = {e['name'] : e['value'] for e in source['corr_params'] if e['type']=='downstream'}
            self.corrs_up = [dict_up['a'], dict_up['b'], dict_up['c'], dict_up['d'], dict_up['e'], dict_up['f']]
            self.corrs_do = [dict_do['a'], dict_do['b'], dict_do['c'], dict_do['d'], dict_do['e'], dict_do['f']]

    def upstream_correction(self, E, E0):
        p0 = self.corrs_up[0]+self.corrs_up[1]/(E-self.corrs_up[2])
        p1 = (self.corrs_up[3]+self.corrs_up[4]/(E-self.corrs_up[5])) * E0
        return p0+p1

    def downstream_correction(self, E, E11):
        p0 = self.corrs_do[0]+self.corrs_do[1]*E
        p1 = (self.corrs_do[2]+self.corrs_do[3]/np.sqrt(E)) * E11
        p2 = (self.corrs_do[4]+self.corrs_do[5]/E) * E11 * E11
        return p0+p1+p2

class MVAUpDownCorrector:
    def __init__(self, up_file, do_file):
        import xgboost as xgb
        self.reg_up = xgb.XGBRegressor(tree_method="hist")
        self.reg_up.load_model(up_file)
        self.reg_do = xgb.XGBRegressor(tree_method="hist")
        self.reg_do.load_model(do_file)

    def upstream_correction(self, Es):
        import xgboost as xgb
        return self.reg_up.predict(Es.T)

    def downstream_correction(self, Es):
        import xgboost as xgb
        return self.reg_do.predict(Es.T)

class LayerCorrector:
    def __init__(self, corr_file):
        with open(corr_file, 'r') as jsonfile:
            self.corrs = json.load(jsonfile)
            self.p0 = []
            self.p1 = []
            self.p2 = []
            self.p3 = []
            for i in range(12):
                self.p0.append(self.corrs[str(i)][0])
                self.p1.append(self.corrs[str(i)][1])
                self.p2.append(self.corrs[str(i)][2])
                self.p3.append(self.corrs[str(i)][3])

    def layer_correction(self, i, E):
        """Correct a single layer"""
        ci = self.corrs[str(i)]
        corr = ci[0] + ci[1]/np.power(E, ci[2]) + ci[3]*E
        return E/corr

    def layers_corrections(self, E):
        """Correct all layers at once
        E: 2D np array [[E0], [E1], [E2],..., [E11]]"""
        corr = self.p0 + np.divide(self.p1, np.power(E.T, self.p2)) + np.multiply(self.p3, E.T)
        return np.divide(E, corr.T)




