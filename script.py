import pandas as pd
import os
import numpy as np
from pycox.evaluation.concordance import concordance_td
from pycox.evaluation import EvalSurv

##BLCA/intermediate_mean/clinical_gex

def gen_metrics( base_folder, ):
    
    for i in range(1,26):
        surv_f_csv = os.join(base_folder, "survival_funcion", "TCGA" ,cancer_type, fusion, modalities, f"split_{i}.csv")
    
    metrics_csv =  f"metrics/BLCA/intermediate_mean/clinical_gex/metrics.csv"
    data_csv = "processed/BLCA_data_preprocessed.csv"
    splits_test_csv = "splits/BLCA_test_splits.csv"
    splits_train_csv = "splits/BLCA_train_splits.csv"
    survival_funct_path = "survival_functions/BLCA/late_mean/clinical_gex/" ## + split_1.csv etc
    
    X = pd.read_csv(data_csv)
    
    #split = pd.read_csv(pthsf + "split_1.csv")
    splits_test = pd.read_csv(splits_test_csv)
    splits_test.shape

   

    event_indicator, event_time = X["OS"].to_numpy(), X["OS_days"].to_numpy()


    pd.read_csv(survival_funct_path + f"split_1.csv").to_numpy().shape

    l = []
    split = pd.read_csv(splits_test_csv)

    for idx in range(25):
        surv = pd.read_csv(survival_funct_path + f"split_{idx+1}.csv").to_numpy().T
        test_idx = split.iloc[idx, :].tolist()

        cidx = concordance_td(
            method='adj_antolini',
            durations = event_time[test_idx],
            events = event_indicator[test_idx],
            surv = surv,
            surv_idx= np.searchsorted(np.unique(test_idx), 
                                      test_idx)
        )
        l.append(cidx)
        # unique già fa il sorting