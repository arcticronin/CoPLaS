import pandas as pd
import numpy as np

from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.preprocessing import OneHotEncoder

from sklearn import set_config
set_config(display="text")  # displays text representation of estimators

from sklearn.model_selection import GridSearchCV, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from pandas.api.types import CategoricalDtype
 
import warnings
from sklearn.exceptions import FitFailedWarning
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from sklearn.base import BaseEstimator, RegressorMixin

# conider replacing base estimator with CoxnetSurvivalAnalysis

class PriorityLassoSurv(BaseEstimator, RegressorMixin):
    
    def __init__(self, modality_list, modality_order = None):
        # Initialize hyperparameters
        assert(modality_list != None)
        
        self.modality_order = modality_order
        # self.coef_
        
        # Initialize other variables if necessary
        self.model = None

    def fit(self, X, y):
        # X: Input features (NumPy array or pandas DataFrame)
        # y: Survival data (structured array or pandas DataFrame with event times and censoring indicators)

        # Implement your model's training here using X and y
        
        if modality_order == None:
            modality_order = _best_priority_order(X, y, modality_list)  
        ## todo check if modality_order has sense or not?
        else:
            for i in modality_order:
                assert(i < len(modality_list))        
        
        for m in modality_order:
            
            Xt = X.iloc[modality_list[modality_order[m]]]
            # Initialize the 5-fold cross-validator
            kf = KFold(n_splits=5, shuffle=True, random_state=42)
            coefficients_df = pd.DataFrame()
            alphas = []

            ## alpha -> best alpha on average -> refit best alpha each fold and take coefficients
            ## mean mean absolute error  of non-censored !!!!!!!
            residual = y
            
            for i, (train_index, test_index) in enumerate(kf.split(Xt)):
                X_train = Xt.iloc[train_index]
                X_test = Xt.iloc[test_index]
                y_train = y[train_index]
                y_test = y[test_index]

                coxnet_pipe = make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(l1_ratio=0.9, alpha_min_ratio=0.01, max_iter=100))
                warnings.simplefilter("ignore", UserWarning)
                warnings.simplefilter("ignore", FitFailedWarning)
                coxnet_pipe.fit(X_train, y_train)
                estimated_alphas = coxnet_pipe.named_steps["coxnetsurvivalanalysis"].alphas_
                cv = KFold(n_splits=5, shuffle=True, random_state=42)
                gcv = GridSearchCV(
                    make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(l1_ratio=0.9)),
                    param_grid={"coxnetsurvivalanalysis__alphas": [[v] for v in estimated_alphas]},
                    cv=cv,
                    error_score=0.5,
                    n_jobs=1,
                ).fit(X_train, y_train)

                #cv_results = pd.DataFrame(gcv.cv_results_)
                best_model = gcv.best_estimator_.named_steps["coxnetsurvivalanalysis"]
                best_coefs = pd.DataFrame(best_model.coef_, index=X_train.columns, columns=["coefficient"])

                alphas.append(best_model)

                if coefficients_df.empty:
                    coefficients_df = best_coefs
                    coefficients_df.rename(columns={'coefficient': 'coeff_0'}, inplace=True)
                else:
                    coefficients_df[f"coeff_{i}"] = best_coefs
                
                # todo set residual (number of rows)
                # todo save coefficients (columns relative to this modality_order -> m-esima)
    
        # coefficents of a cross validation in a single omics
        coefficents_cv = pd.DataFrame(coefficients_df.mean(axis=1))
        
        #self.model = SomeSurvivalModel()
        self.model = CoxnetSurvivalAnalysis()

        return self

    def predict(self, X):
        # X: New data to predict survival on (NumPy array or pandas DataFrame)

        # Implement your model's prediction here using the trained model
        predictions = self.model.predict(X)

        return predictions

    def get_hazard_ratios(self):
        # If you're implementing a Cox PH model, for example, you could define a method to get hazard ratios.
        return self.model.hazard_ratios_
    
    def _best_priority_order(X, y, modality_list):
        # fit a ridge regression for each omic and take mean, of absolute values of coefficients
        scores = []
        for m in modality_list:
            # get alphas from a quickly fit ridge model
            coxnet_pipe = make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(l1_ratio=0.9, alpha_min_ratio=0.01, max_iter=100))
            
            warnings.simplefilter("ignore", UserWarning)
            warnings.simplefilter("ignore", FitFailedWarning)
            
            # maybe iloc(m)
            coxnet_pipe.fit(X[m], y)
            estimated_alphas = coxnet_pipe.named_steps["coxnetsurvivalanalysis"].alphas_

            # try alphas in the model
            cv = KFold(n_splits=5, shuffle=True, random_state=42)
            gcv = GridSearchCV(
                make_pipeline(StandardScaler(), CoxnetSurvivalAnalysis(l1_ratio=0.01)),
                param_grid={"coxnetsurvivalanalysis__alphas": [[v] for v in estimated_alphas]},
                cv=cv,
                error_score=0.5,
                n_jobs=1,
            ).fit(X[m], y)
            cv_results = pd.DataFrame(gcv.cv_results_)
            
            # get best model and append the score to the list
            best_model = gcv.best_estimator_.named_steps["coxnetsurvivalanalysis"]
            scores.append(np.absolute(best_model.coef_).mean())
        return np.argsort(scores)
        
