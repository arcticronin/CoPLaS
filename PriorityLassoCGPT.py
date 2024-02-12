from sksurv.linear_model import CoxnetSurvivalAnalysis
import pandas as pd
from sklearn.base import BaseEstimator, RegressorMixin

class PriorityLASSO(BaseEstimator, RegressorMixin):
    
    def __init__(self, alpha=1.0, l1_ratio=0.5, max_iter=1000):
        # Initialize hyperparameters for CoxnetSurvivalAnalysis
        self.alpha = alpha
        self.l1_ratio = l1_ratio
        self.max_iter = max_iter

        # Initialize other variables if necessary
        self.models = []

    def fit(self, X, y):
        # X: List of pandas DataFrames with the same number of observations
        # y: Survival data (structured array or pandas DataFrame with event times and censoring indicators)

        # Make a copy of y to update the residual at each iteration
        residual = y.copy()

        for df in X:
            # Create and fit a new CoxnetSurvivalAnalysis model
            model = CoxnetSurvivalAnalysis(alpha=self.alpha, l1_ratio=self.l1_ratio, max_iter=self.max_iter)
            model.fit(df, residual)

            # Store the model in the list
            self.models.append(model)

            # Update the residual for the next iteration
            predicted_survival = model.predict(df)
            residual = y.copy()  # Make a copy to preserve the original data
            residual['event'] = residual['event'] - predicted_survival

        return self

    def predict(self, X):
        # X: List of pandas DataFrames with the same number of observations

        # Initialize the final predictions with zeros
        predictions = pd.DataFrame({'event': [0] * len(X[0])})

        for model in self.models:
            # Use the models to predict survival on each dataframe
            predicted_survival = model.predict(X)
            # Update the final predictions with the predicted_survival
            predictions['event'] += predicted_survival

        return predictions
