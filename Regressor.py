import pandas as pd
import sklearn
import scipy
import sklearn.ensemble as forest
import matplotlib.pyplot as plt
import numpy as np
import sklearn.metrics as metrics
from sklearn.gaussian_process import kernels
import seaborn as sns
import tensorflow as tf
import keras
from keras import optimizers
from sklearn import preprocessing
from sklearn.pipeline import Pipeline, TransformerMixin, FeatureUnion
from sklearn.preprocessing import *
from feature_engine import wrappers, outliers
import optuna
from keras.activations import leaky_relu
from sklearn.pipeline import Pipeline   
from optuna.visualization import _param_importances,plot_optimization_history
# from feature_engine import
from sklearn.compose import ColumnTransformer
from sklearn import model_selection
from sklearn.gaussian_process import GaussianProcessRegressor
from feature_engine.creation import MathFeatures, RelativeFeatures
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
from sklearn.feature_selection import RFECV
import json
import pickle
from scikeras.wrappers import KerasRegressor
from sklearn.ensemble import (
    RandomForestRegressor,
    GradientBoostingRegressor,
    AdaBoostRegressor,
)
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.svm import SVR
from xgboost import XGBRegressor
from sklearn.gaussian_process.kernels import RBF
from lightgbm import LGBMRegressor


class Regressor:
    def __init__(self, Data: pd.DataFrame, target: str):
        self.random_state = 5
        

        self.Data = Data
        self.dof = list((set(self.Data.columns) - set([target])))
        self.X = self.Data[self.dof]
        self.y = self.Data[target]
        # self.cv_splits = self.data_splitter_cv()

        # self.X_train = (self.Data.loc[self.cv_splits[0][0]])[self.dof]
        # self.y_train = (self.Data.loc[self.cv_splits[0][0]])[target]
        try:
            with open(r'models_parameters.pickle','rb') as file:
                self.params=pickle.load(file)
        except: None

        # self.X_test = (self.Data.loc[self.cv_splits[0][1]])[self.dof]
        # self.y_test = (self.Data.loc[self.cv_splits[0][1]])[target]
        # self.base_models = self._define_model()
        

        pass


    def param_optimization(self,X,y):
        optuna_param_dist = {
            "RandomForestRegressor": {
                "n_estimators": optuna.distributions.IntDistribution(100, 1500),
                "min_samples_split": optuna.distributions.IntDistribution(2, 25),
                "min_samples_leaf": optuna.distributions.IntDistribution(2, 32),
                "max_depth":optuna.distributions.IntDistribution(2, 32)
            },
            "GradientBoostingRegressor": {
                "subsample": optuna.distributions.FloatDistribution(0.2, 0.8),
                "n_estimators": optuna.distributions.IntDistribution(100, 1500),
                "min_samples_split": optuna.distributions.IntDistribution(2, 35),
                "min_samples_leaf": optuna.distributions.IntDistribution(2, 35),
                "max_depth": optuna.distributions.IntDistribution(2, 35),
                "learning_rate": optuna.distributions.FloatDistribution(0.01, 0.5),
            },
            "XGBRFRegressor": {
                "subsample": optuna.distributions.FloatDistribution(0.2, 0.8),
                "reg_lambda": optuna.distributions.FloatDistribution(0.2, 0.7),
                "reg_alpha": optuna.distributions.FloatDistribution(0.08, 0.5),
                "max_depth": optuna.distributions.IntDistribution(2, 8),
                "learning_rate": optuna.distributions.FloatDistribution(0.01, 0.5),
                "gamma": optuna.distributions.IntDistribution(1, 3),
                "colsample_bytree": optuna.distributions.FloatDistribution(0.2, 0.8),
            },
            "LGBMRegressor": {
                "subsample": optuna.distributions.FloatDistribution(0.2, 0.8),
                "num_leaves": optuna.distributions.IntDistribution(25, 35),
                "n_estimators": optuna.distributions.IntDistribution(100, 1500),
                "min_child_samples": optuna.distributions.IntDistribution(35, 60),
                "max_depth": optuna.distributions.IntDistribution(2, 10),
                "learning_rate": optuna.distributions.FloatDistribution(0.01, 0.5),
                "colsample_bytree": optuna.distributions.FloatDistribution(0.2, 0.8),
            },
            "RNA": {
                "hidden_layer_units": optuna.distributions.IntDistribution(2, 128),
                "num_leaves": optuna.distributions.IntDistribution(25, 35),
                "n_estimators": optuna.distributions.IntDistribution(150, 250),
                "min_child_samples": optuna.distributions.IntDistribution(35, 60),
                "max_depth": optuna.distributions.IntDistribution(2, 10),
                "learning_rate": optuna.distributions.FloatDistribution(0.01, 0.5),
                "colsample_bytree": optuna.distributions.FloatDistribution(0.2, 0.8),
            },
            "GaussianProcessRegressor": {
                "alpha": optuna.distributions.FloatDistribution(1e-10, 1e-1),
                "kernel":optuna.distributions.CategoricalDistribution([kernels.RBF()])
            },
        }
        optuna_params = {}

        models = {
            "RandomForestRegressor": RandomForestRegressor,
            "GradientBoostingRegressor": GradientBoostingRegressor,
            "XGBRFRegressor": XGBRegressor,
            "LGBMRegressor": LGBMRegressor,
            # "RNA": self.optuna_rna,
            "GaussianProcessRegressor": GaussianProcessRegressor
        }
        for name, model in models.items():
            if name == "RNA":
                None
            else:
                optuna_search = optuna.integration.OptunaSearchCV(
                    estimator=model(),refit=True,
                    param_distributions=optuna_param_dist[name],
                    cv=5,
                    verbose=1,
                    scoring=metrics.make_scorer(metrics.r2_score),
                    n_trials=50,

                )
                
                optuna_search.fit(X, np.ravel(y))
                params = optuna_search.study_.best_params
                optuna_params.update({name: params})
        with open(r"models_parameters.pickle", "wb") as file:
            pickle.dump(optuna_params, file)
        

    
        


    def rna_model_base(self):


        rna_model = keras.Sequential()
        n_layers = 4
        rna_model.add(
                    keras.layers.Dense(
                        20, input_shape=(6,)
                    )
                )
        for i in range(n_layers):
            activation = leaky_relu
            num_hidden = 128
            
            rna_model.add(keras.layers.Dense(num_hidden, activation))
            # dropout = 0.3
            # rna_model.add(keras.layers.Dropout(dropout))
        rna_model.add(keras.layers.Dense(1, activation="linear"))
        # lr = params["lr"]
        # optimizer=trial.suggest_categorical('optimizer',[keras.keras.optimizers.SGD(),keras.keras.optimizers.Adam])
        rna_model.compile(loss="mean_squared_error", optimizer='adam', metrics=["mae"])
        return rna_model

    

    # def data_splitter_cv(self):
    #     groups = pd.qcut(self.y, 100)

    #     k_fold = model_selection.KFold(shuffle=True)
    #     fold_nums = np.zeros(len(self.Data))
    #     for n, (t, v) in enumerate(k_fold.split(groups, groups)):
    #         fold_nums[v] = n

    #     cv_splits = []
    #     for i in range(5):
    #         test_indices = np.argwhere(fold_nums == i).flatten()
    #         train_indices = list(set(range(len(self.Data))) - set(test_indices))
    #         cv_splits.append((train_indices, test_indices))
    #     return cv_splits

    def param_models(self):
        models = {}
        for name, mod in self._define_model().items():
            if name == "RNA":
                m = KerasRegressor(model=mod(),epochs=500,batch_size=10)
                pipe=Pipeline([
                    ('StandardScaler',StandardScaler()),
                    (name,m)
                ])
                models.update({name: pipe})
            else:
                m = mod(**self.params[name])
                pipe=Pipeline([
                    ('StandardScaler',StandardScaler()),
                    (name,m)])
                models.update({name: pipe})

        with open(r"models.pickle", "wb") as pickle_file:
            pickle.dump(models, pickle_file)

    def _define_model(self):
        models = {
            "RandomForestRegressor": RandomForestRegressor,
            "GradientBoostingRegressor": GradientBoostingRegressor,
            "XGBRFRegressor": XGBRegressor,
            "LGBMRegressor": LGBMRegressor,
            "RNA": self.rna_model_base,
            "GaussianProcessRegressor": GaussianProcessRegressor,
        }

        return models
