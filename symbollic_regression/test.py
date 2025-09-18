import numpy as np
import pandas as pd
from pysr import PySRRegressor

# Create synthetic data
np.random.seed(0)
X = np.random.randn(100, 2)  # 100 samples, 2 features
y = X[:, 0] * np.sin(X[:, 1]) + X[:, 1]**2  # True relationship

# Define the symbolic regressor
model = PySRRegressor(
    niterations=40,               # Number of evolutionary iterations
    binary_operators=["+", "-", "*", "/"],
    unary_operators=["sin", "cos", "exp", "log"],
    model_selection="best",      # Choose the best model by loss
    loss="loss(x, y) = (x - y)^2",  # Mean squared error
    verbosity=1,
)

# Fit the model
model.fit(X, y)

# Show best equation
print(model)

# Predict using symbolic model
y_pred = model.predict(X)
