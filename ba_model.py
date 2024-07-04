import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.neighbors import KNeighborsRegressor

# Load the dataset
df = pd.read_csv('./network_features/features_kd_243_with_networks.csv')

# Extract features and target variable
x = df[['nonpolar_polar', 'mcsc_hbond', 'Backbone Hbond', 'Solvation Polar', 'Interface Residues VdW Clashing', 'ele_e']].values
y = df['dG'].values

# Split the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42)

# Define the KNN regressor
model = KNeighborsRegressor()

# Define the grid of hyperparameters to search
param_grid = {
    'n_neighbors': [2,3,4, 5,6,7,8, 9,10, 11,12],  # Number of neighbors
    'weights': ['uniform', 'distance'],  # Weight function used in prediction
    'p': [1, 2]  # Power parameter for the Minkowski metric
}

# Perform Grid Search to find the best combination of hyperparameters
grid_search = GridSearchCV(estimator=model, param_grid=param_grid, cv=5, scoring='neg_mean_absolute_error')
grid_search.fit(X_train, y_train)

# Get the best hyperparameters
best_params = grid_search.best_params_
print("Best Hyperparameters:", best_params)

# Use the best model for prediction
best_model = grid_search.best_estimator_

# Predict on training set and test set
y_train_pred = best_model.predict(X_train)
y_test_pred = best_model.predict(X_test)

# Evaluate model performance
mae_train = mean_absolute_error(y_train, y_train_pred)
mae_test = mean_absolute_error(y_test, y_test_pred)
r2_train = r2_score(y_train, y_train_pred)
r2_test = r2_score(y_test, y_test_pred)
pearson_corr_train, _ = pearsonr(y_train, y_train_pred)
pearson_corr_test, _ = pearsonr(y_test, y_test_pred)

print("MAE Train:", mae_train)
print("MAE Test:", mae_test)
print("R^2 Train:", r2_train)
print("R^2 Test:", r2_test)
print("r_train:", pearson_corr_train)
print("r_test:", pearson_corr_test)

# Plot actual vs predicted values
plt.scatter(y_train, y_train_pred, color='blue', label='Training data')
plt.scatter(y_test, y_test_pred, color='red', label='Testing data')
plt.xlabel('Actual ΔG (kcal/mol)')
plt.ylabel('Predicted ΔG (kcal/mol)')
plt.legend()
plt.show()

