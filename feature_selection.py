import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor
import xgboost as xgb
# import lightgbm as lgb
from scipy.stats import pearsonr
from itertools import combinations
from multiprocessing import Pool
from tqdm import tqdm
import warnings

# Ignore all warnings
warnings.filterwarnings("ignore")

# Load your data
data = pd.read_csv('features_kd_243.csv')

# Drop rows with missing values
data = data.dropna()

# Separate features and target
X = data.drop(['dG', 'kd (nM)', 'pdb', 'Antibody name', 'Antigen name', 'Disease', 'disease_type'], axis=1)
y = data['dG']

# Generate combinations of features
def generate_feature_combinations(features, max_combinations):
    return list(combinations(features, max_combinations))

# Define a function to process each combination
def process_combination(combination):
    X_subset = X[list(combination)]
    X_train, X_test, y_train, y_test = train_test_split(X_subset, y, test_size=0.2, random_state=42)
    models = {
        "Linear Regression": LinearRegression(),
        "Ridge Regression": Ridge(),
        "Lasso Regression": Lasso(),
        "ElasticNet Regression": ElasticNet(),
        "Decision Tree": DecisionTreeRegressor(),
        "Random Forest": RandomForestRegressor(),
        "Gradient Boosting": GradientBoostingRegressor(),
        "Support Vector Machine": SVR(),
        "K-Nearest Neighbors": KNeighborsRegressor(),
        "Neural Network": MLPRegressor(),
        "XGBoost": xgb.XGBRegressor(),
#         "LightGBM": lgb.LGBMRegressor()
    }
    results = []
    for name, model in models.items():
        model.fit(X_train, y_train)
        y_train_pred = model.predict(X_train)  # Predict on training set
        y_test_pred = model.predict(X_test)    # Predict on test set
        
        # Evaluate on training set
        mae_train = mean_absolute_error(y_train, y_train_pred)
        mse_train = mean_squared_error(y_train, y_train_pred)
        r2_train = r2_score(y_train, y_train_pred)
        pearson_corr_train, _ = pearsonr(y_train, y_train_pred)
        
        # Evaluate on test set
        mae_test = mean_absolute_error(y_test, y_test_pred)
        mse_test = mean_squared_error(y_test, y_test_pred)
        r2_test = r2_score(y_test, y_test_pred)
        pearson_corr_test, _ = pearsonr(y_test, y_test_pred)
        
        results.append({
            "Model": name,
            "MAE Train": mae_train, "MAE Test": mae_test,
            "MSE Train": mse_train, "MSE Test": mse_test,
            "R^2 Train": r2_train, "R^2 Test": r2_test,
            "Correlation Train": pearson_corr_train,
            "Correlation Test": pearson_corr_test
        })
    return combination, results

# Define the maximum number of features to consider in a combination
max_features = 5
feature_combinations = generate_feature_combinations(X.columns, max_features)

# Process combinations in parallel
if __name__ == '__main__':
    with Pool() as pool:
        results = list(tqdm(pool.imap(process_combination, feature_combinations), total=len(feature_combinations)))

# Collect the results into a DataFrame
all_results = []
for combination, model_results in results:
    for result in model_results:
        all_results.append({**{"Features": combination}, **result})

results_df = pd.DataFrame(all_results)

# sorted_df = results_df.sort_values('Correlation', ascending=False).head(100)

# Save results to a CSV file
results_df.to_csv('regression_results_5comb.csv', index=False)

