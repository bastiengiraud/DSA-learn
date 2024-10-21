import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def remove_duplicates_with_tolerance(df, tol=0.01):
    # Create a boolean mask to mark rows to keep
    rows_to_keep = [True] * len(df)
    duplicates_indices = []

    for i in range(len(df)):
        if not rows_to_keep[i]:
            continue
        
        # Compare the current row with all subsequent rows
        for j in range(i + 1, len(df)):
            # Check if all values between rows are within the tolerance
            if np.all(np.abs(df.iloc[i] - df.iloc[j]) <= tol):
                rows_to_keep[j] = False
                duplicates_indices.append(j)
    
    # Filter the dataframe based on rows_to_keep
    df_cleaned = df[rows_to_keep].reset_index(drop=True)

    return df_cleaned, duplicates_indices

def find_and_remove_duplicates_with_tolerance(df, tol=0.01):
    # Round the values to the nearest 0.01 to account for tolerance
    rounded_df = df.apply(lambda x: np.round(x, decimals=int(-np.log10(tol))))
    
    # Find the indices of duplicates (keeping the first occurrence)
    duplicates = rounded_df.duplicated(keep='first')
    
    # Get the indices of duplicate rows
    duplicate_indices = df.index[duplicates].tolist()
    
    # Remove duplicates from the original dataframe
    df_cleaned = df.drop(duplicate_indices)
    
    return df_cleaned, duplicate_indices

# define stability boundary
stability_boundary = 0.03
margin = 0.0025
lower_bound = stability_boundary - margin
upper_bound = stability_boundary + margin

# Specify the directory and file name
directory = "C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/output/case39/datasets/"
flow_name_method = '39bus_method_flows.csv' # DT method data
flow_name_lhc = 'lhc_flows.csv'
flow_name_imp = 'imp_flows.csv' 

file_name_method = '39bus_method_ops.csv' # DT method data
file_name_lhc = 'lhc_ops.csv'
file_name_imp = 'imp_ops.csv' 

file_name_analysis = 'data_analysis2.csv'

# Create the full file path
flow_path_method = os.path.join(directory, flow_name_method)
flow_path_lhc = os.path.join(directory, flow_name_lhc)
flow_path_imp = os.path.join(directory, flow_name_imp)

file_path_method = os.path.join(directory, file_name_method)
file_path_lhc = os.path.join(directory, file_name_lhc)
file_path_imp = os.path.join(directory, file_name_imp)

file_path_analysis = os.path.join(directory, file_name_analysis)

# Load the dataset
flow_data_method = pd.read_csv(flow_path_method, sep = ';')
flow_data_lhc = pd.read_csv(flow_path_lhc, sep = ';')
flow_data_imp = pd.read_csv(flow_path_imp, sep = ';')

op_data_method = pd.read_csv(file_path_method, sep = ';')
op_data_lhc = pd.read_csv(file_path_lhc, sep = ';')
op_data_imp = pd.read_csv(file_path_imp, sep = ';')

# remove duplicates
method = op_data_method.drop(columns = ['N0', 'flow_viol', 'over_volt', 'under_volt', 'N1', 'damping', 'distance'], axis=1)
lhc = op_data_lhc.drop(columns = ['N0', 'flow_viol', 'over_volt', 'under_volt', 'N1', 'damping', 'distance'], axis=1)
imp = op_data_imp.drop(columns = ['N0', 'flow_viol', 'over_volt', 'under_volt', 'N1', 'damping', 'distance'], axis=1)

_, method_duplicates = find_and_remove_duplicates_with_tolerance(method, tol=0.01)
_, lhc_duplicates = find_and_remove_duplicates_with_tolerance(lhc, tol=0.01)
_, imp_duplicates = find_and_remove_duplicates_with_tolerance(imp, tol=0.01)

flow_data_method.drop(method_duplicates, inplace=True)#.reset_index(drop=True)
flow_data_lhc.drop(lhc_duplicates, inplace=True)#.reset_index(drop=True)
flow_data_imp.drop(imp_duplicates, inplace=True)#.reset_index(drop=True)

op_data_method.drop(method_duplicates, inplace=True)#.reset_index(drop=True)
op_data_lhc.drop(lhc_duplicates, inplace=True)#.reset_index(drop=True)
op_data_imp.drop(imp_duplicates, inplace=True)#.reset_index(drop=True)

# for every set
# number of samples
method_samples = flow_data_method.shape[0]
lhc_samples = flow_data_lhc.shape[0]
imp_samples = flow_data_imp.shape[0]

print("number of method samples: ", method_samples)
print("number of lhc samples: ", lhc_samples)
print("number of importance samples: ", imp_samples)
print(" ")

# number of AC feasible samples
feasible_method = flow_data_method['feasible'].value_counts()[1]
feasible_method_percentage = (feasible_method/method_samples)*100
feasible_lhc = flow_data_lhc['feasible'].value_counts()[1]
feasible_lhc_percentage = (feasible_lhc/lhc_samples)*100
feasible_imp = flow_data_imp['feasible'].value_counts()[1]
feasible_imp_percentage = (feasible_imp/imp_samples)*100

print("number of AC feasible method samples: ", feasible_method)
print("number of AC feasible lhc samples: ", feasible_lhc)
print("number of AC feasible importance samples: ", feasible_imp)
print(" ")

# number of small signal stable samples
stable_method = flow_data_method['stable'].value_counts()[1]
stable_method_percentage = (stable_method/method_samples)*100
stable_lhc = flow_data_lhc['stable'].value_counts()[1]
stable_lhc_percentage = (stable_lhc/lhc_samples)*100
stable_imp = flow_data_imp['stable'].value_counts()[1]
stable_imp_percentage = (stable_imp/imp_samples)*100

print("number of small signal stable method samples: ", stable_method)
print("number of small signal stable lhc samples: ", stable_lhc)
print("number of small signal stable importance samples: ", stable_imp)
print(" ")

# number of feasible and small signal stable samples
secure_method = flow_data_method[(flow_data_method['feasible'] == 1) & (flow_data_method['stable'] == 1)].shape[0]
secure_method_percentage = (secure_method/method_samples)*100
secure_lhc = flow_data_lhc[(flow_data_lhc['feasible'] == 1) & (flow_data_lhc['stable'] == 1)].shape[0]
secure_lhc_percentage = (secure_lhc/lhc_samples)*100
secure_imp = flow_data_imp[(flow_data_imp['feasible'] == 1) & (flow_data_imp['stable'] == 1)].shape[0]
secure_imp_percentage = (secure_imp/imp_samples)*100

print("number of secure method samples: ", secure_method)
print("number of secure lhc samples: ", secure_lhc)
print("number of secure importance samples: ", secure_imp)
print(" ")

# number of operating points in the stability boundary
boundary_op_method = op_data_method[(op_data_method['damping'] > lower_bound) & (op_data_method['damping'] < upper_bound)].shape[0]
boundary_method_percentage = (boundary_op_method/method_samples)*100
boundary_op_lhc = op_data_lhc[(op_data_lhc['damping'] > lower_bound) & (op_data_lhc['damping'] < upper_bound)].shape[0]
boundary_lhc_percentage = (boundary_op_lhc/lhc_samples)*100
boundary_op_imp = op_data_imp[(op_data_imp['damping'] > lower_bound) & (op_data_imp['damping'] < upper_bound)].shape[0]
boundary_imp_percentage = (boundary_op_imp/imp_samples)*100

print("number of method samples in boundary region: ", boundary_op_method)
print("number of lhc samples in boundary region: ", boundary_op_lhc)
print("number of importance samples in boundary region: ", boundary_op_imp)
print(" ")

data = {
    'method': ['method', 'lhc', 'importance'],
    'number of samples': [method_samples, lhc_samples, imp_samples],
    'AC feasible samples': [feasible_method, feasible_lhc, feasible_imp],
    'AC feasible percentage [%]': [feasible_method_percentage, feasible_lhc_percentage, feasible_imp_percentage],
    'stable samples': [stable_method, stable_lhc, stable_imp],
    'stable samples percentage [%]': [stable_method_percentage, stable_lhc_percentage, stable_imp_percentage],
    'secure samples': [secure_method, secure_lhc, secure_imp],
    'secure samples percentage [%]': [secure_method_percentage, secure_lhc_percentage, secure_imp_percentage],
    'samples in HIC region': [boundary_op_method, boundary_op_lhc, boundary_op_imp],
    'samples in HIC region percentage [%]': [boundary_method_percentage, boundary_lhc_percentage, boundary_imp_percentage],
}
    
df_analysis = pd.DataFrame.from_dict(data)
df_analysis.to_csv(file_path_analysis, sep = ';')





