import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import classification_report, accuracy_score, f1_score, confusion_matrix, precision_score, recall_score, det_curve
from sklearn import tree
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.model_selection import GridSearchCV
plt.style.use(['C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/0) DTU Admin/5) Templates/thesis.mplstyle'])
plt.rcParams['text.usetex'] = False


def fpr_score(y_actual, y_hat):
    """
    Function to determine the FP score, and to get the indices of the OPs
    where there is missclassification.

    y_actual = pandas series, comes from train-test split
    y_hat = np.array, comes from prediction

    """
    TP = 0
    FP = 0
    TN = 0
    FN = 0

    FP_list = []
    FN_list = []

    counter = 0
    for i, val in y_actual.items():
        if y_actual[i]==y_hat[counter]==1:
            TP += 1
        if y_hat[counter]==1 and y_actual[i]!=y_hat[counter]:
            FP += 1
            FP_list.append(i)
        if y_actual[i]==y_hat[counter]==0:
            TN += 1
        if y_hat[counter]==0 and y_actual[i]!=y_hat[counter]:
            FN += 1
            FN_list.append(i)
        counter += 1

    try:
        FPR = FP / (FP + TN)
    except:
        FPR = 33

    try:
        FNR = FN / (TP + FN)
    except:
        FNR = 33
  

    return FPR, FNR, FP_list, FN_list

def plot_damping_histograms(x_lb, bar_width, damping_data, training_data, sec_training_data, insec_training_data, title):
    """
    Function to plot the damping histograms for training and testing data.
    """
    # Define bin edges with steps of 'bar_width' starting from 'x_lb'
    bins = np.arange(x_lb, max(damping_data.max(), sec_training_data.max()) + bar_width, bar_width)

    # Create a figure with two subplots, arranged horizontally (1 row, 2 columns)
    fig, axs = plt.subplots(1, 4, figsize=(14, 6))

    max_y = 0  # Track max y-value to set consistent y-limits

    # Plot the histogram for 'training_data'
    axs[0].hist(training_data, bins=bins, alpha=0.7, color='red', edgecolor='black')
    axs[0].set_title('Training Data', fontsize=14)
    axs[0].set_xlabel('Damping Values', fontsize=12)
    axs[0].set_ylabel('Amount of OPs (log scale)', fontsize=12)
    #axs[0].set_yscale('log')  # Set y-axis to logarithmic scale
    max_y = max(max_y, axs[0].get_ylim()[1])

    axs[1].hist(sec_training_data, bins=bins, alpha=0.7, color='purple', edgecolor='black')
    axs[1].set_title('Secure training Data', fontsize=14)
    axs[1].set_xlabel('Damping Values', fontsize=12)
    #axs[1].set_yscale('log')  # Set y-axis to logarithmic scale
    max_y = max(max_y, axs[1].get_ylim()[1])

    # Plot the histogram for 'training_data'
    axs[2].hist(insec_training_data, bins=bins, alpha=0.7, color='green', edgecolor='black')
    axs[2].set_title('Insecure training Data', fontsize=14)
    axs[2].set_xlabel('Damping Values', fontsize=12)
    #axs[2].set_yscale('log')  # Set y-axis to logarithmic scale
    max_y = max(max_y, axs[2].get_ylim()[1])

    # Plot the histogram for 'damping_data'
    axs[3].hist(damping_data, bins=bins, alpha=0.7, color='blue', edgecolor='black')
    axs[3].set_title('Testing Data', fontsize=14)
    axs[3].set_xlabel('Damping Values', fontsize=12)
    #axs[3].set_yscale('log')  # Set y-axis to logarithmic scale
    max_y = max(max_y, axs[3].get_ylim()[1])

    # Set the same y-limits for both plots
    for ax in axs:
        ax.set_ylim(1e-1, max_y)  # Adjust lower limit as needed

        # Set x-tick labels to show bin range
        bin_labels = [f"{bins[i]:.3f}-{bins[i + 1]:.3f}" for i in range(len(bins) - 1)]
        ax.set_xticks(bins[:-1] + bar_width / 2)  # Center ticks between bin edges
        ax.set_xticklabels(bin_labels, rotation=45, fontsize=8)
        
        # Style adjustments
        ax.grid(True, linestyle='--', alpha=0.6)  # Add grid lines for readability
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Add a title to the entire figure
    fig.suptitle(title, fontsize=16)

    # Adjust spacing to prevent label overlap
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()


# specify test set: 'all', 'HIC'
case_size = "39"
method_type = "SD10"

test_index = "all"
equal_split = 'false'

# Specify the directory and file name
directory = f"C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/output/case{case_size}/datasets/"
flow_name_method = f'{case_size}bus_{method_type}_method_flows.csv' # DT method data
flow_name_lhc = f'{case_size}bus_lhc_flows.csv'
flow_name_imp = f'{case_size}bus_imp_flows.csv' 

file_name_method = f'{case_size}bus_{method_type}_method_ops.csv' # DT method data
file_name_lhc = f'{case_size}bus_lhc_ops.csv'
file_name_imp = f'{case_size}bus_imp_ops.csv' 

file_name_result = 'test.csv'
seed = 42

# define stability boundary
stability_boundary = 0.03
lower_bound = 0.0275
upper_bound = 0.0325

# number of samples
total_samples = 10000
secure_percentage = 0.5
total_secure = int(secure_percentage*total_samples)
total_insecure = int((1-secure_percentage)*total_samples)

# Create the full file path
flow_path_method = os.path.join(directory, flow_name_method)
flow_path_lhc = os.path.join(directory, flow_name_lhc)
flow_path_imp = os.path.join(directory, flow_name_imp)

file_path_method = os.path.join(directory, file_name_method)
file_path_lhc = os.path.join(directory, file_name_lhc)
file_path_imp = os.path.join(directory, file_name_imp)

file_path_result = os.path.join(directory, file_name_result)

# Load the OP dataset

###############
op_data_method = pd.read_csv(file_path_method, sep = ';')
op_data_method_original = op_data_method
#### only HIC
if test_index == 'HIC':
    op_data_method = op_data_method[(op_data_method['damping'] >= lower_bound) & (op_data_method['damping'] <= upper_bound)]

if equal_split == 'only_method' or equal_split == 'true':
    ####################### filter to get 50/50 secure insecure split
    # Filter the DataFrame based on the two conditions
    data_above_lower_bound = op_data_method[(op_data_method['damping'] >= stability_boundary) & (op_data_method['N0'] == 1)]
    data_below_lower_bound = op_data_method[(op_data_method['damping'] < stability_boundary) | (op_data_method['N0'] != 1)]

    # Sample 10,000 rows from each subset (indices are preserved)
    sample_above = data_above_lower_bound.sample(n=total_secure, random_state=seed).reset_index()
    sample_below = data_below_lower_bound.sample(n=total_insecure, random_state=seed).reset_index()

    # Concatenate the two sampled DataFrames
    op_data_method = pd.concat([sample_above, sample_below])
    damping_data_method = op_data_method['damping']
    ##############################################################
else:
    op_data_method = op_data_method.sample(n = total_samples, random_state=seed).reset_index(drop=False)
    damping_data_method = op_data_method['damping']


###############
op_data_method_hic = op_data_method_original[(op_data_method_original['damping'] >= 0.029) & (op_data_method_original['damping'] <= 0.031)]
op_data_method_hic = op_data_method_hic.sample(n = int(0.25*total_samples), random_state=seed).reset_index(drop=False)
damping_data_method_hic = op_data_method_hic['damping']

###############
op_data_lhc = pd.read_csv(file_path_lhc, sep = ';')
op_data_lhc_original = op_data_lhc
#### only HIC
#op_data_lhc = op_data_lhc[(op_data_lhc['damping'] >= 0.0275) & (op_data_lhc['damping'] <= 0.0325)]

if equal_split == 'true':
    ####################### filter to get 50/50 secure insecure split
    # Filter the DataFrame based on the two conditions
    data_above_lower_bound = op_data_lhc[(op_data_lhc['damping'] >= stability_boundary) & (op_data_lhc['N0'] == 1)]
    data_below_lower_bound = op_data_lhc[(op_data_lhc['damping'] < stability_boundary) | (op_data_lhc['N0'] != 1)]

    # Sample 10,000 rows from each subset (indices are preserved)
    sample_above = data_above_lower_bound.sample(n=total_secure, random_state=seed).reset_index()
    sample_below = data_below_lower_bound.sample(n=total_insecure, random_state=seed).reset_index()

    # Concatenate the two sampled DataFrames
    op_data_lhc = pd.concat([sample_above, sample_below])
    damping_data_lhc = op_data_lhc['damping']
    ##############################################################
elif equal_split == 'only_method':
#     ####################### filter to get 50/50 secure insecure split
#     # Filter the DataFrame based on the two conditions
#     data_above_lower_bound = op_data_lhc[(op_data_lhc['damping'] >= stability_boundary) & (op_data_lhc['N0'] == 1)]
#     data_below_lower_bound = op_data_lhc[(op_data_lhc['damping'] < stability_boundary) | (op_data_lhc['N0'] != 1)]

#     # Sample 10,000 rows from each subset (indices are preserved)
#     sample_above = data_above_lower_bound.sample(n=int(0.5*total_samples), random_state=seed).reset_index()
#     sample_below = data_below_lower_bound.sample(n=int(0.5*total_samples), random_state=seed).reset_index()

#     # Concatenate the two sampled DataFrames
#     op_data_lhc = pd.concat([sample_above, sample_below])
#     damping_data_lhc = op_data_lhc['damping']
#     ##############################################################
    op_data_lhc = op_data_lhc.sample(n = total_samples, random_state=seed).reset_index(drop=False)
    damping_data_lhc = op_data_lhc['damping']
else:
    op_data_lhc = op_data_lhc.sample(n = total_samples, random_state=seed).reset_index(drop=False)
    damping_data_lhc = op_data_lhc['damping']


###############
op_data_lhc_far = op_data_lhc_original[(op_data_lhc_original['damping'] < 0.02) | (op_data_lhc_original['damping'] > 0.04)]
op_data_lhc_far = op_data_lhc_far.sample(n = int(0.25*total_samples), random_state=seed).reset_index(drop=False)
damping_data_lhc_far = op_data_lhc_far['damping']

##############
op_data_imp = pd.read_csv(file_path_imp, sep = ';')
op_data_imp_original = op_data_imp
#### only HIC
#op_data_imp = op_data_imp[(op_data_imp['damping'] >= 0.0275) & (op_data_imp['damping'] <= 0.0325)]

if equal_split == 'true':
    ####################### filter to get 50/50 secure insecure split
    # Filter the DataFrame based on the two conditions
    data_above_lower_bound = op_data_imp[(op_data_imp['damping'] >= stability_boundary) & (op_data_imp['N0'] == 1)]
    data_below_lower_bound = op_data_imp[(op_data_imp['damping'] < stability_boundary) | (op_data_imp['N0'] != 1)]

    # Sample 10,000 rows from each subset (indices are preserved)
    sample_above = data_above_lower_bound.sample(n=total_secure, random_state=seed).reset_index()
    sample_below = data_below_lower_bound.sample(n=total_insecure, random_state=seed).reset_index()

    # Concatenate the two sampled DataFrames
    op_data_imp = pd.concat([sample_above, sample_below])
    damping_data_imp = op_data_imp['damping']
    ##############################################################
elif equal_split == 'only_method':
#     ####################### filter to get 50/50 secure insecure split
#     # Filter the DataFrame based on the two conditions
#     data_above_lower_bound = op_data_imp[(op_data_imp['damping'] >= stability_boundary) & (op_data_imp['N0'] == 1)]
#     data_below_lower_bound = op_data_imp[(op_data_imp['damping'] < stability_boundary) | (op_data_imp['N0'] != 1)]

#     # Sample 10,000 rows from each subset (indices are preserved)
#     sample_above = data_above_lower_bound.sample(n=int(0.5*total_samples), random_state=seed).reset_index()
#     sample_below = data_below_lower_bound.sample(n=int(0.5*total_samples), random_state=seed).reset_index()

#     # Concatenate the two sampled DataFrames
#     op_data_imp = pd.concat([sample_above, sample_below])
#     damping_data_imp = op_data_imp['damping']
#     ##############################################################
    op_data_imp = op_data_imp.sample(n = total_samples, random_state=seed).reset_index(drop=False)
    damping_data_imp = op_data_imp['damping']
else:
    op_data_imp = op_data_imp.sample(n = total_samples, random_state=seed).reset_index(drop=False)
    damping_data_imp = op_data_imp['damping']


# print indices
print(op_data_method['index'])
print(op_data_imp['index'])
print(op_data_lhc['index'])


# Load the flow dataset
############# method
data_method = pd.read_csv(flow_path_method, sep = ';')
data_method = data_method.iloc[op_data_method['index']]
X_method = data_method.drop(columns = ['feasible', 'stable'], axis=1) # drop columns you won't use like ['Feas', 'N1 ....]
y_method = ((data_method['stable'] == 1) & (data_method['feasible'] == 1)).astype(int)
X_train_method, X_test_method, y_train_method, y_test_method, damping_data_method_train, damping_data_method_test, op_data_method_train, op_data_method_test = train_test_split(X_method, y_method, damping_data_method, op_data_method, test_size=0.25, random_state=seed)
method_secure = (y_method == 1).sum()
method_insecure = (y_method == 0).sum()

method_stable = (op_data_method['damping'] >= 0.03).sum()
method_instable = (op_data_method['damping'] <= 0.03).sum()

method_hic = ((op_data_method['damping'] >= 0.0275) & (op_data_method['damping'] <= 0.0325)).sum()

method_feasible = (op_data_method['N0'] == 1).sum()
method_infeasible = (op_data_method['N0'] == 0).sum()
print("method secure: ", method_secure)
print("method insecure: ", method_insecure)

################## method hic
data_method = pd.read_csv(flow_path_method, sep = ';')
data_method_hic = data_method.iloc[op_data_method_hic['index']]
X_method_hic = data_method_hic.drop(columns = ['feasible', 'stable'], axis=1) # drop columns you won't use like ['Feas', 'N1 ....]
y_method_hic = ((data_method_hic['stable'] == 1) & (data_method_hic['feasible'] == 1)).astype(int)
X_train_method_hic, X_test_method_hic, y_train_method_hic, y_test_method_hic, damping_data_method_train_hic, damping_data_method_test_hic, op_data_method_train_hic, op_data_method_test_hic = train_test_split(X_method_hic, y_method_hic, damping_data_method_hic, op_data_method_hic, test_size=0.25, random_state=seed)
print((y_method_hic == 0).sum())
print((y_method_hic == 1).sum())

################### lhc
data_lhc = pd.read_csv(flow_path_lhc, sep = ';')
data_lhc = data_lhc.iloc[op_data_lhc['index']]
X_lhc = data_lhc.drop(columns = ['feasible', 'stable'], axis=1) # drop columns you won't use like ['Feas', 'N1 ....]
y_lhc = ((data_lhc['stable'] == 1) & (data_lhc['feasible'] == 1)).astype(int)
X_train_lhc, X_test_lhc, y_train_lhc, y_test_lhc, damping_data_lhc_train, damping_data_lhc_test, op_data_lhc_train, op_data_lhc_test = train_test_split(X_lhc, y_lhc, damping_data_lhc, op_data_lhc, test_size=0.25, random_state=seed)
lhc_secure = (y_lhc == 1).sum()
lhc_insecure = (y_lhc == 0).sum()

lhc_stable = (op_data_lhc['damping'] >= 0.03).sum()
lhc_instable = (op_data_lhc['damping'] <= 0.03).sum()

lhc_hic = ((op_data_lhc['damping'] >= 0.0275) & (op_data_lhc['damping'] <= 0.0325)).sum()

lhc_feasible = (op_data_lhc['N0'] == 1).sum()
lhc_infeasible = (op_data_lhc['N0'] == 0).sum()
print("lhc secure: ", lhc_secure)
print("lhc insecure: ", lhc_insecure)

##################### lhc far
data_lhc = pd.read_csv(flow_path_lhc, sep = ';')
data_lhc_far = data_lhc.iloc[op_data_lhc_far['index']]
X_lhc_far = data_lhc_far.drop(columns = ['feasible', 'stable'], axis=1) # drop columns you won't use like ['Feas', 'N1 ....]
y_lhc_far = ((data_lhc_far['stable'] == 1) & (data_lhc_far['feasible'] == 1)).astype(int)
X_train_lhc_far, X_test_lhc_far, y_train_lhc_far, y_test_lhc_far, damping_data_lhc_far_train, damping_data_lhc_far_test, op_data_lhc_far_train, op_data_lhc_far_test = train_test_split(X_lhc_far, y_lhc_far, damping_data_lhc_far, op_data_lhc_far, test_size=0.25, random_state=seed)
lhc_secure_far = (y_lhc_far == 1).sum()
lhc_insecure_far = (y_lhc_far == 0).sum()


################# imp
data_imp = pd.read_csv(flow_path_imp, sep = ';')
data_imp = data_imp.iloc[op_data_imp['index']]
X_imp = data_imp.drop(columns = ['feasible', 'stable'], axis=1) # drop columns you won't use like ['Feas', 'N1 ....]
y_imp = ((data_imp['stable'] == 1) & (data_imp['feasible'] == 1)).astype(int)
X_train_imp, X_test_imp, y_train_imp, y_test_imp,  damping_data_imp_train, damping_data_imp_test, op_data_imp_train, op_data_imp_test = train_test_split(X_imp, y_imp, damping_data_imp, op_data_imp, test_size=0.25, random_state=seed)
imp_secure = (y_imp == 1).sum()
imp_insecure = (y_imp == 0).sum()

imp_stable = (op_data_imp['damping'] >= 0.03).sum()
imp_instable = (op_data_imp['damping'] <= 0.03).sum()

imp_hic = ((op_data_imp['damping'] >= 0.0275) & (op_data_imp['damping'] <= 0.0325)).sum()

imp_feasible = (op_data_imp['N0'] == 1).sum()
imp_infeasible = (op_data_imp['N0'] == 0).sum()
print("imp secure: ", imp_secure)
print("imp insecure: ", imp_insecure)

# Clear the CSV file before starting the run (optional)
if os.path.exists(file_path_result):
    with open(file_path_result, 'w') as f:
        f.write('')  # Write an empty string to clear the file

# iteratively reduce training set size
reduce_train = 1
initial_fraction = 1.0  # Start with the full dataset
reduction_step = 0.0#0.05    # Reduce by this fraction in each iteration

for i in range(reduce_train):
    # Calculate the new fraction for the current iteration
    current_fraction = initial_fraction - (i * reduction_step)
    
    # Ensure the fraction doesn't go below a certain limit
    current_fraction = max(current_fraction, 0.1)  # Keep at least 10% of the data

    # Randomly sample a subset of the training sets
    X_train_method_subset = X_train_method.sample(frac=current_fraction, random_state=seed)
    y_train_method_subset = y_train_method.loc[X_train_method_subset.index]  # Ensure labels match

    # X_train_method_hic_subset = X_train_method_hic.sample(frac=current_fraction, random_state=42)
    # y_train_method_hic_subset = y_train_method_hic.loc[X_train_method_hic_subset.index]  # Ensure labels match
    
    X_train_lhc_subset = X_train_lhc.sample(frac=current_fraction, random_state=seed)
    y_train_lhc_subset = y_train_lhc.loc[X_train_lhc_subset.index]  # Ensure labels match
    
    X_train_imp_subset = X_train_imp.sample(frac=current_fraction, random_state=seed)
    y_train_imp_subset = y_train_imp.loc[X_train_imp_subset.index]  # Ensure labels match

    ccp_alpha = [0.001, 0.005, 0.01, 0.02, 0.03]

    # Initialize the DecisionTreeClassifier with gini impurity and max depth of 5. clf = classifier
    clf_method = DecisionTreeClassifier(criterion='gini', max_depth=5, ccp_alpha=0.01, random_state=seed) # splitter='best' or 'random' , ccp_alpha=0.01 for tree pruning
    clf_lhc = DecisionTreeClassifier(criterion='gini', max_depth=5, ccp_alpha=0.01, random_state=seed) # criterion = 'gini' or 'entropy'
    clf_imp = DecisionTreeClassifier(criterion='gini', max_depth=5, ccp_alpha=0.01, random_state=seed)

    # Perform 10-fold cross-validation and compute accuracy for each fold
    cv_scores_method = cross_val_score(clf_method, X_train_method_subset, y_train_method_subset, cv=10, scoring='accuracy')
    cv_scores_lhc = cross_val_score(clf_lhc, X_train_lhc_subset, y_train_lhc_subset, cv=10, scoring='accuracy')
    cv_scores_imp = cross_val_score(clf_imp, X_train_imp_subset, y_train_imp_subset, cv=10, scoring='accuracy')

    print(f"10-Fold CV Accuracy for clf_method: {cv_scores_method.mean():.2f} ± {cv_scores_method.std():.2f}")
    print(f"10-Fold CV Accuracy for clf_lhc: {cv_scores_lhc.mean():.2f} ± {cv_scores_lhc.std():.2f}")
    print(f"10-Fold CV Accuracy for clf_imp: {cv_scores_imp.mean():.2f} ± {cv_scores_imp.std():.2f}")

    # Train the classifier on the reduced training sets
    clf_method.fit(X_train_method_subset, y_train_method_subset)
    clf_lhc.fit(X_train_lhc_subset, y_train_lhc_subset)
    clf_imp.fit(X_train_imp_subset, y_train_imp_subset)

    ########### Make predictions on the training and test sets and check for overfitting
    y_train_pred_method = clf_method.predict(X_train_method_subset)
    y_test_pred_method = clf_method.predict(X_test_method)
    y_pred_method = y_test_pred_method

    y_train_pred_lhc = clf_lhc.predict(X_train_lhc_subset)
    y_test_pred_lhc = clf_lhc.predict(X_test_lhc)
    y_pred_lhc = y_test_pred_lhc

    y_train_pred_imp = clf_imp.predict(X_train_imp_subset)
    y_test_pred_imp = clf_imp.predict(X_test_imp)
    y_pred_imp = y_test_pred_imp

    # Calculate train and test accuracies for each classifier
    train_accuracy_method = accuracy_score(y_train_method_subset, y_train_pred_method)
    test_accuracy_method = accuracy_score(y_test_method, y_test_pred_method)

    train_accuracy_lhc = accuracy_score(y_train_lhc_subset, y_train_pred_lhc)
    test_accuracy_lhc = accuracy_score(y_test_lhc, y_test_pred_lhc)

    train_accuracy_imp = accuracy_score(y_train_imp_subset, y_train_pred_imp)
    test_accuracy_imp = accuracy_score(y_test_imp, y_test_pred_imp)

    # Print accuracy results to compare train and test accuracy
    print("\nTraining and Test Accuracies:")
    print(f"Method Classifier - Train Accuracy: {train_accuracy_method:.2f}, Test Accuracy: {test_accuracy_method:.2f}")
    print(f"LHC Classifier - Train Accuracy: {train_accuracy_lhc:.2f}, Test Accuracy: {test_accuracy_lhc:.2f}")
    print(f"IMP Classifier - Train Accuracy: {train_accuracy_imp:.2f}, Test Accuracy: {test_accuracy_imp:.2f}")

    # Interpretation for overfitting/underfitting
    def interpret_overfitting(train_acc, test_acc, name):
        if train_acc > test_acc:
            print(f"{name} might be overfitting (Train Accuracy > Test Accuracy).")
        elif test_acc > train_acc:
            print(f"{name} might be underfitting (Test Accuracy > Train Accuracy).")
        else:
            print(f"{name} is performing well with balanced Train and Test accuracies.")

    # Interpret results for each classifier
    print("\nOverfitting Interpretation:")
    interpret_overfitting(train_accuracy_method, test_accuracy_method, "Method Classifier")
    interpret_overfitting(train_accuracy_lhc, test_accuracy_lhc, "LHC Classifier")
    interpret_overfitting(train_accuracy_imp, test_accuracy_imp, "IMP Classifier")

    ################# cross tests
    y_pred_method_lhc = clf_method.predict(X_test_lhc)
    y_pred_method_lhc_far = clf_method.predict(X_test_lhc_far)
    y_pred_method_imp = clf_method.predict(X_test_imp)
    y_pred_method_hic = clf_method.predict(X_test_method_hic)
    # y_prob_method = clf_method.predict_proba(X_test_method) # due to early stopping, there is a probability of being wrong. Not all leafs are pure. 

    y_pred_lhc_method = clf_lhc.predict(X_test_method)
    y_pred_lhc_method_hic = clf_lhc.predict(X_test_method_hic)
    y_pred_lhc_imp = clf_lhc.predict(X_test_imp)
    y_pred_lhc_far = clf_lhc.predict(X_test_lhc_far)
    # y_prob_lhc = clf_lhc.predict_proba(X_test_lhc)

    y_pred_imp_method = clf_imp.predict(X_test_method)
    y_pred_imp_method_hic = clf_imp.predict(X_test_method_hic)
    y_pred_imp_lhc = clf_imp.predict(X_test_lhc)
    y_pred_imp_lhc_far = clf_imp.predict(X_test_lhc_far)
    # y_prob_imp = clf_imp.predict_proba(X_test_imp)


    # Evaluate the classifier accuracy
    accuracy_method = accuracy_score(y_test_method, y_pred_method)
    accuracy_method_hic = accuracy_score(y_test_method_hic, y_pred_method_hic)
    accuracy_method_lhc = accuracy_score(y_test_lhc, y_pred_method_lhc)
    accuracy_method_lhc_far = accuracy_score(y_test_lhc_far, y_pred_method_lhc_far)
    accuracy_method_imp = accuracy_score(y_test_imp, y_pred_method_imp)
    print("Test accuracy train  method, test method:", accuracy_method)
    print("Test accuracy train  method, test lhc:", accuracy_method_lhc)
    print("Test accuracy train  method, test imp:", accuracy_method_imp)

    accuracy_lhc = accuracy_score(y_test_lhc, y_pred_lhc)
    accuracy_lhc_far = accuracy_score(y_test_lhc_far, y_pred_lhc_far)
    accuracy_lhc_method = accuracy_score(y_test_method, y_pred_lhc_method)
    accuracy_lhc_method_hic = accuracy_score(y_test_method_hic, y_pred_lhc_method_hic)
    accuracy_lhc_imp = accuracy_score(y_test_imp, y_pred_lhc_imp)
    print("Test accuracy train lhc, test lhc:", accuracy_lhc)
    print("Test accuracy train lhc, test method:", accuracy_lhc_method)
    print("Test accuracy train lhc, test imp:", accuracy_lhc_imp)

    accuracy_imp = accuracy_score(y_test_imp, y_pred_imp)
    accuracy_imp_method = accuracy_score(y_test_method, y_pred_imp_method)
    accuracy_imp_method_hic = accuracy_score(y_test_method_hic, y_pred_imp_method_hic)
    accuracy_imp_lhc = accuracy_score(y_test_lhc, y_pred_imp_lhc)
    accuracy_imp_lhc_far = accuracy_score(y_test_lhc_far, y_pred_imp_lhc_far)
    print("Test accuracy train imp, test imp:", accuracy_imp)
    print("Test accuracy train imp, test method:", accuracy_imp_method)
    print("Test accuracy train imp, test lhc:", accuracy_imp_lhc)

    # Evaluate the classifier precision
    precision_method = precision_score(y_test_method, y_pred_method)
    precision_method_hic = precision_score(y_test_method_hic, y_pred_method_hic)
    precision_method_lhc = precision_score(y_test_lhc, y_pred_method_lhc)
    precision_method_imp = precision_score(y_test_imp, y_pred_method_imp)
    print("Test precision train  method, test method:", precision_method)
    print("Test precision train  method, test lhc:", precision_method_lhc)
    print("Test precision train  method, test imp:", precision_method_imp)

    precision_lhc = precision_score(y_test_lhc, y_pred_lhc)
    precision_lhc_method = precision_score(y_test_method, y_pred_lhc_method)
    precision_lhc_method_hic = precision_score(y_test_method_hic, y_pred_lhc_method_hic)
    precision_lhc_imp = precision_score(y_test_imp, y_pred_lhc_imp)
    print("Test precision train lhc, test lhc:", precision_lhc)
    print("Test precision train lhc, test method:", precision_lhc_method)
    print("Test precision train lhc, test imp:", precision_lhc_imp)

    precision_imp = precision_score(y_test_imp, y_pred_imp)
    precision_imp_method = precision_score(y_test_method, y_pred_imp_method)
    precision_imp_method_hic = precision_score(y_test_method_hic, y_pred_imp_method_hic)
    precision_imp_lhc = precision_score(y_test_lhc, y_pred_imp_lhc)
    print("Test precision train imp, test imp:", precision_imp)
    print("Test precision train imp, test method:", precision_imp_method)
    print("Test precision train imp, test lhc:", precision_imp_lhc)


    # Evaluate the classifier recall
    recall_method = recall_score(y_test_method, y_pred_method)
    recall_method_hic = recall_score(y_test_method_hic, y_pred_method_hic)
    recall_method_lhc = recall_score(y_test_lhc, y_pred_method_lhc)
    recall_method_imp = recall_score(y_test_imp, y_pred_method_imp)
    print("Test recall train  method, test method:", recall_method)
    print("Test recall train  method, test lhc:", recall_method_lhc)
    print("Test recall train  method, test imp:", recall_method_imp)

    recall_lhc = recall_score(y_test_lhc, y_pred_lhc)
    recall_lhc_method = recall_score(y_test_method, y_pred_lhc_method)
    recall_lhc_method_hic = recall_score(y_test_method_hic, y_pred_lhc_method_hic)
    recall_lhc_imp = recall_score(y_test_imp, y_pred_lhc_imp)
    print("Test recall train lhc, test lhc:", recall_lhc)
    print("Test recall train lhc, test method:", recall_lhc_method)
    print("Test recall train lhc, test imp:", recall_lhc_imp)

    recall_imp = recall_score(y_test_imp, y_pred_imp)
    recall_imp_method = recall_score(y_test_method, y_pred_imp_method)
    recall_imp_method_hic = recall_score(y_test_method_hic, y_pred_imp_method_hic)
    recall_imp_lhc = recall_score(y_test_lhc, y_pred_imp_lhc)
    print("Test recall train imp, test imp:", recall_imp)
    print("Test recall train imp, test method:", recall_imp_method)
    print("Test recall train imp, test lhc:", recall_imp_lhc)


    # Evaluate the classifier F1 score
    f1_method = f1_score(y_test_method, y_pred_method)
    f1_method_hic = f1_score(y_test_method_hic, y_pred_method_hic)
    f1_method_lhc = f1_score(y_test_lhc, y_pred_method_lhc)
    f1_method_lhc_far = f1_score(y_test_lhc_far, y_pred_method_lhc_far)
    f1_method_imp = f1_score(y_test_imp, y_pred_method_imp)
    print("Test f1_score train  method, test method:", f1_method)
    print("Test f1_score train  method, test lhc:", f1_method_lhc)
    print("Test f1_score train  method, test imp:", f1_method_imp)

    f1_lhc = f1_score(y_test_lhc, y_pred_lhc)
    f1_lhc_far = f1_score(y_test_lhc_far, y_pred_lhc_far)
    f1_lhc_method = f1_score(y_test_method, y_pred_lhc_method)
    f1_lhc_method_hic = f1_score(y_test_method_hic, y_pred_lhc_method_hic)
    f1_lhc_imp = f1_score(y_test_imp, y_pred_lhc_imp)
    print("Test f1_score train lhc, test lhc:", f1_lhc)
    print("Test f1_score train lhc, test method:", f1_lhc_method)
    print("Test f1_score train lhc, test imp:", f1_lhc_imp)

    f1_imp = f1_score(y_test_imp, y_pred_imp)
    f1_imp_method = f1_score(y_test_method, y_pred_imp_method)
    f1_imp_method_hic = f1_score(y_test_method_hic, y_pred_imp_method_hic)
    f1_imp_lhc = f1_score(y_test_lhc, y_pred_imp_lhc)
    f1_imp_lhc_far = f1_score(y_test_lhc_far, y_pred_imp_lhc_far)
    print("Test f1_score train imp, test imp:", f1_imp)
    print("Test f1_score train imp, test method:", f1_imp_method)
    print("Test f1_score train imp, test lhc:", f1_imp_lhc)

    # Evaluate the classifier false positive rate. det_curve = fpr, fnr, threshold
    fpr_method, fnr_method, FP_list_method, FN_list_method = fpr_score(y_test_method, y_pred_method)
    fpr_method_hic, fnr_method_hic, FP_list_method_hic, FN_list_method_hic = fpr_score(y_test_method_hic, y_pred_method_hic)
    fpr_method_lhc, fnr_method_lhc, FP_list_method_lhc, FN_list_method_lhc = fpr_score(y_test_lhc, y_pred_method_lhc)
    fpr_method_lhc_far, fnr_method_lhc_far, FP_list_method_lhc_far, FN_list_method_lhc_far = fpr_score(y_test_lhc_far, y_pred_method_lhc_far)
    fpr_method_imp, fnr_method_imp, FP_list_method_imp, FN_list_method_imp = fpr_score(y_test_imp, y_pred_method_imp)
    print("Test fpr train  method, test method:", fpr_method)
    print("Test fpr train  method, test lhc:", fpr_method_lhc)
    print("Test fpr train  method, test imp:", fpr_method_imp)

    fpr_lhc, fnr_lhc, FP_list_lhc, FN_list_lhc = fpr_score(y_test_lhc, y_pred_lhc)
    fpr_lhc_far, fnr_lhc_far, FP_list_lhc_far, FN_list_lhc_far = fpr_score(y_test_lhc_far, y_pred_lhc_far)
    fpr_lhc_method, fnr_lhc_method, FP_list_lhc_method, FN_list_lhc_method = fpr_score(y_test_method, y_pred_lhc_method)
    fpr_lhc_method_hic, fnr_lhc_method_hic, FP_list_lhc_method_hic, FN_list_lhc_method_hic = fpr_score(y_test_method_hic, y_pred_lhc_method_hic)
    fpr_lhc_imp, fnr_lhc_imp, FP_list_lhc_imp, FN_list_lhc_imp = fpr_score(y_test_imp, y_pred_lhc_imp)
    print("Test fpr train lhc, test lhc:", fpr_lhc)
    print("Test fpr train lhc, test method:", fpr_lhc_method)
    print("Test fpr train lhc, test imp:", fpr_lhc_imp)

    fpr_imp, fnr_imp, FP_list_imp, FN_list_imp = fpr_score(y_test_imp, y_pred_imp)
    fpr_imp_method, fnr_imp_method, FP_list_imp_method, FN_list_imp_method = fpr_score(y_test_method, y_pred_imp_method)
    fpr_imp_method_hic, fnr_imp_method_hic, FP_list_imp_method_hic, FN_list_imp_method_hic = fpr_score(y_test_method_hic, y_pred_imp_method_hic)
    fpr_imp_lhc, fnr_imp_lhc, FP_list_imp_lhc, FN_list_imp_lhc = fpr_score(y_test_lhc, y_pred_imp_lhc)
    fpr_imp_lhc_far, fnr_imp_lhc_far, FP_list_imp_lhc_far, FN_list_imp_lhc_far = fpr_score(y_test_lhc_far, y_pred_imp_lhc_far)
    print("Test fpr train imp, test imp:", fpr_imp)
    print("Test fpr train imp, test method:", fpr_imp_method)
    print("Test fpr train imp, test lhc:", fpr_imp_lhc)

    # Determine feature importance
    feature_names_method = X_method.columns
    feature_importance_method = clf_method.feature_importances_ # importance of each feature w.r.t. prediction

    feature_names_lhc = X_lhc.columns
    feature_importance_lhc = clf_lhc.feature_importances_

    feature_names_imp = X_imp.columns
    feature_importance_imp = clf_imp.feature_importances_

    # Create a DataFrame to store the feature importance values
    importance_df_method = pd.DataFrame({
        'Feature': feature_names_method,
        'Importance': feature_importance_method
    })

    importance_df_lhc = pd.DataFrame({
        'Feature': feature_names_lhc,
        'Importance': feature_importance_lhc
    })

    importance_df_imp = pd.DataFrame({
        'Feature': feature_names_imp,
        'Importance': feature_importance_imp
    })

    # Sort by importance (descending order)
    importance_df_method = importance_df_method.sort_values(by='Importance', ascending=False)
    importance_df_lhc = importance_df_lhc.sort_values(by='Importance', ascending=False)
    importance_df_imp = importance_df_imp.sort_values(by='Importance', ascending=False)

    # Select features with importance greater than 0
    features_method = list(importance_df_method[importance_df_method['Importance'] > 0]['Feature'])
    features_lhc = list(importance_df_lhc[importance_df_lhc['Importance'] > 0]['Feature'])
    features_imp = list(importance_df_imp[importance_df_imp['Importance'] > 0]['Feature'])


    data_result = {
        'training data': ['method', 'lhc', 'importance'],
        'train size': [len(X_train_method_subset), len(X_train_lhc_subset), len(X_train_imp_subset)],
        'test size': [len(y_test_method), len(y_test_lhc), len(y_test_imp)],
        'space0': ['', '', ''],
        'feasible': [method_feasible, lhc_feasible, imp_feasible],
        'stable': [method_stable, lhc_stable, imp_stable],
        'secure': [method_secure, lhc_secure, imp_secure],
        'HIC': [method_hic, lhc_hic, imp_hic],
        'spacey': ['', '', ''],
        'method test accuracy': [accuracy_method, accuracy_lhc_method, accuracy_imp_method],
        'method hic test accuracy': [accuracy_method_hic, accuracy_lhc_method_hic, accuracy_imp_method_hic],
        'lhc test accuracy': [accuracy_method_lhc, accuracy_lhc, accuracy_imp_lhc],
        'lhc far test accuracy': [accuracy_method_lhc_far, accuracy_lhc_far, accuracy_imp_lhc_far],
        'imp test accuracy': [accuracy_method_imp, accuracy_lhc_imp, accuracy_imp],
        'space': ['', '', ''],
        'method test recall': [recall_method, recall_lhc_method, recall_imp_method],
        'method hic test recall': [recall_method_hic, recall_lhc_method_hic, recall_imp_method_hic],
        'lhc test recall': [recall_method_lhc, recall_lhc, recall_imp_lhc],
        'imp test recall': [recall_method_imp, recall_lhc_imp, recall_imp],
        'space1': ['', '', ''],
        'method test precision': [precision_method, precision_lhc_method, precision_imp_method],
        'method hic test precision': [precision_method_hic, precision_lhc_method_hic, precision_imp_method_hic],
        'lhc test precision': [precision_method_lhc, precision_lhc, precision_imp_lhc],
        'imp test precision': [precision_method_imp, precision_lhc_imp, precision_imp],
        'space2': ['', '', ''],
        'method test f1': [f1_method, f1_lhc_method, f1_imp_method],
        'method hic test f1': [f1_method_hic, f1_lhc_method_hic, f1_imp_method_hic],
        'lhc test f1': [f1_method_lhc, f1_lhc, f1_imp_lhc],
        'lhc far test f1': [f1_method_lhc_far, f1_lhc_far, f1_imp_lhc_far],
        'imp test f1': [f1_method_imp, f1_lhc_imp, f1_imp],
        'space3': ['', '', ''],
        'method test fpr': [fpr_method, fpr_lhc_method, fpr_imp_method],
        'method hic test fpr': [fpr_method_hic, fpr_lhc_method_hic, fpr_imp_method_hic],
        'lhc test fpr': [fpr_method_lhc, fpr_lhc, fpr_imp_lhc],
        'lhc far test fpr': [fpr_method_lhc_far, fpr_lhc_far, fpr_imp_lhc_far],
        'imp test fpr': [fpr_method_imp, fpr_lhc_imp, fpr_imp],
        'space4': ['', '', ''],
        'method test fnr': [fnr_method, fnr_lhc_method, fnr_imp_method],
        'method hic test fnr': [fnr_method_hic, fnr_lhc_method_hic, fnr_imp_method_hic],
        'lhc test fnr': [fnr_method_lhc, fnr_lhc, fnr_imp_lhc],
        'lhc test fnr': [fnr_method_lhc_far, fnr_lhc_far, fnr_imp_lhc_far],
        'imp test fnr': [fnr_method_imp, fnr_lhc_imp, fnr_imp],

    }
        
    df_result = pd.DataFrame.from_dict(data_result)

    # Append a newline before writing new results
    with open(file_path_result, 'a', newline='') as f:
        if f.tell() > 0:  # Check if the file already has content
            f.write('\n')  # Add a blank line before appending the new DataFrame
        df_result.to_csv(f, sep=';', header=f.tell() == 0, index=False)





y_method_sec = (y_train_method == 1).reset_index(drop=True)#.astype(int)
y_lhc_sec = (y_train_lhc == 1).reset_index(drop=True)#.astype(int)
y_imp_sec = (y_train_imp == 1).reset_index(drop=True)#.astype(int)

print(y_method_sec.sort_index(ascending=True), damping_data_method_train.sort_index(ascending=True))

damping_data_method_train= damping_data_method_train.reset_index(drop=True)
damping_data_lhc_train= damping_data_lhc_train.reset_index(drop=True)
damping_data_imp_train= damping_data_imp_train.reset_index(drop=True)

damping_secure_data_method = damping_data_method_train[y_method_sec]
damping_secure_data_lhc = damping_data_lhc_train[y_lhc_sec]
damping_secure_data_imp = damping_data_imp_train[y_imp_sec]

y_method_ins = (y_train_method == 0).reset_index(drop=True)
y_lhc_ins = (y_train_lhc == 0).reset_index(drop=True)
y_imp_ins = (y_train_imp == 0).reset_index(drop=True)

damping_insecure_data_method = damping_data_method_train[y_method_ins]
damping_insecure_data_lhc = damping_data_lhc_train[y_lhc_ins]
damping_insecure_data_imp = damping_data_imp_train[y_imp_ins]

damping_all_test_data = pd.concat([damping_data_method_test, damping_data_lhc_test, damping_data_imp_test])


# plot the damping of the missclassified OPs
x_lb = -0.02
bar_width = 0.005
plot_damping_histograms(x_lb, bar_width, damping_all_test_data, damping_data_method_train, damping_secure_data_method, damping_insecure_data_method, "Proposed method samples - damping vs number of OPs")
plot_damping_histograms(x_lb, bar_width, damping_all_test_data, damping_data_lhc_train, damping_secure_data_lhc, damping_insecure_data_lhc, "LHC samples - damping vs number of OPs")
plot_damping_histograms(x_lb, bar_width, damping_all_test_data, damping_data_imp_train, damping_secure_data_imp, damping_insecure_data_imp, "Importance samples - damping vs number of OPs")
# #damping_all_test_data



##############################################################################

fp_method_ops = pd.concat([op_data_method_original.iloc[FP_list_method], op_data_lhc_original.iloc[FP_list_method_lhc], op_data_imp_original.iloc[FP_list_method_imp]])
fp_lhc_ops = pd.concat([op_data_lhc_original.loc[FP_list_lhc], op_data_method_original.loc[FP_list_lhc_method], op_data_imp_original.loc[FP_list_lhc_imp]])
fp_imp_ops = pd.concat([op_data_imp_original.loc[FP_list_imp], op_data_method_original.loc[FP_list_imp_method], op_data_lhc_original.loc[FP_list_imp_lhc]])


fn_method_ops = pd.concat([op_data_method_original.loc[FN_list_method], op_data_lhc_original.loc[FN_list_method_lhc], op_data_imp_original.loc[FN_list_method_imp]])
fn_lhc_ops = pd.concat([op_data_lhc_original.loc[FN_list_lhc], op_data_method_original.loc[FN_list_lhc_method], op_data_imp_original.loc[FN_list_lhc_imp]])
fn_imp_ops = pd.concat([op_data_imp_original.loc[FN_list_imp], op_data_method_original.loc[FN_list_imp_method], op_data_lhc_original.loc[FN_list_imp_lhc]])

miss_method = pd.concat([fp_method_ops, fn_method_ops])
miss_lhc = pd.concat([fp_lhc_ops, fn_lhc_ops])
miss_imp = pd.concat([fp_imp_ops, fn_imp_ops])

fp_method_cf = fp_method_ops[['N0', 'damping']] # miss_method
fp_lhc_cf = fp_lhc_ops[['N0', 'damping']] # miss_lhc
fp_imp_cf = fp_imp_ops[['N0', 'damping']] # miss_imp

fn_method_cf = fn_method_ops[['N0', 'damping']] # miss_method
fn_lhc_cf = fn_lhc_ops[['N0', 'damping']] # miss_lhc
fn_imp_cf = fn_imp_ops[['N0', 'damping']] # miss_imp

total_test = pd.concat([op_data_method_test, op_data_lhc_test, op_data_imp_test])
total_test = total_test[['N0', 'damping']]

# Calculate frequencies for the three conditions in each dataframe
def calculate_fp_frequencies(df):
    event1 = (df['N0'] == 0).sum()
    event2 = (df['damping'] < stability_boundary).sum()
    event3 = ((df['N0'] == 0) & (df['damping'] < stability_boundary)).sum()
    return [event1, event2, event3]

def calculate_fn_frequencies(df):
    event1 = (df['N0'] == 1).sum()
    event2 = (df['damping'] > stability_boundary).sum()
    event3 = ((df['N0'] == 1) & (df['damping'] > stability_boundary)).sum()
    return [event1, event2, event3]

# Get frequencies for each dataframe
method_fp_counts = calculate_fp_frequencies(fp_method_cf)
lhc_fp_counts = calculate_fp_frequencies(fp_lhc_cf)
imp_fp_counts = calculate_fp_frequencies(fp_imp_cf)

# Prepare data for plotting
data_fp = [method_fp_counts, lhc_fp_counts, imp_fp_counts]

# Get frequencies for each dataframe
method_fn_counts = calculate_fn_frequencies(fn_method_cf)
lhc_fn_counts = calculate_fn_frequencies(fp_lhc_cf)
imp_fn_counts = calculate_fn_frequencies(fp_imp_cf)

# Prepare data for plotting
data_fn = [method_fn_counts, lhc_fn_counts, imp_fn_counts]

##################### false positive

# Labels and colors
labels = ['N0 == 0', 'damping < 0.0', 'N0 == 0 & damping < 0.0']
df_labels = ['Method', 'LHC', 'IMP']
colors = ['#f4a261', '#e76f51', '#e63946']

# Bar plot with spacing between groups
fig, ax = plt.subplots(figsize=(12, 7))
width = 0.25  # Width of individual bars
x = np.arange(len(labels))  # x positions for groups of bars

# Plotting bars with edges
for i, (d, color) in enumerate(zip(data_fp, colors)):
    bars = ax.bar(x + i * width, d, width, label=df_labels[i], color=color, edgecolor='black')

    # Adding data labels above the bars
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, yval, int(yval), 
                ha='center', va='bottom', fontsize=10, fontweight='bold')

# Customize axes and styling
ax.set_xticks(x + width)
ax.set_xticklabels(labels)
ax.set_ylabel('Frequency of Occurrence', fontsize=12)
ax.set_title('Causes for False Positive Classification', fontsize=14, fontweight='bold')
ax.legend(title="Data Source", fontsize=10, title_fontsize=12)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(True, axis='y', linestyle='--', alpha=0.7)

# Increase bottom margin for better layout
plt.subplots_adjust(bottom=0.15)

# Display plot
plt.tight_layout()
plt.show()

################### false negative 

# Labels and colors
labels = ['N0 == 0', 'damping < 0.03', 'N0 == 0 & damping < 0.03']
df_labels = ['Method', 'LHC', 'IMP']
colors = ['#f4a261', '#e76f51', '#e63946']

# Bar plot with spacing between groups
fig, ax = plt.subplots(figsize=(12, 7))
width = 0.25  # Width of individual bars
x = np.arange(len(labels))  # x positions for groups of bars

# Plotting bars with edges
for i, (d, color) in enumerate(zip(data_fn, colors)):
    bars = ax.bar(x + i * width, d, width, label=df_labels[i], color=color, edgecolor='black')

    # Adding data labels above the bars
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, yval, int(yval), 
                ha='center', va='bottom', fontsize=10, fontweight='bold')

# Customize axes and styling
ax.set_xticks(x + width)
ax.set_xticklabels(labels)
ax.set_ylabel('Frequency of Occurrence', fontsize=12)
ax.set_title('Causes for False Negative Classification', fontsize=14, fontweight='bold')
ax.legend(title="Data Source", fontsize=10, title_fontsize=12)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(True, axis='y', linestyle='--', alpha=0.7)

# Increase bottom margin for better layout
plt.subplots_adjust(bottom=0.15)

# Display plot
plt.tight_layout()
plt.show()



####################################################################


fp_method_violations = fp_method_ops[['N0', 'damping']] # miss_method
fp_lhc_violations = fp_lhc_ops[['N0', 'damping']] # miss_lhc
fp_imp_violations = fp_imp_ops[['N0', 'damping']] # miss_imp

# Parameters
bins = np.arange(0.0, 0.05, 0.0025)
titles = ['Method', 'LHC', 'Importance']
dataframes = [fp_method_violations, fp_lhc_violations, fp_imp_violations]

# Define professional colors
colors = ["#1b9e77", "#d95f02", "#7570b3"]

# Initialize maximum count tracker
max_count = 0

# Precompute maximum count across all datasets
for df in dataframes:
    for i in range(len(bins) - 1):
        bin_df = df[(df['damping'] >= bins[i]) & (df['damping'] < bins[i + 1])]
        both_conditions = ((bin_df['damping'] < 0.03) & (bin_df['N0'] == 0)).sum()
        stability_only = ((bin_df['damping'] < 0.03) & ~(bin_df['N0'] == 0)).sum()
        feasibility_only = ((bin_df['N0'] == 0) & ~(bin_df['damping'] < 0.03)).sum()
        count = both_conditions + stability_only + feasibility_only
        max_count = max(max_count, count)

# Set consistent y-axis limit with some padding
y_limit = max_count + 10  # Add padding for aesthetics

# Create a figure with subplots
fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(8, 10))

for ax, df, title in zip(axs, dataframes, titles):
    # Initialize lists to store counts for each condition
    stability_only = []
    feasibility_only = []
    both_conditions = []

    # Calculate counts for each bin
    for i in range(len(bins) - 1):
        bin_df = df[(df['damping'] >= bins[i]) & (df['damping'] < bins[i + 1])]
        # Check for both conditions first
        both = (bin_df['damping'] < 0.03) & (bin_df['N0'] == 0)
        both_conditions.append(both.sum())
        
        # Check for stability-only condition in the remaining data
        stability = (bin_df['damping'] < 0.03) & ~both
        stability_only.append(stability.sum())
        
        # Check for feasibility-only condition in the remaining data
        feasibility = (bin_df['N0'] == 0) & ~both
        feasibility_only.append(feasibility.sum())

    # Create stacked bar plot
    x = np.arange(len(bins) - 1)
    ax.bar(x, stability_only, width=0.8, label='Stability', color=colors[0], edgecolor='black')
    ax.bar(x, feasibility_only, width=0.8, bottom=stability_only, label='Feasibility', color=colors[1], edgecolor='black')
    ax.bar(x, both_conditions, width=0.8,
           bottom=np.array(stability_only) + np.array(feasibility_only),
           label='Security', color=colors[2], edgecolor='black')

    # Customize each subplot
    ax.set_title(f'{title} Violations')
    ax.set_ylabel('Count')
    ax.set_ylim([0, y_limit])  # Use calculated maximum count
    ax.grid(axis='y')
    ax.tick_params(axis='both', which='major')

# Add x-axis labels to the bottom subplot
axs[-1].set_xticks(x)
axs[-1].set_xticklabels([f"{bins[i] * 100:.2f}-{bins[i + 1] * 100:.2f}%" for i in range(len(bins) - 1)], rotation=45)
axs[-1].set_xlabel('Damping Range (%)')

# Add a shared legend for all subplots
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', title="Condition", bbox_to_anchor=(0.95, 0.95), frameon=False)

# Adjust layout
plt.tight_layout(rect=[0, 0, 0.93, 1])

# Show plot
plt.show()



####################################################################


fn_method_violations = fn_method_ops[['N0', 'damping']] # miss_method
fn_lhc_violations = fn_lhc_ops[['N0', 'damping']] # miss_lhc
fn_imp_violations = fn_imp_ops[['N0', 'damping']] # miss_imp

# Parameters
bins = np.arange(0.0, 0.05, 0.0025)
titles = ['Method', 'LHC', 'Importance']
dataframes = [fn_method_violations, fn_lhc_violations, fn_imp_violations]

# Define professional colors
colors = ["#1b9e77", "#d95f02", "#7570b3"]

# Initialize maximum count tracker
max_count = 0

# Precompute maximum count across all datasets
for df in dataframes:
    for i in range(len(bins) - 1):
        bin_df = df[(df['damping'] >= bins[i]) & (df['damping'] < bins[i + 1])]
        both_conditions = ((bin_df['damping'] > 0.03) & (bin_df['N0'] == 1)).sum()
        stability_only = ((bin_df['damping'] > 0.03) & ~(bin_df['N0'] == 1)).sum()
        feasibility_only = ((bin_df['N0'] == 1) & ~(bin_df['damping'] > 0.03)).sum()
        count = both_conditions + stability_only + feasibility_only
        max_count = max(max_count, count)

# Set consistent y-axis limit with some padding
y_limit = max_count + 10  # Add padding for aesthetics

# Create a figure with subplots
fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(8, 10))

for ax, df, title in zip(axs, dataframes, titles):
    # Initialize lists to store counts for each condition
    stability_only = []
    feasibility_only = []
    both_conditions = []

    # Calculate counts for each bin
    for i in range(len(bins) - 1):
        bin_df = df[(df['damping'] >= bins[i]) & (df['damping'] < bins[i + 1])]
        # Check for both conditions first
        both = (bin_df['damping'] > 0.03) & (bin_df['N0'] == 1)
        both_conditions.append(both.sum())
        
        # Check for stability-only condition in the remaining data
        stability = (bin_df['damping'] > 0.03) & ~both
        stability_only.append(stability.sum())
        
        # Check for feasibility-only condition in the remaining data
        feasibility = (bin_df['N0'] == 1) & ~both
        feasibility_only.append(feasibility.sum())

    # Create stacked bar plot
    x = np.arange(len(bins) - 1)
    ax.bar(x, stability_only, width=0.8, label='Stability', color=colors[0], edgecolor='black')
    ax.bar(x, feasibility_only, width=0.8, bottom=stability_only, label='Feasibility', color=colors[1], edgecolor='black')
    ax.bar(x, both_conditions, width=0.8,
           bottom=np.array(stability_only) + np.array(feasibility_only),
           label='Security', color=colors[2], edgecolor='black')

    # Customize each subplot
    ax.set_title(f'{title} Violations')
    ax.set_ylabel('Count')
    ax.set_ylim([0, y_limit])  # Use calculated maximum count
    ax.grid(axis='y')
    ax.tick_params(axis='both', which='major')

# Add x-axis labels to the bottom subplot
axs[-1].set_xticks(x)
axs[-1].set_xticklabels([f"{bins[i] * 100:.2f}-{bins[i + 1] * 100:.2f}%" for i in range(len(bins) - 1)], rotation=45)
axs[-1].set_xlabel('Damping Range (%)')

# Add a shared legend for all subplots
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', title="Condition", bbox_to_anchor=(0.95, 0.95), frameon=False)

# Adjust layout
plt.tight_layout(rect=[0, 0, 0.93, 1])

# Show plot
plt.show()


########################## combined ###################

# False Positive and False Negative DataFrames
fp_method_violations = fp_method_ops[['N0', 'damping']]
fp_lhc_violations = fp_lhc_ops[['N0', 'damping']]
fp_imp_violations = fp_imp_ops[['N0', 'damping']]

fn_method_violations = fn_method_ops[['N0', 'damping']]
fn_lhc_violations = fn_lhc_ops[['N0', 'damping']]
fn_imp_violations = fn_imp_ops[['N0', 'damping']]

# Parameters
bins = np.arange(0.0, 0.05, 0.0025)
titles = ['Method', 'LHC', 'Importance']
dataframes_fp = [fp_method_violations, fp_lhc_violations, fp_imp_violations]
dataframes_fn = [fn_method_violations, fn_lhc_violations, fn_imp_violations]

# Define professional colors
colors = ["#1b9e77", "#d95f02", "#7570b3"]

# Initialize maximum count tracker
max_count = np.zeros([len(bins),len(dataframes_fp)])
width = 0.8
df_num = 0

# Precompute maximum count across all datasets (False Positives and False Negatives)
for df in dataframes_fp:
    for i in range(len(bins) - 1):
        bin_df = df[(df['damping'] >= bins[i]) & (df['damping'] < bins[i + 1])]
        
        # False Positives conditions
        both_conditions_fp = ((bin_df['damping'] < 0.03) & (bin_df['N0'] == 0)).sum()
        stability_only_fp = ((bin_df['damping'] < 0.03) & ~(bin_df['N0'] == 0)).sum()
        feasibility_only_fp = ((bin_df['N0'] == 0) & ~(bin_df['damping'] < 0.03)).sum()

        count = both_conditions_fp + stability_only_fp + feasibility_only_fp
        
        # Add count to the corresponding bin
        max_count[i,df_num] += count
    df_num += 1

df_num = 0
# Second loop: Handle False Negatives (dataframes_fn)
for df in dataframes_fn:
    for i in range(len(bins) - 1):
        bin_df = df[(df['damping'] >= bins[i]) & (df['damping'] < bins[i + 1])]

        # False Negatives conditions
        both_conditions_fn = ((bin_df['damping'] > 0.03) & (bin_df['N0'] == 1)).sum()
        stability_only_fn = ((bin_df['damping'] > 0.03) & ~(bin_df['N0'] == 1)).sum()
        feasibility_only_fn = ((bin_df['N0'] == 1) & ~(bin_df['damping'] > 0.03)).sum()

        count = both_conditions_fn + stability_only_fn + feasibility_only_fn
        
        # Add count to the corresponding bin
        max_count[i,df_num] += count
    df_num += 1

# Set consistent y-axis limit with some padding
y_limit = max_count.max()
print(max_count)

# Create a figure with subplots
fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(6, 8))

for ax, df_fp, df_fn, title in zip(axs, dataframes_fp, dataframes_fn, titles):
    # Initialize lists to store counts for each condition (False Positives and False Negatives)
    stability_only_fp, feasibility_only_fp, both_conditions_fp = [], [], []
    stability_only_fn, feasibility_only_fn, both_conditions_fn = [], [], []

    # Calculate counts for each bin (combining both False Positive and False Negative)
    for i in range(len(bins) - 1):
        bin_df_fp = df_fp[(df_fp['damping'] >= bins[i]) & (df_fp['damping'] < bins[i + 1])]
        bin_df_fn = df_fn[(df_fn['damping'] >= bins[i]) & (df_fn['damping'] < bins[i + 1])]

        # False Positive Conditions
        both_fp = (bin_df_fp['damping'] < 0.03) & (bin_df_fp['N0'] == 0)
        stability_fp = (bin_df_fp['damping'] < 0.03) & ~both_fp
        feasibility_fp = (bin_df_fp['N0'] == 0) & ~both_fp
        both_conditions_fp.append(both_fp.sum())
        stability_only_fp.append(stability_fp.sum())
        feasibility_only_fp.append(feasibility_fp.sum())

        # False Negative Conditions
        both_fn = (bin_df_fn['damping'] > 0.03) & (bin_df_fn['N0'] == 1)
        stability_fn = (bin_df_fn['damping'] > 0.03) & ~both_fn
        feasibility_fn = (bin_df_fn['N0'] == 1) & ~both_fn
        both_conditions_fn.append(both_fn.sum())
        stability_only_fn.append(stability_fn.sum())
        feasibility_only_fn.append(feasibility_fn.sum())

    # Combine both False Positive and False Negative for each category
    stability_only = np.array(stability_only_fp) + np.array(stability_only_fn)
    feasibility_only = np.array(feasibility_only_fp) + np.array(feasibility_only_fn)
    both_conditions = np.array(both_conditions_fp) + np.array(both_conditions_fn)

    # Create stacked bar plot for each subplot
    x = np.arange(len(bins) - 1)
    x_centered = x - width / 2 

    ax.bar(x_centered, stability_only, width=width, label='Stability', color=colors[0], edgecolor='black')
    ax.bar(x_centered, feasibility_only, width=width, bottom=stability_only, label='Feasibility', color=colors[1], edgecolor='black')
    ax.bar(x_centered, both_conditions, width=width,
           bottom=np.array(stability_only) + np.array(feasibility_only),
           label='Security', color=colors[2], edgecolor='black')

    # Customize each subplot
    ax.set_title(f'{title}')
    ax.set_ylabel('Count')
    ax.set_ylim([0, y_limit*1.05])
    ax.tick_params(axis='both', which='major')

    # Manually set y-ticks (fewer ticks)
    y_ticks = np.arange(0, y_limit * 1.05, step=100)  # Adjust step size
    ax.set_yticks(y_ticks)

    # set grid lines
    ax.grid(True)

    # plot vertical lines
    xticks_positions = [np.where(np.isclose(bins, 0.0))[0][0], np.where(np.isclose(bins, 0.0275))[0][0], np.where(np.isclose(bins, 0.03))[0][0], np.where(np.isclose(bins, 0.0325))[0][0]]
    for pos in x[xticks_positions]:
        ax.axvline(x= pos - width - (1-width)/2, color='grey', linestyle='--', alpha=0.6, linewidth=1)
    

# Add x-axis labels to the bottom subplot
# axs[-1].set_xticks(x_centered)
# axs[-1].set_xticklabels([f"{bins[i] * 100:.2f}-{bins[i + 1] * 100:.2f}" for i in range(len(bins) - 1)], rotation=90)

# only set a few ticks on the x axis
axs[-1].set_xticks(x[xticks_positions]  - width - (1-width)/2 )  
axs[-1].set_xticklabels(['0.0', '2.75', '3.0', '3.25'], rotation=45)
axs[-1].set_xlabel('Damping range ζ [%]')

# Add a general title over the entire figure
# fig.suptitle('Missclassified operating points')

# Add a legend to the first subplot (axs[0])
handles, labels = axs[0].get_legend_handles_labels()
axs[0].legend(handles, labels, loc='upper right', bbox_to_anchor=(1, 1))

# Adjust layout
plt.tight_layout(rect=[0, 0, 0.93, 1])

# Show plot
plt.show()






############## missclassification

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Define bins for damping values
bins = np.arange(0.0, 0.045, 0.0025)
labels = ['Method', 'LHC', 'Importance']
titles = ['Method', 'LHC', 'Importance']

# Set a professional color palette
colors = sns.color_palette("YlOrRd", len(labels) - 1)
colors.append('#1f77b4')  # Dark blue color for the third event

# Combine the dataframes into a list for easier plotting
miss_dataframes = [miss_method, miss_lhc, miss_imp]

# Function to calculate frequencies in each bin
def calculate_frequencies(df, bins):
    counts = []
    for i in range(len(bins) - 1):
        bin_count = ((df['damping'] >= bins[i]) & (df['damping'] < bins[i + 1])).sum()
        counts.append(bin_count)
    return counts

# Calculate frequencies for each miss dataframe
method_miss_counts = calculate_frequencies(miss_method, bins)
lhc_miss_counts = calculate_frequencies(miss_lhc, bins)
imp_miss_counts = calculate_frequencies(miss_imp, bins)

# compute number of test samples in each bin
#total_test_counts = calculate_frequencies(total_test, bins)
total_test_counts = [len(total_test)] * len(imp_miss_counts) # list with total number of test samples in every bin

# Calculate fractions for each method
method_fractions = [m / t if t > 0 else 0 for m, t in zip(method_miss_counts, total_test_counts)]
lhc_fractions = [l / t if t > 0 else 0 for l, t in zip(lhc_miss_counts, total_test_counts)]
imp_fractions = [i / t if t > 0 else 0 for i, t in zip(imp_miss_counts, total_test_counts)]

# Prepare data for the bars
miss_counts_data = [method_miss_counts, lhc_miss_counts, imp_miss_counts]
titles = ['Method', 'LHC', 'Importance']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
labels = ['Method', 'LHC', 'Importance']

# X positions for each bin (we'll shift these positions for each dataset)
x = np.arange(len(bins) - 1)

# Width of each bar in the group
width = 0.25

# Plot setup with academic styling
fig, ax = plt.subplots(figsize=(8, 4))

# Plot each miss dataset as a set of bars grouped at each x instance
for i, (counts, color, label) in enumerate(zip(miss_counts_data, colors, labels)):
    ax.bar(x + i * width - width / 2, counts, width, color=color, edgecolor='black', label=label, alpha=0.8)

# Customize the plot appearance
ax.grid(True, axis='y', linestyle='--', alpha=0.6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_ylim([0, 750])  # Set y-limits, max(max(method_miss_counts), max(lhc_miss_counts), max(imp_miss_counts)) * 1.1

# Set labels and title
ax.set_xlabel('Damping of OP ζ', fontsize=16, labelpad=10)
ax.set_ylabel('# missclassified OPs', fontsize=16, labelpad=10)

# Set x-ticks at 0 and 3%
xticks_positions = [np.where(np.isclose(bins, 0.0275))[0][0], np.where(np.isclose(bins, 0.03))[0][0], np.where(np.isclose(bins, 0.0325))[0][0]]
ax.set_xticks(x[xticks_positions]  - width - width / 2 )  # Add specific x-tick positions (at 0% and 3%) # 
ax.set_xticklabels(['2.75%', '3%', '3.25%'], rotation=30, fontsize=16)

for pos in x[xticks_positions]:
    ax.axvline(x= pos - width - width / 2, color='grey', linestyle=':', alpha=0.6, linewidth=2)

# ax.set_xticks(x + width - width / 2)
# ax.set_xticklabels([f"{bins[i] * 100:.2f}%-{bins[i + 1] * 100:.2f}%" for i in range(len(bins) - 1)], 
#                    rotation=30, ha='right', fontsize=12)

# Manually set y-ticks (fewer ticks)
y_ticks = np.arange(0, max(max(method_miss_counts), max(lhc_miss_counts), max(imp_miss_counts)) * 1.15, step=250)  # Adjust step size
ax.set_yticks(y_ticks)

# Add a box around the plot by enabling the top and right spines
ax.spines['top'].set_visible(True)   # Show the top spine (axis)
ax.spines['right'].set_visible(True) # Show the right spine (axis)
ax.spines['top'].set_linewidth(1)   # Set the thickness of the top spine
ax.spines['right'].set_linewidth(1) # Set the thickness of the right spine
ax.spines['top'].set_color('black')  # Set the color of the top spine
ax.spines['right'].set_color('black') # Set the color of the right spine

# Optionally adjust the bottom and left spines if needed
ax.spines['left'].set_color('black')  # Set the color of the left spine
ax.spines['left'].set_linewidth(1)    # Set the thickness of the left spine
ax.spines['bottom'].set_color('black') # Set the color of the bottom spine
ax.spines['bottom'].set_linewidth(1)   # Set the thickness of the bottom spine

# Add a legend
ax.legend(fontsize=16, loc="upper right", frameon=True)

# Adjust layout for better spacing
plt.tight_layout(rect=[0, 0, 1, 0.97])

# Show plot
plt.show()


#################### second option

# Plot setup: three rows for three datasets
fig, axes = plt.subplots(3, 1, figsize=(8, 12), sharex=True)  # Share x-axis for consistency

# Titles, colors, and labels
titles = ['Method', 'LHC', 'Importance']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
labels = ['Method', 'LHC', 'Importance']

# Plot each miss dataset as a separate bar plot
for i, (ax, counts, color, title) in enumerate(zip(axes, miss_counts_data, colors, titles)):
    # Bar plot for the current dataset
    ax.bar(x, counts, width=0.6, color=color, edgecolor='black', alpha=0.8)

    # Customize the subplot
    ax.set_title(title, fontsize=16, pad=10)
    ax.grid(True, axis='y', linestyle='--', alpha=0.6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='y', which='major', labelsize=14)
    ax.set_ylabel('# Missclassified OPs', fontsize=14, labelpad=10)

    # Set y-ticks with fewer intervals
    y_ticks = np.arange(0, max(counts) * 1.1, step=250)
    ax.set_yticks(y_ticks)

# Customize the shared x-axis for the bottom subplot
axes[-1].set_xlabel('Damping of OP ζ', fontsize=16, labelpad=10)

# Set x-ticks for damping bins
xticks_positions = [
    np.where(np.isclose(bins, 0.0275))[0][0],
    np.where(np.isclose(bins, 0.03))[0][0],
    np.where(np.isclose(bins, 0.0325))[0][0],
]
axes[-1].set_xticks(x[xticks_positions] - 2*width)  # Positions for specific damping values
axes[-1].set_xticklabels(['2.75%', '3%', '3.25%'], rotation=30, fontsize=14)

# Optionally add vertical lines to highlight specific damping bins
for pos in x[xticks_positions]:
    for ax in axes:
        ax.axvline(x=pos - 2*width , color='grey', linestyle=':', alpha=0.6, linewidth=2)

# Adjust layout to avoid overlaps and improve spacing
plt.tight_layout(h_pad=3)  # Horizontal padding between subplots

# Show the final plot
plt.show()





########################## in percentages....

# # Calculate "non-miss" fractions
# method_non_miss = [1 - f for f in method_fractions]
# lhc_non_miss = [1 - f for f in lhc_fractions]
# imp_non_miss = [1 - f for f in imp_fractions]

# # Data for each subplot
# fractions_data = [(method_fractions, method_non_miss), 
#                   (lhc_fractions, lhc_non_miss), 
#                   (imp_fractions, imp_non_miss)]

# titles = ['Method Miss Fractions', 'LHC Miss Fractions', 'IMP Miss Fractions']
# colors_miss = ['red', 'blue', 'green']
# colors_non_miss = ['pink', 'lightblue', 'lightgreen']
# labels = ['Misses', 'Non-misses']

# # Plot setup
# fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(12, 12), sharex=True)

# # X positions for each bin
# x = np.arange(len(bins) - 1)
# width = 0.8

# # Function to plot each dataset as a normalized stacked bar chart
# def plot_normalized_subplot(ax, fractions, title, x, width, color_miss, color_non_miss):
#     miss_fractions, non_miss_fractions = fractions
#     ax.bar(x, miss_fractions, width, color=color_miss, edgecolor='black', label='Misses')
#     ax.bar(x, non_miss_fractions, width, bottom=miss_fractions, color=color_non_miss, edgecolor='black', label='Non-misses')
    
#     # Customize subplot appearance
#     ax.grid(True, axis='y', linestyle='--', alpha=0.7)
#     ax.spines['top'].set_visible(False)
#     ax.spines['right'].set_visible(False)
#     ax.tick_params(axis='both', which='major', labelsize=12)
#     ax.set_ylim([0, 1])  # Normalize each bin to 100%
#     ax.legend(fontsize=12, loc="upper right", frameon=False)
#     ax.set_title(title, fontsize=14, fontweight='bold')

# # Plot each dataset in its subplot
# for ax, fractions, title, color_miss, color_non_miss in zip(axs, fractions_data, titles, colors_miss, colors_non_miss):
#     plot_normalized_subplot(ax, fractions, title, x, width, color_miss, color_non_miss)

# # Set x-axis and y-axis labels
# axs[-1].set_xlabel('Damping Range (%)', fontsize=12, fontweight='bold', labelpad=10)
# axs[1].set_ylabel('Fraction of Operations', fontsize=12, fontweight='bold', labelpad=10)

# # Improve x-axis labels format
# axs[-1].set_xticks(x)
# axs[-1].set_xticklabels([f"{bins[i] * 100:.1f}%-{bins[i + 1] * 100:.1f}%" for i in range(len(bins) - 1)], 
#                         rotation=30, ha='right', fontsize=12)

# plt.tight_layout(rect=[0, 0, 1, 0.97])  # Adjust layout to accommodate title
# plt.show()






# from sklearn.model_selection import GridSearchCV

# # Hyperparameter to fine tune
# param_grid = {
#     'max_depth': range(1, 10, 1),
#     'min_samples_leaf': range(1, 20, 2),
#     'min_samples_split': range(2, 20, 2),
#     'criterion': ["entropy", "gini"]
# }
# # Decision tree classifier
# tree = DecisionTreeClassifier(random_state=1)
# # GridSearchCV
# grid_search = GridSearchCV(estimator=tree, param_grid=param_grid, 
#                            cv=5, verbose=True)
# grid_search.fit(X_train, y_train)

# # Best score and estimator
# print("best accuracy", grid_search.best_score_)
# print(grid_search.best_estimator_)



# # plot confusion matrix
# confusion_mat_method = confusion_matrix(y_test_method, y_pred_method, labels = [0,1])
# sns.heatmap(confusion_mat_method, annot=True, cmap='Paired', 
#             cbar=False, fmt="d", xticklabels=[
#             'Secure', 'Insecure'], yticklabels=['Secure', 'Insecure'])

# confusion_mat_lhc = confusion_matrix(y_test_lhc, y_pred_lhc, labels = [0,1])
# sns.heatmap(confusion_mat_lhc, annot=True, cmap='Paired', 
#             cbar=False, fmt="d", xticklabels=[
#             'Secure', 'Insecure'], yticklabels=['Secure', 'Insecure'])

# confusion_mat_imp = confusion_matrix(y_test_imp, y_pred_imp, labels = [0,1])
# sns.heatmap(confusion_mat_imp, annot=True, cmap='Paired', 
#             cbar=False, fmt="d", xticklabels=[
#             'Secure', 'Insecure'], yticklabels=['Secure', 'Insecure'])



# # Visualize the decision tree
# plt.figure(figsize=(20,10))
# tree.plot_tree(clf_method, filled=True, feature_names=[f'feature_{i}' for i in range(X_method.shape[1])], class_names=['insecure', 'secure'])
# plt.show()

# plt.figure(figsize=(20,10))
# tree.plot_tree(clf_lhc, filled=True, feature_names=[f'feature_{i}' for i in range(X_lhc.shape[1])], class_names=['insecure', 'secure'])
# plt.show()

# plt.figure(figsize=(20,10))
# tree.plot_tree(clf_imp, filled=True, feature_names=[f'feature_{i}' for i in range(X_imp.shape[1])], class_names=['insecure', 'secure'])
# plt.show()




# # show classification report
# classification_method = classification_report(y_test_method, y_pred_method) # shows precision, recall, f1-score 
# print("Classification Report:\n", classification_method)

# classification_lhc = classification_report(y_test_lhc, y_pred_lhc) # shows precision, recall, f1-score 
# print("Classification Report:\n", classification_lhc)

# classification_imp = classification_report(y_test_imp, y_pred_imp) # shows precision, recall, f1-score 
# print("Classification Report:\n", classification_imp)




# Plot the top 10 features
# importance_df_method.head(10).plot(kind='bar', x='Feature', y='Importance', legend=False)
# plt.show()

# importance_df_lhc.head(10).plot(kind='bar', x='Feature', y='Importance', legend=False)
# plt.show()

# importance_df_imp.head(10).plot(kind='bar', x='Feature', y='Importance', legend=False)
# plt.show()