import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import classification_report, accuracy_score, f1_score, confusion_matrix, precision_score, recall_score, det_curve
from sklearn import tree
import matplotlib.pyplot as plt
import seaborn as sns
import os

"""

Al-amin had three test datasets: historical, generic and rare.
I can maybe use the test sets of the three methods, and additionally make my own rare test set.

The rare test set could maybe be all samples in the three test sets that are 2% < damping < 4%.
Al-amin: rare data is data deviating from the historical distribution. He picked random samples from
his unified sampling approach. Maybe I can bundle all test samples, draw a distribution over it,
and randomply pick out points which are 'rare'. The goal is to see how well the method generalizes
to 'rare' points.

We could also try to see if we can prove if the HIC region indeed is important. Check where the 
missclassification occurs. Is this close to the HIC region, or far from the HIC region?
We could maybe make a dataset that consists only of points in this HIC region.

A lot of false positives for the method are at samples with -1 damping. This is because a lot of 
data from LHC and importance is centered at this point, and the training data of method contains
no samples with -1 damping. For the rest, the method contains almost only samples around HIC region
and classifies a lot of points along the whole spectrum correctly. The unbalance around -1 damping 
messes it up. Also, maybe HIC region should be made larger.


"""

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

    FPR = FP/(FP+TN)

    return FPR, FP_list, FN_list

def plot_damping_histograms(x_lb, bar_width, damping_data, training_data, FP_data, FN_data, title):
    """
    Function to plot the damping of the missclassified points.
    
    """
    # Define bin edges with steps of 0.0025, starting from 0
    bins = np.arange(x_lb, max(damping_data.max(), FP_data.max(), FN_data.max()) + bar_width, bar_width)

    # Create a figure with three subplots, arranged horizontally (1 row, 3 columns)
    fig, axs = plt.subplots(1, 4, figsize=(20, 5))

    max_y = 0

    # Plot the histogram for 'additional_data'
    axs[0].hist(training_data, bins=bins, alpha=0.7, color='purple')
    axs[0].set_title('Training Data')
    axs[0].set_xlabel('Damping Values')
    axs[0].set_ylabel('Amount of OPs (log scale)')
    axs[0].set_yscale('log')  # Set y-axis to logarithmic scale
    max_y = max(max_y, axs[0].get_ylim()[1]) 

    # Plot the histogram for 'damping_data'
    axs[1].hist(damping_data, bins=bins, alpha=0.7, color='blue')
    axs[1].set_title('Testing Data')
    axs[1].set_xlabel('Damping Values')
    axs[1].set_yscale('log')  # Set y-axis to logarithmic scale
    max_y = max(max_y, axs[1].get_ylim()[1])

    # Plot the histogram for 'FP_data'
    axs[2].hist(FP_data, bins=bins, alpha=0.7, color='red')
    axs[2].set_title('FP Data')
    axs[2].set_xlabel('Damping Values')
    axs[2].set_yscale('log')  # Set y-axis to logarithmic scale
    max_y = max(max_y, axs[2].get_ylim()[1])

    # Plot the histogram for 'FN_data'
    axs[3].hist(FN_data, bins=bins, alpha=0.7, color='green')
    axs[3].set_title('FN Data')
    axs[3].set_xlabel('Damping Values')
    axs[3].set_yscale('log')  # Set y-axis to logarithmic scale
    max_y = max(max_y, axs[3].get_ylim()[1])

    # Set the same y-limits for all plots
    for ax in axs:
        #ax.set_ylim(1e-1, max_y)  # You can adjust the lower limit (1e-1) as needed
        ax.set_ylim(1e-1, 15000)

    # Add a title to the entire figure
    fig.suptitle(title, fontsize=16)

    # Add grid lines to all subplots
    for ax in axs:
        ax.grid(True)

    # Show the plot
    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust spacing, leaving room for the title
    plt.show()

# specify test set: 'all', 'HIC', 'noHIC'
test_index = "all"

# Specify the directory and file name
directory = "C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/output/case39/datasets/"
flow_name_method = '39bus_method_flows.csv' # DT method data
flow_name_lhc = '39bus_lhc_flows.csv'
flow_name_imp = '39bus_imp_flows.csv' 

file_name_method = '39bus_method_ops.csv' # DT method data
file_name_lhc = '39bus_lhc_ops.csv'
file_name_imp = '39bus_imp_ops.csv' 

file_name_result = 'dt_result.csv'

# Create the full file path
flow_path_method = os.path.join(directory, flow_name_method)
flow_path_lhc = os.path.join(directory, flow_name_lhc)
flow_path_imp = os.path.join(directory, flow_name_imp)

file_path_method = os.path.join(directory, file_name_method)
file_path_lhc = os.path.join(directory, file_name_lhc)
file_path_imp = os.path.join(directory, file_name_imp)

file_path_result = os.path.join(directory, file_name_result)

# Load the OP dataset
op_data_method = pd.read_csv(file_path_method, sep = ';')
op_index_method = np.ones(len(op_data_method), dtype=bool)
op_index_method_HIC = (op_data_method['damping'] >= 0.0275) & (op_data_method['damping'] <= 0.0325)
op_index_method_not_HIC = (op_data_method['damping'] < 0.0275) | (op_data_method['damping'] > 0.0325)

op_data_lhc = pd.read_csv(file_path_lhc, sep = ';')
op_index_lhc = np.ones(len(op_data_lhc), dtype=bool)
op_index_lhc_HIC = (op_data_lhc['damping'] >= 0.0275) & (op_data_lhc['damping'] <= 0.0325)
op_index_lhc_not_HIC = (op_data_lhc['damping'] < 0.0275) | (op_data_lhc['damping'] > 0.0325)

op_data_imp = pd.read_csv(file_path_imp, sep = ';')
op_index_imp = np.ones(len(op_data_imp), dtype=bool)
op_index_imp_HIC = (op_data_imp['damping'] >= 0.0275) & (op_data_imp['damping'] <= 0.0325)
op_index_imp_not_HIC = (op_data_imp['damping'] < 0.0275) | (op_data_imp['damping'] > 0.0325)

if test_index == 'all':
    op_index_method = op_index_method
    op_index_lhc = op_index_lhc
    op_index_imp = op_index_imp
elif test_index == 'HIC': 
    op_index_method = op_index_method_HIC
    op_index_lhc = op_index_lhc_HIC
    op_index_imp = op_index_imp_HIC
elif test_index == 'noHIC':
    op_index_method = op_index_method_not_HIC
    op_index_lhc = op_index_lhc_not_HIC
    op_index_imp = op_index_imp_not_HIC


# Load the flow dataset
data_method = pd.read_csv(flow_path_method, sep = ';')
X_method = data_method.drop(columns = ['feasible', 'stable'], axis=1) # drop columns you won't use like ['Feas', 'N1 ....]
y_method_stable = data_method['stable'] 
y_method_feasible = data_method['feasible']
y_method = y_method_stable*y_method_feasible
X_train_method, X_test_method, y_train_method, y_test_method, method_index_train, method_index_test = train_test_split(X_method, y_method, op_index_method, test_size=0.25, random_state=42)
X_test_method = X_test_method[method_index_test]
y_test_method = y_test_method[method_index_test]

data_lhc = pd.read_csv(flow_path_lhc, sep = ';')
X_lhc = data_lhc.drop(columns = ['feasible', 'stable'], axis=1) # drop columns you won't use like ['Feas', 'N1 ....]
y_lhc_stable = data_lhc['stable'] 
y_lhc_feasible = data_lhc['feasible'] 
y_lhc = y_lhc_stable*y_lhc_feasible
X_train_lhc, X_test_lhc, y_train_lhc, y_test_lhc, lhc_index_train, lhc_index_test = train_test_split(X_lhc, y_lhc, op_index_lhc, test_size=0.25, random_state=42)
X_test_lhc = X_test_lhc[lhc_index_test]
y_test_lhc = y_test_lhc[lhc_index_test]

data_imp = pd.read_csv(flow_path_imp, sep = ';')
X_imp = data_imp.drop(columns = ['feasible', 'stable'], axis=1) # drop columns you won't use like ['Feas', 'N1 ....]
y_imp_stable = data_imp['stable'] 
y_imp_feasible = data_imp['feasible'] 
y_imp = y_imp_stable*y_imp_feasible
X_train_imp, X_test_imp, y_train_imp, y_test_imp, imp_index_train, imp_index_test = train_test_split(X_imp, y_imp, op_index_imp, test_size=0.25, random_state=42)
X_test_imp = X_test_imp[imp_index_test]
y_test_imp = y_test_imp[imp_index_test]

# Clear the CSV file before starting the run (optional)
if os.path.exists(file_path_result):
    with open(file_path_result, 'w') as f:
        f.write('')  # Write an empty string to clear the file

# iteratively reduce training set size
reduce_train = 1
initial_fraction = 1.0  # Start with the full dataset
reduction_step = 0.0    # Reduce by this fraction in each iteration

for i in range(reduce_train):
    # Calculate the new fraction for the current iteration
    current_fraction = initial_fraction - (i * reduction_step)
    
    # Ensure the fraction doesn't go below a certain limit
    current_fraction = max(current_fraction, 0.1)  # Keep at least 10% of the data

    # Randomly sample a subset of the training sets
    X_train_method_subset = X_train_method.sample(frac=current_fraction, random_state=42)
    y_train_method_subset = y_train_method.loc[X_train_method_subset.index]  # Ensure labels match
    
    X_train_lhc_subset = X_train_lhc.sample(frac=current_fraction, random_state=42)
    y_train_lhc_subset = y_train_lhc.loc[X_train_lhc_subset.index]  # Ensure labels match
    
    X_train_imp_subset = X_train_imp.sample(frac=current_fraction, random_state=42)
    y_train_imp_subset = y_train_imp.loc[X_train_imp_subset.index]  # Ensure labels match

    # Initialize the DecisionTreeClassifier with gini impurity and max depth of 5. clf = classifier
    clf_method = DecisionTreeClassifier(criterion='gini', max_depth=5, random_state=42) # splitter='best' or 'random' , ccp_alpha=0.01 for tree pruning
    clf_lhc = DecisionTreeClassifier(criterion='gini', max_depth=5, random_state=42)
    clf_imp = DecisionTreeClassifier(criterion='gini', max_depth=5, random_state=42)

    # Perform 10-fold cross-validation and compute accuracy for each fold
    cv_scores_method = cross_val_score(clf_method, X_method, y_method, cv=10, scoring='accuracy')
    cv_scores_lhc = cross_val_score(clf_lhc, X_lhc, y_lhc, cv=10, scoring='accuracy')
    cv_scores_imp = cross_val_score(clf_imp, X_imp, y_imp, cv=10, scoring='accuracy')

    # # Output the results
    # print("Cross-validation scores:", cv_scores)
    # print("Mean cross-validation accuracy:", np.mean(cv_scores))
    # print("Standard deviation of cross-validation accuracy:", np.std(cv_scores))

    # Train the classifier on the reduced training sets
    clf_method.fit(X_train_method_subset, y_train_method_subset)
    clf_lhc.fit(X_train_lhc_subset, y_train_lhc_subset)
    clf_imp.fit(X_train_imp_subset, y_train_imp_subset)

    # Make predictions on the test sets
    y_pred_method = clf_method.predict(X_test_method)
    y_pred_method_lhc = clf_method.predict(X_test_lhc)
    y_pred_method_imp = clf_method.predict(X_test_imp)
    # y_prob_method = clf_method.predict_proba(X_test_method) # due to early stopping, there is a probability of being wrong. Not all leafs are pure. 

    y_pred_lhc = clf_lhc.predict(X_test_lhc)
    y_pred_lhc_method = clf_lhc.predict(X_test_method)
    y_pred_lhc_imp = clf_lhc.predict(X_test_imp)
    # y_prob_lhc = clf_lhc.predict_proba(X_test_lhc)

    y_pred_imp = clf_imp.predict(X_test_imp)
    y_pred_imp_method = clf_imp.predict(X_test_method)
    y_pred_imp_lhc = clf_imp.predict(X_test_lhc)
    # y_prob_imp = clf_imp.predict_proba(X_test_imp)

    # Evaluate the classifier accuracy
    accuracy_method = accuracy_score(y_test_method, y_pred_method)
    accuracy_method_lhc = accuracy_score(y_test_lhc, y_pred_method_lhc)
    accuracy_method_imp = accuracy_score(y_test_imp, y_pred_method_imp)
    print("Test accuracy train  method, test method:", accuracy_method)
    print("Test accuracy train  method, test lhc:", accuracy_method_lhc)
    print("Test accuracy train  method, test imp:", accuracy_method_imp)

    accuracy_lhc = accuracy_score(y_test_lhc, y_pred_lhc)
    accuracy_lhc_method = accuracy_score(y_test_method, y_pred_lhc_method)
    accuracy_lhc_imp = accuracy_score(y_test_imp, y_pred_lhc_imp)
    print("Test accuracy train lhc, test lhc:", accuracy_lhc)
    print("Test accuracy train lhc, test method:", accuracy_lhc_method)
    print("Test accuracy train lhc, test imp:", accuracy_lhc_imp)

    accuracy_imp = accuracy_score(y_test_imp, y_pred_imp)
    accuracy_imp_method = accuracy_score(y_test_method, y_pred_imp_method)
    accuracy_imp_lhc = accuracy_score(y_test_lhc, y_pred_imp_lhc)
    print("Test accuracy train imp, test imp:", accuracy_imp)
    print("Test accuracy train imp, test method:", accuracy_imp_method)
    print("Test accuracy train imp, test lhc:", accuracy_imp_lhc)


    # Evaluate the classifier precision
    precision_method = precision_score(y_test_method, y_pred_method)
    precision_method_lhc = precision_score(y_test_lhc, y_pred_method_lhc)
    precision_method_imp = precision_score(y_test_imp, y_pred_method_imp)
    print("Test precision train  method, test method:", precision_method)
    print("Test precision train  method, test lhc:", precision_method_lhc)
    print("Test precision train  method, test imp:", precision_method_imp)

    precision_lhc = precision_score(y_test_lhc, y_pred_lhc)
    precision_lhc_method = precision_score(y_test_method, y_pred_lhc_method)
    precision_lhc_imp = precision_score(y_test_imp, y_pred_lhc_imp)
    print("Test precision train lhc, test lhc:", precision_lhc)
    print("Test precision train lhc, test method:", precision_lhc_method)
    print("Test precision train lhc, test imp:", precision_lhc_imp)

    precision_imp = precision_score(y_test_imp, y_pred_imp)
    precision_imp_method = precision_score(y_test_method, y_pred_imp_method)
    precision_imp_lhc = precision_score(y_test_lhc, y_pred_imp_lhc)
    print("Test precision train imp, test imp:", precision_imp)
    print("Test precision train imp, test method:", precision_imp_method)
    print("Test precision train imp, test lhc:", precision_imp_lhc)


    # Evaluate the classifier recall
    recall_method = recall_score(y_test_method, y_pred_method)
    recall_method_lhc = recall_score(y_test_lhc, y_pred_method_lhc)
    recall_method_imp = recall_score(y_test_imp, y_pred_method_imp)
    print("Test recall train  method, test method:", recall_method)
    print("Test recall train  method, test lhc:", recall_method_lhc)
    print("Test recall train  method, test imp:", recall_method_imp)

    recall_lhc = recall_score(y_test_lhc, y_pred_lhc)
    recall_lhc_method = recall_score(y_test_method, y_pred_lhc_method)
    recall_lhc_imp = recall_score(y_test_imp, y_pred_lhc_imp)
    print("Test recall train lhc, test lhc:", recall_lhc)
    print("Test recall train lhc, test method:", recall_lhc_method)
    print("Test recall train lhc, test imp:", recall_lhc_imp)

    recall_imp = recall_score(y_test_imp, y_pred_imp)
    recall_imp_method = recall_score(y_test_method, y_pred_imp_method)
    recall_imp_lhc = recall_score(y_test_lhc, y_pred_imp_lhc)
    print("Test recall train imp, test imp:", recall_imp)
    print("Test recall train imp, test method:", recall_imp_method)
    print("Test recall train imp, test lhc:", recall_imp_lhc)


    # Evaluate the classifier F1 score
    f1_method = f1_score(y_test_method, y_pred_method)
    f1_method_lhc = f1_score(y_test_lhc, y_pred_method_lhc)
    f1_method_imp = f1_score(y_test_imp, y_pred_method_imp)
    print("Test f1_score train  method, test method:", f1_method)
    print("Test f1_score train  method, test lhc:", f1_method_lhc)
    print("Test f1_score train  method, test imp:", f1_method_imp)

    f1_lhc = f1_score(y_test_lhc, y_pred_lhc)
    f1_lhc_method = f1_score(y_test_method, y_pred_lhc_method)
    f1_lhc_imp = f1_score(y_test_imp, y_pred_lhc_imp)
    print("Test f1_score train lhc, test lhc:", f1_lhc)
    print("Test f1_score train lhc, test method:", f1_lhc_method)
    print("Test f1_score train lhc, test imp:", f1_lhc_imp)

    f1_imp = f1_score(y_test_imp, y_pred_imp)
    f1_imp_method = f1_score(y_test_method, y_pred_imp_method)
    f1_imp_lhc = f1_score(y_test_lhc, y_pred_imp_lhc)
    print("Test f1_score train imp, test imp:", f1_imp)
    print("Test f1_score train imp, test method:", f1_imp_method)
    print("Test f1_score train imp, test lhc:", f1_imp_lhc)


    # Evaluate the classifier false positive rate. det_curve = fpr, fnr, threshold
    fpr_method, FP_list_method, FN_list_method = fpr_score(y_test_method, y_pred_method)
    fpr_method_lhc, FP_list_method_lhc, FN_list_method_lhc = fpr_score(y_test_lhc, y_pred_method_lhc)
    fpr_method_imp, FP_list_method_imp, FN_list_method_imp = fpr_score(y_test_imp, y_pred_method_imp)
    print("Test fpr train  method, test method:", fpr_method)
    print("Test fpr train  method, test lhc:", fpr_method_lhc)
    print("Test fpr train  method, test imp:", fpr_method_imp)

    fpr_lhc, FP_list_lhc, FN_list_lhc = fpr_score(y_test_lhc, y_pred_lhc)
    fpr_lhc_method, FP_list_lhc_method, FN_list_lhc_method = fpr_score(y_test_method, y_pred_lhc_method)
    fpr_lhc_imp, FP_list_lhc_imp, FN_list_lhc_imp = fpr_score(y_test_imp, y_pred_lhc_imp)
    print("Test fpr train lhc, test lhc:", fpr_lhc)
    print("Test fpr train lhc, test method:", fpr_lhc_method)
    print("Test fpr train lhc, test imp:", fpr_lhc_imp)

    fpr_imp, FP_list_imp, FN_list_imp = fpr_score(y_test_imp, y_pred_imp)
    fpr_imp_method, FP_list_imp_method, FN_list_imp_method = fpr_score(y_test_method, y_pred_imp_method)
    fpr_imp_lhc, FP_list_imp_lhc, FN_list_imp_lhc = fpr_score(y_test_lhc, y_pred_imp_lhc)
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
        'method test accuracy': [accuracy_method, accuracy_lhc_method, accuracy_imp_method],
        'lhc test accuracy': [accuracy_method_lhc, accuracy_lhc, accuracy_imp_lhc],
        'imp test accuracy': [accuracy_method_imp, accuracy_lhc_imp, accuracy_imp],
        'space': ['', '', ''],
        'method test recall': [recall_method, recall_lhc_method, recall_imp_method],
        'lhc test recall': [recall_method_lhc, recall_lhc, recall_imp_lhc],
        'imp test recall': [recall_method_imp, recall_lhc_imp, recall_imp],
        'space1': ['', '', ''],
        'method test precision': [precision_method, precision_lhc_method, precision_imp_method],
        'lhc test precision': [precision_method_lhc, precision_lhc, precision_imp_lhc],
        'imp test precision': [precision_method_imp, precision_lhc_imp, precision_imp],
        'space2': ['', '', ''],
        'method test f1': [f1_method, f1_lhc_method, f1_imp_method],
        'lhc test f1': [f1_method_lhc, f1_lhc, f1_imp_lhc],
        'imp test f1': [f1_method_imp, f1_lhc_imp, f1_imp],
        'space3': ['', '', ''],
        'method test fpr': [fpr_method, fpr_lhc_method, fpr_imp_method],
        'lhc test fpr': [fpr_method_lhc, fpr_lhc, fpr_imp_lhc],
        'imp test fpr': [fpr_method_imp, fpr_lhc_imp, fpr_imp],

    }
        
    df_result = pd.DataFrame.from_dict(data_result)

    # Append a newline before writing new results
    with open(file_path_result, 'a', newline='') as f:
        if f.tell() > 0:  # Check if the file already has content
            f.write('\n')  # Add a blank line before appending the new DataFrame
        df_result.to_csv(f, sep=';', header=f.tell() == 0, index=False)





" determine distance to HIC region of missclassified points "

# get damping OPs
damping_data_method = op_data_method["damping"]
damping_data_lhc = op_data_lhc["damping"]
damping_data_imp = op_data_imp["damping"]

damping_test_data_method = damping_data_method[y_test_method.index]
damping_test_data_lhc = damping_data_lhc[y_test_lhc.index]
damping_test_data_imp = damping_data_imp[y_test_imp.index]

damping_train_data_method = damping_data_method[y_train_method.index]
damping_train_data_lhc = damping_data_lhc[y_train_lhc.index]
damping_train_data_imp = damping_data_imp[y_train_imp.index]

damping_all_test_data = pd.concat([damping_test_data_method, damping_test_data_lhc, damping_test_data_imp])

# get damping of FPs and FNs
FP_damping_method = pd.concat([damping_data_method[FP_list_method], damping_data_lhc[FP_list_method_lhc], damping_data_imp[FP_list_method_imp]])
FP_damping_lhc = pd.concat([damping_data_lhc[FP_list_lhc], damping_data_method[FP_list_lhc_method], damping_data_imp[FP_list_lhc_imp]])
FP_damping_imp = pd.concat([damping_data_imp[FP_list_imp], damping_data_method[FP_list_imp_method], damping_data_lhc[FP_list_imp_lhc]])

FN_damping_method = pd.concat([damping_data_method[FN_list_method], damping_data_lhc[FN_list_method_lhc], damping_data_imp[FN_list_method_imp]])
FN_damping_lhc = pd.concat([damping_data_lhc[FN_list_lhc], damping_data_method[FN_list_lhc_method], damping_data_imp[FN_list_lhc_imp]])
FN_damping_imp = pd.concat([damping_data_imp[FN_list_imp], damping_data_method[FN_list_imp_method], damping_data_lhc[FN_list_imp_lhc]])

# plot the damping of the missclassified OPs
x_lb = 0
bar_width = 0.005
plot_damping_histograms(x_lb, bar_width, damping_all_test_data, damping_train_data_method, FP_damping_method, FN_damping_method, "Proposed method samples - damping vs number of OPs")
plot_damping_histograms(x_lb, bar_width, damping_all_test_data, damping_train_data_lhc, FP_damping_lhc, FN_damping_lhc, "LHC samples - damping vs number of OPs")
plot_damping_histograms(x_lb, bar_width, damping_all_test_data, damping_train_data_imp, FP_damping_imp, FN_damping_imp, "Importance samples - damping vs number of OPs")



##############################################################################

fp_method_ops = pd.concat([op_data_method.iloc[FP_list_method], op_data_lhc.iloc[FP_list_method_lhc], op_data_imp.iloc[FP_list_method_imp]])
fp_lhc_ops = pd.concat([op_data_lhc.iloc[FP_list_lhc], op_data_method.iloc[FP_list_lhc_method], op_data_imp.iloc[FP_list_lhc_imp]])
fp_imp_ops = pd.concat([op_data_imp.iloc[FP_list_imp], op_data_method.iloc[FP_list_imp_method], op_data_lhc.iloc[FP_list_imp_lhc]])

fn_method_ops = pd.concat([op_data_method.iloc[FN_list_method], op_data_lhc.iloc[FN_list_method_lhc], op_data_imp.iloc[FN_list_method_imp]])
fn_lhc_ops = pd.concat([op_data_lhc.iloc[FN_list_lhc], op_data_method.iloc[FN_list_lhc_method], op_data_imp.iloc[FN_list_lhc_imp]])
fn_imp_ops = pd.concat([op_data_imp.iloc[FN_list_imp], op_data_method.iloc[FN_list_imp_method], op_data_lhc.iloc[FN_list_imp_lhc]])

miss_method = pd.concat([fp_method_ops, fn_method_ops])
miss_lhc = pd.concat([fp_lhc_ops, fn_lhc_ops])
miss_imp = pd.concat([fp_imp_ops, fn_imp_ops])

fp_method_cf = fp_method_ops[['N0', 'damping']] # miss_method
fp_lhc_cf = fp_lhc_ops[['N0', 'damping']] # miss_lhc
fp_imp_cf = fp_imp_ops[['N0', 'damping']] # miss_imp

# Calculate frequencies for the three conditions in each dataframe
def calculate_event_frequencies(df):
    event1 = (df['N0'] == 0).sum()
    event2 = (df['damping'] < 0.03).sum()
    event3 = ((df['N0'] == 0) & (df['damping'] < 0.03)).sum()
    return [event1, event2, event3]

# Get frequencies for each dataframe
method_counts = calculate_event_frequencies(fp_method_cf)
lhc_counts = calculate_event_frequencies(fp_lhc_cf)
imp_counts = calculate_event_frequencies(fp_imp_cf)

# Prepare data for plotting
data = [method_counts, lhc_counts, imp_counts]

# Labels and colors
labels = ['N0 == 0', 'damping < 0.03', 'N0 == 0 & damping < 0.03']
df_labels = ['Method', 'LHC', 'IMP']
colors = ['#f4a261', '#e76f51', '#e63946']

# Bar plot with spacing between groups
fig, ax = plt.subplots(figsize=(12, 7))
width = 0.25  # Width of individual bars
x = np.arange(len(labels))  # x positions for groups of bars

# Plotting bars with edges
for i, (d, color) in enumerate(zip(data, colors)):
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


####################################################################

# Assuming fp_method_ops, fp_lhc_ops, and fp_imp_ops are already defined
fp_method_violations = fp_method_ops[['N0P', 'N0Q', 'N0L', 'N0OV', 'N0UV', 'damping']] # miss_method
fp_lhc_violations = fp_lhc_ops[['N0P', 'N0Q', 'N0L', 'N0OV', 'N0UV', 'damping']] # miss_lhc
fp_imp_violations = fp_imp_ops[['N0P', 'N0Q', 'N0L', 'N0OV', 'N0UV', 'damping']] # miss_imp

# Define bins for damping values
bins = np.arange(-0.07, 0.07, 0.005)
labels = ['N0P', 'N0Q', 'N0L', 'N0OV', 'N0UV', 'All Zero']  # Adding 'All Zero' for the extra bar

# Create a figure and axes for the subplots
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(16, 24))  # Adjust the figure size as needed

# Loop through the DataFrames and axes to create the plots
dataframes = [fp_method_violations, fp_lhc_violations, fp_imp_violations]
titles = ['Method Violations', 'LHC Violations', 'Importance Violations']

for ax, df, title in zip(axs, dataframes, titles):
    # Initialize a DataFrame to store the frequency counts per bin
    bin_counts = pd.DataFrame(0, index=np.arange(len(bins) - 1), columns=labels)

    # Calculate frequencies of non-zero occurrences in each bin
    for i in range(len(bins) - 1):
        # Filter rows where 'damping' falls into the current bin range
        bin_df = df[(df['damping'] >= bins[i]) & (df['damping'] < bins[i + 1])]
        # Count non-zero occurrences for each column in the current bin
        for label in labels[:-1]:  # Count for N0P, N0Q, N0L, N0OV, N0UV
            bin_counts.loc[i, label] = (bin_df[label] != 0).sum()
        
        # Count occurrences where all specified columns are zero and damping is less than 0.03
        all_zero_count = ((bin_df[['N0P', 'N0Q', 'N0L', 'N0OV', 'N0UV']] == 0).all(axis=1) & 
                          (bin_df['damping'] < 0.03)).sum()
        bin_counts.loc[i, 'All Zero'] = all_zero_count

    # Plotting
    width = 0.13  # Width of individual bars within each bin

    # Define colors using a color palette
    colors = sns.color_palette("YlOrRd", len(labels) - 1)  # Exclude the 'All Zero' bar from this palette
    colors.append('#1f77b4')  # Append dark blue color for 'All Zero' bar

    # X positions for each set of bars, centered on the bin edges
    x = np.arange(len(bins) - 1)

    # Plot each label as a separate set of bars within each bin
    for j, (label, color) in enumerate(zip(labels, colors)):
        ax.bar(x + j * width - (width * len(labels) / 2), bin_counts[label], width, 
               label=label, color=color, edgecolor='black')

    # Customize plot
    ax.set_xticks(x)  # Set x-ticks centered on the bins
    ax.set_xticklabels([f"{bins[i] * 100:.1f}%-{bins[i + 1] * 100:.1f}%" for i in range(len(bins) - 1)], rotation=30, fontsize=8)
    ax.set_xlabel('Damping Range (%)', fontsize=12, labelpad=10)  # Updated the x-label to indicate percentage
    ax.set_ylabel('Frequency of Violations', fontsize=12)
    ax.set_title(f'False Positive Cases - {title}', fontsize=14, pad=15)
    ax.legend(title="Violation Type", title_fontsize=12, fontsize=10)
    ax.grid(True, axis='y', linestyle='--', alpha=0.6)

    # Style adjustments
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.show()




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
