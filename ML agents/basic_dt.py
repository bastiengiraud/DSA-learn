import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import classification_report, accuracy_score, f1_score, confusion_matrix, precision_score, recall_score
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


"""

# Specify the directory and file name
directory = "C:/Users/bagir/OneDrive - Danmarks Tekniske Universitet/Dokumenter/1) Projects/2) Datasets/2) Datasets code/output/case39/datasets/"
file_name_method = 'flows2.csv' # DT method data
file_name_lhc = 'lhc_flows.csv'
file_name_imp = 'imp_flows.csv' 

file_name_result = 'dt_result2.csv' 

# Create the full file path
file_path_method = os.path.join(directory, file_name_method)
file_path_lhc = os.path.join(directory, file_name_lhc)
file_path_imp = os.path.join(directory, file_name_imp)

file_path_result = os.path.join(directory, file_name_result)

# Load the dataset
data_method = pd.read_csv(file_path_method, sep = ';')
X_method = data_method.drop(columns = ['feasible', 'stable'], axis=1) # drop columns you won't use like ['Feas', 'N1 ....]
y_method_stable = data_method['stable'] 
y_method_feasible = data_method['feasible']
y_method = y_method_stable*y_method_feasible
X_train_method, X_test_method, y_train_method, y_test_method = train_test_split(X_method, y_method, test_size=0.25, random_state=42)

data_lhc = pd.read_csv(file_path_lhc, sep = ';')
X_lhc = data_lhc.drop(columns = ['feasible', 'stable'], axis=1) # drop columns you won't use like ['Feas', 'N1 ....]
y_lhc_stable = data_lhc['stable'] 
y_lhc_feasible = data_lhc['feasible'] 
y_lhc = y_lhc_stable*y_lhc_feasible
X_train_lhc, X_test_lhc, y_train_lhc, y_test_lhc = train_test_split(X_lhc, y_lhc, test_size=0.25, random_state=42)

data_imp = pd.read_csv(file_path_imp, sep = ';')
X_imp = data_imp.drop(columns = ['feasible', 'stable'], axis=1) # drop columns you won't use like ['Feas', 'N1 ....]
y_imp_stable = data_imp['stable'] 
y_imp_feasible = data_imp['feasible'] 
y_imp = y_imp_stable*y_imp_feasible
X_train_imp, X_test_imp, y_train_imp, y_test_imp = train_test_split(X_imp, y_imp, test_size=0.25, random_state=42)

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

# Train the classifier on the full training set
clf_method.fit(X_train_method, y_train_method) 
clf_method.get_params()

clf_lhc.fit(X_train_lhc, y_train_lhc) 
clf_lhc.get_params()

clf_imp.fit(X_train_imp, y_train_imp) 
clf_imp.get_params()

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

# Evaluate the classifier
accuracy_method = f1_score(y_test_method, y_pred_method)
accuracy_method_lhc = f1_score(y_test_lhc, y_pred_method_lhc)
accuracy_method_imp = f1_score(y_test_imp, y_pred_method_imp)
print("Test f1_score train  method, test method:", accuracy_method)
print("Test f1_score train  method, test lhc:", accuracy_method_lhc)
print("Test f1_score train  method, test imp:", accuracy_method_imp)

accuracy_lhc = f1_score(y_test_lhc, y_pred_lhc)
accuracy_lhc_method = f1_score(y_test_method, y_pred_lhc_method)
accuracy_lhc_imp = f1_score(y_test_imp, y_pred_lhc_imp)
print("Test f1_score train lhc, test lhc:", accuracy_lhc)
print("Test f1_score train lhc, test method:", accuracy_lhc_method)
print("Test f1_score train lhc, test imp:", accuracy_lhc_imp)

accuracy_imp = f1_score(y_test_imp, y_pred_imp)
accuracy_imp_method = f1_score(y_test_method, y_pred_imp_method)
accuracy_imp_lhc = f1_score(y_test_lhc, y_pred_imp_lhc)
print("Test f1_score train imp, test imp:", accuracy_imp)
print("Test f1_score train imp, test method:", accuracy_imp_method)
print("Test f1_score train imp, test lhc:", accuracy_imp_lhc)

# plot confusion matrix
confusion_mat_method = confusion_matrix(y_test_method, y_pred_method, labels = [0,1])
sns.heatmap(confusion_mat_method, annot=True, cmap='Paired', 
            cbar=False, fmt="d", xticklabels=[
            'Secure', 'Insecure'], yticklabels=['Secure', 'Insecure'])

confusion_mat_lhc = confusion_matrix(y_test_lhc, y_pred_lhc, labels = [0,1])
sns.heatmap(confusion_mat_lhc, annot=True, cmap='Paired', 
            cbar=False, fmt="d", xticklabels=[
            'Secure', 'Insecure'], yticklabels=['Secure', 'Insecure'])

confusion_mat_imp = confusion_matrix(y_test_imp, y_pred_imp, labels = [0,1])
sns.heatmap(confusion_mat_imp, annot=True, cmap='Paired', 
            cbar=False, fmt="d", xticklabels=[
            'Secure', 'Insecure'], yticklabels=['Secure', 'Insecure'])

# show classification report
classification_method = classification_report(y_test_method, y_pred_method) # shows precision, recall, f1-score 
print("Classification Report:\n", classification_method)

classification_lhc = classification_report(y_test_lhc, y_pred_lhc) # shows precision, recall, f1-score 
print("Classification Report:\n", classification_lhc)

classification_imp = classification_report(y_test_imp, y_pred_imp) # shows precision, recall, f1-score 
print("Classification Report:\n", classification_imp)

# Visualize the decision tree
plt.figure(figsize=(20,10))
tree.plot_tree(clf_method, filled=True, feature_names=[f'feature_{i}' for i in range(X_method.shape[1])], class_names=['insecure', 'secure'])
plt.show()

plt.figure(figsize=(20,10))
tree.plot_tree(clf_lhc, filled=True, feature_names=[f'feature_{i}' for i in range(X_lhc.shape[1])], class_names=['insecure', 'secure'])
plt.show()

plt.figure(figsize=(20,10))
tree.plot_tree(clf_imp, filled=True, feature_names=[f'feature_{i}' for i in range(X_imp.shape[1])], class_names=['insecure', 'secure'])
plt.show()

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

# Plot the top 10 features
importance_df_method.head(10).plot(kind='bar', x='Feature', y='Importance', legend=False)
plt.show()

importance_df_lhc.head(10).plot(kind='bar', x='Feature', y='Importance', legend=False)
plt.show()

importance_df_imp.head(10).plot(kind='bar', x='Feature', y='Importance', legend=False)
plt.show()


data_result = {
    'training data': ['method', 'lhc', 'importance'],
    'method test data': [accuracy_method, accuracy_lhc_method, accuracy_imp_method],
    'lhc test data': [accuracy_method_lhc, accuracy_lhc, accuracy_imp_lhc],
    'imp test data': [accuracy_method_imp, accuracy_lhc_imp, accuracy_imp],
}
    
df_result = pd.DataFrame.from_dict(data_result)
df_result.to_csv(file_path_result, sep = ';')










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

