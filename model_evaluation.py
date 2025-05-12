import pickle
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
import numpy as np

# Load the trained Random Forest model
model = pickle.load(open('random_forest_ddi_model.pkl', 'rb'))

# Example test data (to be replaced this with actual test data)
y_true = [1, 0, 1, 1, 0, 1, 0, 1]  # True labels (1 = Interaction, 0 = No Interaction)
y_pred = [1, 0, 1, 0, 0, 1, 0, 1]  # Predicted labels by model

# Confusion Matrix
cm = confusion_matrix(y_true, y_pred)
print("Confusion Matrix:\n", cm)

# Calculate Accuracy
accuracy = accuracy_score(y_true, y_pred)
print(f"Accuracy: {accuracy:.2f}")

# Calculate Precision
precision = precision_score(y_true, y_pred)
print(f"Precision: {precision:.2f}")

# Calculate Recall
recall = recall_score(y_true, y_pred)
print(f"Recall: {recall:.2f}")

# Calculate F1 Score
f1 = f1_score(y_true, y_pred)
print(f"F1 Score: {f1:.2f}")

# Calculate ROC-AUC
auc = roc_auc_score(y_true, y_pred)
print(f"ROC-AUC: {auc:.2f}")
