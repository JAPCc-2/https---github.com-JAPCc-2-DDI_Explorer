import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, roc_auc_score

# Load features
X = np.load('X_features.npy')

# Dummy labels for now (temporary)
# Let's assume: first half = no interaction (0), second half = interaction (1)
y = np.array([0] * (len(X)//2) + [1] * (len(X) - len(X)//2))

# Split into train and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Initialize Random Forest
model = RandomForestClassifier(n_estimators=100, random_state=42)

# Train the model
model.fit(X_train, y_train)

# Predict
y_pred = model.predict(X_test)

# Evaluate
print("✅ Model Evaluation Results:")
print(classification_report(y_test, y_pred))

# ROC AUC Score
y_pred_proba = model.predict_proba(X_test)[:,1]
roc_auc = roc_auc_score(y_test, y_pred_proba)
print(f"ROC-AUC Score: {roc_auc:.3f}")

# Save model (optional)
import pickle
with open('random_forest_ddi_model.pkl', 'wb') as f:
    pickle.dump(model, f)

print("✅ Random Forest Model saved as 'random_forest_ddi_model.pkl'")
