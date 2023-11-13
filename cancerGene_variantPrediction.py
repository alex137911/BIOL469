import os
import pandas as pd
import joblib
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score
from sklearn.model_selection import cross_val_score

# -------------------------------------------------------------------
# Set working directory
os.chdir("C:/Users/acale/OneDrive/Documents/Waterloo BME/4A/BIOL 469/Final Project/BIOL469/Data/")

# Read in the data
cancerGenes_olapFinal = pd.read_csv('cancerGenes_olapFinal.tsv', sep='\t')

# Process data
# Converts categorical text data into a numerical format so that machine learning algorithms can process them
# For example, if you have categories like ["Benign", "Malignant"], they might be converted to [0, 1].
label_encoder = LabelEncoder()
cancerGenes_olapFinal['generalOntology']          = label_encoder.fit_transform(cancerGenes_olapFinal['generalOntology'])
cancerGenes_olapFinal['reclassifiedSignificance'] = label_encoder.fit_transform(cancerGenes_olapFinal['reclassifiedSignificance'])

# Split data into training and testing sets (70% training, 30% testing)
x = cancerGenes_olapFinal[['generalOntology']]         # Features
y = cancerGenes_olapFinal['reclassifiedSignificance']  # Target variable

x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.3, random_state=42)

# Build the model
# 100 decision trees created by random forest model
clf = RandomForestClassifier(n_estimators=100)
clf.fit(x_train, y_train)

# Evaluate the model
y_pred = clf.predict(x_test)
print(classification_report(y_test, y_pred))
print("Accuracy:", accuracy_score(y_test, y_pred))

# Cross validation (5-fold cross-validation)
scores = cross_val_score(clf, x, y, cv=5)
print("Cross-validated scores:", scores)

# Provides the average score across all folds
print("Average score:", scores.mean())

# Save model
joblib.dump(clf, 'random_forest_model.pkl')