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

# -------------------------------------------------------------------
# Process data
# Convert categorical text data into a numerical format so that machine learning algorithms can process them
# Create separate label encoders for each column
label_encoder_ontology     = LabelEncoder()
label_encoder_significance = LabelEncoder()

# Fit and transform the data
# cancerGenes_olapFinal['generalOntology'] = label_encoder_ontology.fit_transform(cancerGenes_olapFinal['generalOntology'])
cancerGenes_olapFinal['generalOntology_encoded'] = label_encoder_ontology.fit_transform(cancerGenes_olapFinal['generalOntology'])
cancerGenes_olapFinal['reclassifiedSignificance_encoded'] = label_encoder_significance.fit_transform(cancerGenes_olapFinal['reclassifiedSignificance'])

# Retrieve the mapping of original category names to encoded numbers
ontology_categories = dict(zip(label_encoder_ontology.classes_, range(len(label_encoder_ontology.classes_))))
significance_categories = dict(zip(label_encoder_significance.classes_, range(len(label_encoder_significance.classes_))))

# Print the mapping
print("General Ontology Categories Mapping:", ontology_categories)
print("Reclassified Significance Categories:", significance_categories)

# Split data into training and testing sets (70% training, 30% testing)
x = cancerGenes_olapFinal[['generalOntology_encoded']]  # Features
y = cancerGenes_olapFinal['reclassifiedSignificance_encoded']   # Target variable

# Ensure that the same rows are always used for training and testing with random_state
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.3, random_state=42)

# Build the model
# 100 decision trees created by random forest model
clf = RandomForestClassifier(n_estimators=100)
clf.fit(x_train, y_train)

# Evaluate the model
y_pred = clf.predict(x_test)
print(classification_report(y_test, y_pred))
print("Accuracy:", accuracy_score(y_test, y_pred))

# -------------------------------------------------------------------
# Cross validation (5-fold cross-validation)
scores = cross_val_score(clf, x, y, cv=5)
print("Cross-validated scores:", scores)

# Provides the average score across all folds
print("Average score:", scores.mean())

# Save model
joblib.dump(clf, 'random_forest_model.pkl')

# -------------------------------------------------------------------
# Ensure that label_encoder_significance is fitted with all categories
label_encoder_significance.fit(cancerGenes_olapFinal['reclassifiedSignificance'])

# Check the fitted labels
print("Fitted Labels:", label_encoder_significance.classes_)

# Predict on the entire dataset
cancerGenes_olapFinal['predictedSignificance'] = clf.predict(cancerGenes_olapFinal[['generalOntology_encoded']])

# Add names of the categories for easier interpretation using inverse transformation
cancerGenes_olapFinal['predictedSignificance_category'] = label_encoder_significance.inverse_transform(cancerGenes_olapFinal['predictedSignificance'])

# Analyze specific cases
# For example, filter for synonymous mutations
synonymous_mutations = cancerGenes_olapFinal[cancerGenes_olapFinal['generalOntology'] == 'Synonymous Variant']
print(synonymous_mutations[['generalOntology', 'reclassifiedSignificance', 'predictedSignificance_category']])

# Statistical Summary
# Group by 'generalOntology' and see the distribution of predicted categories
summary = cancerGenes_olapFinal.groupby('generalOntology')['predictedSignificance_category'].value_counts().unstack()
print(summary)