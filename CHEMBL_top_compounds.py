import pandas as pd

# Load the data
targets_df = pd.read_csv('chembl_targets.csv')
compounds_df = pd.read_csv('chembl_compounds.csv')

# Inspect the data
print("Targets DataFrame:")
print(targets_df.head())

print("\nCompounds DataFrame:")
print(compounds_df.head())

# Filter compounds based on activity type and value (Example: IC50, EC50)
filtered_compounds = compounds_df[
    (compounds_df['activity_type'].str.contains('IC50|EC50', na=False)) &
    (compounds_df['activity_value'].astype(float) < 1000)  # Example threshold for activity value
]

# Rank compounds by activity value (lower values indicate higher potency)
ranked_compounds = filtered_compounds.sort_values(by='activity_value')

# Display top 10 compounds
print("\nTop 10 Compounds:")
print(ranked_compounds.head(10))

# Save the filtered and ranked compounds to a CSV file
ranked_compounds.to_csv('ranked_compounds.csv', index=False)