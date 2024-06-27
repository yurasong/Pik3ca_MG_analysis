############################################################
# Please activate conda environment: conda activate pyscenic
############################################################

import pandas as pd

# Import data using pandas
file_names = [f"auc_mtx{i}.csv" for i in range(9)]
data_frames = [pd.read_csv(file) for file in file_names]

# Extract column names from each data frame and merge them
all_columns = set()
for df in data_frames:
    all_columns.update(df.columns)

# Remove the 'Cell' column since it is used as a key for further manipulation
all_columns.discard("Cell")
columns_list = list(all_columns)

# Initialize result dataframe with 'Cell' column
result = pd.DataFrame(data_frames[0]["Cell"])

# Run the loop to calculate the average for each column
for column in columns_list:
    count = 0
    result[column] = 0
    
    # Add values from each data frame if the column exists
    for df in data_frames:
        if column in df.columns:
            count += 1
            result[column] += df[column].astype("float")
    
    # Calculate the average
    if count > 0:
        result[column] /= count

# Export the result to a CSV file
result.to_csv("AUC_Average.csv", index=False)