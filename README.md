# Cancer-Detection

Detailed Explanation and Steps of the project can be viewed at: https://youtu.be/UL7GO2DI5TQ

To run the project model, following are the steps:

Step 1: Generating the csv files for different cancer types

Run the command: python matrix_generation.py
Input to this code requires the miRNA biomarker data and file IDs to each patient provided as the data  folder and file_id.tsv respectively. 
Output of this code is the Breast_matrix.csv, Lung_matrix.csv, Kidney_matrix.csv, Stomach_matrix.csv, Head and Neck_matrix.csv, Prostate_matrix.csv, Liver_matrix.csv, Thyroid_matrix.csv files.

Step 2: Obtaining important features and classification accuracy for each cancer type

Run the command: python predict.py
Input to this code are the matrix.csv files generated in step 1.
Output of this code is the cancertype_accuracy.txt (data ensemble classification accuracy), cancertype_features.txt (important features), cancertype_graph.png (feature importance graph) files where cancertype is Breast, Lung, Kidney, Stomach, Head and Neck, Prostate, Liver and Thyroid cancer

Project Members: Tamara Anne Fernandes, Khushboo Korani & Janhavi Karekar
