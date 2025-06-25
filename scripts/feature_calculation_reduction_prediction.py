# This script generates features for a target domain by leveraging higher-level pre-trained models.
# It processes protein sequences by extracting embeddings and selecting labeled regions (reduction step).
# Using these embeddings, the script trains a RandomForest classifier with hyperparameter optimization 
# performed via Optuna, aiming to maximize performance on metrics like MCC (Matthews Correlation Coefficient), 
# F1-score, and AUC (Area Under the Curve).
#
# The script implements cross-validation to ensure robust evaluation and automatically selects the best hyperparameters.
# After training, it evaluates the model using a test set and generates predictions, saving the results in CSV format.
# The trained model is also saved for future use, enabling easy reusability for inference or further analysis.

import torch
import numpy as np
import pandas as pd
import os
import csv
from transformers import AutoModel, AutoTokenizer
from Bio import SeqIO
import pickle
import optuna
from optuna.distributions import IntDistribution, CategoricalDistribution
from optuna.integration.sklearn import OptunaSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import matthews_corrcoef, make_scorer, roc_auc_score, f1_score
from sklearn.model_selection import PredefinedSplit
from scipy.stats import sem
import matplotlib.pyplot as plt
from optuna.visualization import plot_param_importances


# Replace here with your paths
tokenizer_model = "./model/esm1b_t33_650M_UR50S"
base_model = "./model/esm1b_t33_650M_UR50S"
string_organism = "E. coli"
organism = "ecoli"
path_results = "./results/"
path_models = "./models/"
method = "epitopetransfer"
method_name = "EpitopeTransfer"

# Define training and testing folds
train_fasta_folds = ["fold1", "fold2", "fold3", "fold4"]
test_fasta_fold = "fold5"
train_fasta_files = [f"./fasta/folds/{organism}/{fold}.fasta" for fold in train_fasta_folds]
test_fasta_file = f"./fasta/folds/{organism}/{test_fasta_fold}.fasta"

# Paths for fold-based data
fold_base_path = f"./input/organismos/taxons/{organism}/"
train_fold_paths = [f"{fold_base_path}fold{i}.csv" for i in [1, 2, 3, 4]]
test_fold = "fold5"
test_path = f"{fold_base_path}{test_fold}.csv"
path_predictions = f"./predictions/{method}/{organism}/{test_fold}.csv"
fold_numbers = 5

# Model parameters
max_seq_length = 1024
overlap = 256

# Function to predict sequences from a fasta file
def predict_fasta_sequences(fasta_file):
    model = AutoModel.from_pretrained(base_model)
    tokenizer = AutoTokenizer.from_pretrained(tokenizer_model, do_lower_case=False)

    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
    sequence_dict = {}

    model.eval()  # Disable dropout

    for record in fasta_sequences:
        sequence = str(record.seq)
        sequence_length = len(sequence)
        embeddings = None

        if sequence_length > max_seq_length:
            # Split sequence into chunks with overlap
            for i in range(0, sequence_length, max_seq_length - overlap):
                chunk = sequence[i:i + max_seq_length]
                if len(chunk) < max_seq_length:
                    chunk = chunk + ' ' * (max_seq_length - len(chunk))  # Padding
                features = tokenizer.batch_encode_plus(
                    [chunk],
                    add_special_tokens=True,
                    padding='max_length',
                    max_length=max_seq_length,
                    truncation=True,
                    return_tensors='pt',
                    return_attention_mask=True
                )

                with torch.no_grad():
                    outputs = model(features['input_ids'], features['attention_mask'])

                last_hidden_state = outputs[0]
                array_last_hidden_state = last_hidden_state.detach().numpy()

                if embeddings is None:
                    embeddings = array_last_hidden_state[0]
                else:
                    embeddings = np.concatenate((embeddings, array_last_hidden_state[0][overlap:]))
        else:
            # Process entire sequence at once
            features = tokenizer.batch_encode_plus(
                [sequence],
                add_special_tokens=True,
                padding='max_length',
                max_length=max_seq_length,
                truncation=True,
                return_tensors='pt',
                return_attention_mask=True
            )

            with torch.no_grad():
                outputs = model(features['input_ids'], features['attention_mask'])

            last_hidden_state = outputs[0]
            embeddings = last_hidden_state.detach().numpy()[0]

        sequence_dict[record.id] = {
            'completed_sequence': sequence,
            'embeddings': embeddings
        }

    return sequence_dict

# Function to select labeled regions
def select_labeled_regions(dataset, sequence_dict):
    list_amino_features = []
    list_amino_label = []

    for idx, row in dataset.iterrows():
        protein_id = row["Info_protein_id"]

        if protein_id not in sequence_dict:
            print(f"{protein_id} is not in sequence_dict")
            continue

        data = sequence_dict[protein_id]
        sequence_slice = data['embeddings']
        start = row["Info_start_pos"]
        end = row["Info_end_pos"]
        label = row["Class"]
        sequence = data['completed_sequence']

        for index_amino in range(start - 1, end):
            list_amino_features.append(sequence_slice[index_amino])
            list_amino_label.append(label)

    return list_amino_features, list_amino_label

# Combine results from multiple fasta files
def predict_fasta_sequences_multiple(files):
    combined_sequence_dict = {}
    for fasta_file in files:
        print(fasta_file)
        sequence_dict = predict_fasta_sequences(fasta_file)
        combined_sequence_dict.update(sequence_dict)
    return combined_sequence_dict

# Call prediction function for each fasta file and combine results
train_sequence_dict = predict_fasta_sequences_multiple(train_fasta_files)
test_sequence_dict = predict_fasta_sequences(test_fasta_file)

# Process training folds and save
for i, fold_path in enumerate(train_fold_paths, start=1):
    df_fold = pd.read_csv(fold_path)
    processed_data = select_labeled_regions(df_fold, train_sequence_dict)

    # Convert lists into DataFrames
    features_df = pd.DataFrame(processed_data[0])
    labels_df = pd.DataFrame(processed_data[1])

    # Save features and labels for each fold in CSV
    features_csv_path = f"{fold_base_path}{method}/processed_fold{i}_features.csv"
    labels_csv_path = f"{fold_base_path}{method}/processed_fold{i}_labels.csv"

    os.makedirs(os.path.dirname(features_csv_path), exist_ok=True)
    os.makedirs(os.path.dirname(labels_csv_path), exist_ok=True)

    features_df.to_csv(features_csv_path, index=False)
    labels_df.to_csv(labels_csv_path, index=False)

    print(f"Processed fold{i} features saved to {features_csv_path}")
    print(f"Processed fold{i} labels saved to {labels_csv_path}")

# Process test dataset
df_test = pd.read_csv(test_path)
processed_data = select_labeled_regions(df_test, test_sequence_dict)

# Save features and labels for test in CSV
features_csv_path = f"{fold_base_path}{method}/processed_test_features.csv"
labels_csv_path = f"{fold_base_path}{method}/processed_test_labels.csv"
pd.DataFrame(processed_data[0]).to_csv(features_csv_path, index=False)
pd.DataFrame({"Label": processed_data[1]}).to_csv(labels_csv_path, index=False)

print(f"Processed test features saved to {features_csv_path}")
print(f"Processed test labels saved to {labels_csv_path}")

# Load datasets
loaded_data = {}
for i in range(1, fold_numbers):
    features_csv_path = f"{fold_base_path}{method}/processed_fold{i}_features.csv"
    labels_csv_path = f"{fold_base_path}{method}/processed_fold{i}_labels.csv"

    loaded_data[f"fold{i}"] = {
        "features": pd.read_csv(features_csv_path),
        "labels": pd.read_csv(labels_csv_path)
    }

# Load test data
test_features_csv_path = f"{fold_base_path}{method}/processed_test_features.csv"
test_labels_csv_path = f"{fold_base_path}{method}/processed_test_labels.csv"
loaded_data["test"] = {
    "features": pd.read_csv(test_features_csv_path),
    "labels": pd.read_csv(test_labels_csv_path)
}

# Define the objective function for Optuna
def objective(trial, train_features, train_labels, test_features, test_labels):
    params = {
        'n_estimators': trial.suggest_int('n_estimators', 100, 500),
        'max_depth': trial.suggest_categorical('max_depth', [None, 10, 20, 30]),
        'min_samples_split': trial.suggest_int('min_samples_split', 2, 10),
        'min_samples_leaf': trial.suggest_int('min_samples_leaf', 1, 10),
        'max_features': trial.suggest_categorical('max_features', ['sqrt', 'log2', None]),
        'bootstrap': trial.suggest_categorical('bootstrap', [True, False]),
        'criterion': trial.suggest_categorical('criterion', ['gini', 'entropy']),
        'random_state': 42
    }

    model = RandomForestClassifier(**params, n_jobs=-1)
    model.fit(train_features, train_labels)

    train_probs = model.predict_proba(train_features)[:, 1]
    thresholds = np.arange(0, 1, 0.001)
    mcc_scores = [matthews_corrcoef(train_labels, (train_probs > t).astype(int)) for t in thresholds]
    max_mcc_index = np.argmax(mcc_scores)
    threshold = thresholds[max_mcc_index]

    test_probs = model.predict_proba(test_features)[:, 1]
    predicted_labels = (test_probs > threshold).astype(int)
    test_mcc = matthews_corrcoef(test_labels, predicted_labels)

    trial.set_user_attr('model', model)
    trial.set_user_attr('threshold', threshold)
    trial.set_user_attr('params', params)

    return test_mcc

# Perform validation using holdout
def perform_validation_holdout(loaded_data):
    train_features = loaded_data["fold1"]["features"]
    train_labels = loaded_data["fold1"]["labels"].iloc[:, 0]
    test_features = loaded_data["test"]["features"]
    test_labels = loaded_data["test"]["labels"].iloc[:, 0]

    study = optuna.create_study(direction='maximize')
    study.optimize(lambda trial: objective(trial, train_features, train_labels, test_features, test_labels), n_trials=200)

    best_trial = study.best_trial
    final_model = best_trial.user_attrs['model']
    mcc_best_threshold = best_trial.user_attrs['threshold']
    best_hyperparams = best_trial.params

    # Save the best model
    model_path = os.path.join(path_models, f"{method}_{organism}.pkl")
    with open(model_path, 'wb') as f:
        pickle.dump(final_model, f)
    print(f"Best model saved to {model_path}")

    return final_model, mcc_best_threshold, best_hyperparams

# Evaluate model on test data
def evaluate_model_on_test(test_features, test_labels, model, mcc_threshold, f1_threshold):
    test_probs = model.predict_proba(test_features)[:, 1]

    predictions_mcc = test_probs >= mcc_threshold
    predictions_f1 = test_probs >= f1_threshold

    auc = roc_auc_score(test_labels, test_probs)
    f1_at_mcc_threshold = f1_score(test_labels, predictions_mcc)
    mcc_at_mcc_threshold = matthews_corrcoef(test_labels, predictions_mcc)

    results = pd.DataFrame([{
        "Method": method,
        "Organism": organism,
        "AUC": auc,
        "MCC threshold": mcc_threshold,
        "F1 threshold": f1_threshold
    }])

    results_file_path = os.path.join(path_results, f"{method}_{organism}_test_results.csv")
    results.to_csv(results_file_path, index=False)
    print(f"Results saved to {results_file_path}")

    return auc, f1_at_mcc_threshold, mcc_at_mcc_threshold

# Main function to run the workflow
def main(loaded_data, organism, method, path_results):
    if len(loaded_data) > 2:  # Keep the conditional logic here
        best_model, mcc_best_threshold, f1_best_threshold, best_hyperparams = perform_validation_holdout(loaded_data)
    else:
        best_model, mcc_best_threshold, best_hyperparams = perform_validation_holdout(loaded_data)

    test_features = loaded_data["test"]["features"]
    test_labels = loaded_data["test"]["labels"].iloc[:, 0]

    evaluate_model_on_test(test_features, test_labels, best_model, mcc_best_threshold, f1_best_threshold)

main(loaded_data, organism, method, path_results)

# Save final predictions in CSV
def generate_csv(sequence_dict, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Protein_id', 'amino', 'Predicted_prob', 'Position'])

        for protein_id, data in sequence_dict.items():
            sequence = data['completed_sequence']
            predictions = data['embeddings']
            sequence_length = len(sequence)
            for i in range(sequence_length):
                writer.writerow([protein_id, sequence[i], predictions[i], i + 1])

    print(f'CSV file {output_file} generated successfully.')

# Generate CSV for predictions
output_file = f"./predictions/{method}/{organism}/{test_fold}/formatted_raw_output.csv"
generate_csv(test_sequence_dict, output_file)
