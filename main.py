import sys
import pandas as pd
from joblib import load
from sklearn.metrics import matthews_corrcoef, f1_score, roc_auc_score, balanced_accuracy_score, precision_score, recall_score, confusion_matrix
import os

if __name__ == "__main__":
    taxa = [
        "bpertussis", "corynebacterium", "orthopoxvirus", "ecoli", "enterobacteriaceae", "lentivirus", "mtuberculosis",
        "paeruginosa", "smansoni", "tgondii", "pfalciparum", "ctrachomatis", "human_gammaherpesvirus_4", "influenza_a",
        "cdifficile", "filoviridae", "measles_morbilivirus", "ovolvulus", "mononegavirales"
    ]
    taxa_string = [
        "B. pertussis", "Corynebacterium", "Orthopoxvirus", "E. coli", "Enterobacteriaceae", "Lentivirus",
        "M. tuberculosis", "P. aeruginosa", "S. mansoni", "T. gondii", "P. falciparum", "C. trachomatis",
        "Human Gammaherpesvirus 4", "Influenza A", "C. difficile", "Filoviridae", "Measles morbilivirus",
        "Ovolvulus", "Mononegavirales"
    ]
    thresholds_esm2 = [0.422, 0.5, 0.5, 0.492, 0.432, 0.6, 0.507, 0.445, 0.262, 0.443, 0.571, 0.516, 0.435, 0.559, 0.241, 0.365, 0.5, 0.501, 0.325]
    thresholds_esm1b = [0.801, 0.468, 0.497, 0.481, 0.201, 0.841, 0.226, 0.7, 0.169, 0.274, 0.69, 0.879, 0.555, 0.358, 0.53, 0.217, 0.313, 0.484, 0.194]

    print("Working directory:", os.getcwd())

    # Argumentos: modelo (esm1b ou esm2) e lista de taxa
    if len(sys.argv) < 2:
        print("Usage: python script.py <model_type: esm1b|esm2> [taxa1 taxa2 ...|all]")
        sys.exit(1)

    model_type = sys.argv[1].lower()  # 'esm1b' ou 'esm2'
    if model_type not in ['esm1b', 'esm2']:
        print("First argument must be 'esm1b' or 'esm2'")
        sys.exit(1)

    taxa_inputs = sys.argv[2:]

    # Para esm1b, "all" NÃO é permitido:
    if model_type == 'esm1b':
        if "all" in taxa_inputs:
            print("The option 'all' is not allowed for esm1b. Please specify the taxa explicitly.")
            sys.exit(1)
        thresholds = thresholds_esm1b
        model_prefix = 'epitopetransfer_esm1b'
        input_prefix = 'epitopetransfer_esm1b'
        taxa_to_process = [t for t in taxa_inputs if t in taxa]
    else:
        thresholds = thresholds_esm2
        model_prefix = 'epitopetransfer_esm2'
        input_prefix = 'epitopetransfer_esm2'
        if "all" in taxa_inputs:
            taxa_to_process = taxa
        else:
            taxa_to_process = [t for t in taxa_inputs if t in taxa]

    for taxa_input in taxa_to_process:
        index = taxa.index(taxa_input)
        string_taxa = taxa_string[index]
        best_threshold = thresholds[index]

        model_path = f"./models/{model_prefix}_{taxa_input}.pkl"
        test_features_csv_path = f"./input/{input_prefix}/{taxa_input}/processed_test_features.csv"
        test_labels_csv_path = f"./input/{input_prefix}/{taxa_input}/test_labels.csv"

        if not (os.path.exists(model_path) and os.path.exists(test_features_csv_path) and os.path.exists(test_labels_csv_path)):
            print(f"Files not found for {taxa_input}:")
            print(f"  Model: {model_path}")
            print(f"  Features: {test_features_csv_path}")
            print(f"  Labels: {test_labels_csv_path}")
            continue

        model = load(model_path)
        test_features = pd.read_csv(test_features_csv_path)
        test_labels = pd.read_csv(test_labels_csv_path).iloc[:, 0]

        test_probs = model.predict_proba(test_features)[:, 1]
        test_predictions = (test_probs > best_threshold).astype(int)

        # Métricas
        test_mcc = matthews_corrcoef(test_labels, test_predictions)
        f1 = f1_score(test_labels, test_predictions)
        auc_score = roc_auc_score(test_labels, test_probs)
        bacc = balanced_accuracy_score(test_labels, test_predictions)
        ppv = precision_score(test_labels, test_predictions, zero_division=0)
        sens = recall_score(test_labels, test_predictions, zero_division=0)

        tn, fp, fn, tp = confusion_matrix(test_labels, test_predictions).ravel()
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        npv = tn / (tn + fn) if (tn + fn) > 0 else 0

        title = f"{string_taxa}"
        results = (
            f"AUC: {auc_score:.3f}\n"
            f"F1 : {f1:.3f}\n"
            f"MCC: {test_mcc:.3f}\n"
            f"BACC: {bacc:.3f}\n"
            f"PPV: {ppv:.3f}\n"
            f"NPV: {npv:.3f}\n"
            f"SENS: {sens:.3f}\n"
            f"SPEC: {specificity:.3f}"
        )
        print(f"{title}\n{'-' * len(title)}\n{results}\n")
