import numpy as np
import pandas as pd
import argparse
from sklearn.metrics import accuracy_score, f1_score


parser = argparse.ArgumentParser()
parser.add_argument('-true_label_path', type=str)
parser.add_argument('-prediction_label_path', type=str)


def main():
    args = parser.parse_args()
    test_y = pd.read_csv(args.true_label_path)
    pred_y = pd.read_csv(args.prediction_label_path)
    acc_ctp = accuracy_score(np.array(test_y['celltype']), np.array(pred_y['celltype']))
    acc_subctp = accuracy_score(np.array(test_y['subcelltype']), np.array(pred_y['subcelltype']))
    f1_ctp = f1_score(np.array(test_y['celltype']), np.array(pred_y['celltype']), average='macro')
    f1_subctp = f1_score(np.array(test_y['subcelltype']), np.array(pred_y['subcelltype']), average='macro')

    print("Accuracy for big cell type:")
    print(acc_ctp)
    print("Accuracy for sub-cell type:")
    print(acc_subctp)
    print("f1 for big cell type:")
    print(f1_ctp)
    print("f1 for sub-cell type:")
    print(f1_subctp)


if __name__ == "__main__":
    main()
