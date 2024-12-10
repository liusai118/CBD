import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import pandas as pd
import numpy as np
from sklearn.metrics import accuracy_score, f1_score


from models import BCL_Network
import myDataSet as ms


def getFullDataLoader(X, y, batch_size):
    full_DataSet = ms.MyDataSet(input=X.reset_index(drop=True), label=y.reset_index(drop=True))
    full_DataLoader = DataLoader(dataset=full_DataSet, batch_size=batch_size, shuffle=False)
    return full_DataLoader


def predict(model_path, test_loader):

    model = BCL_Network().to(device)
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.eval()  

    all_outputs = []

    with torch.no_grad():
        for inputs, _ in test_loader:
            inputs = inputs.to(device)
            outputs = model(inputs)

            binary_outputs = (outputs > 0.5).int()
            all_outputs.append(binary_outputs.cpu().numpy())


    return np.concatenate(all_outputs).flatten()


def save_predictions(sequences, y_pred, file_path):
    predictions_df = pd.DataFrame({
        "Sequence": sequences,
        "Predicted Label": y_pred
    })
    predictions_df.to_csv(file_path, index=False)
    print(f"Predictions saved to {file_path}")

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

all_data = pd.read_csv('D:\\Programming\\python\\PycharmProjects\\ProteinDNABinding\\data\\rawdata\\NRF22\\test3', sep='\t')


print(all_data.head())
print(all_data.columns)


sequences = all_data.iloc[:, 0]
X = sequences
y = all_data.iloc[:, 1]


batch_size = 128
test_loader = getFullDataLoader(X, y, batch_size)


model_path = 'D:/Programming/python/PycharmProjects/ProteinDNABinding/NRF22/fold_3/validate_params_epoch_2.pkl'


y_pred = predict(model_path, test_loader)


output_csv_path = 'D:\\Programming\\python\\PycharmProjects\\ProteinDNABinding\\results\\predictions3.csv'
save_predictions(sequences, y_pred, output_csv_path)
