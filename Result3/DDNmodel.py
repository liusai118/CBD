# === chemberta_train_model.py ===
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from transformers import RobertaModel, RobertaTokenizer
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import BRICS
import matplotlib.pyplot as plt

# === Hyperparameters ===
EMBED_DIM = 32
LR = 1e-3
EPOCHS_AUTOENCODER = 200
EPOCHS_PREDICTOR = 100

# === Load data ===
data = pd.read_csv("pearson_relation.csv")

# === Load ChemBERTa ===
device = "cuda" if torch.cuda.is_available() else "cpu"
tokenizer = RobertaTokenizer.from_pretrained("seyonec/ChemBERTa-zinc-base-v1")
chemberta = RobertaModel.from_pretrained("seyonec/ChemBERTa-zinc-base-v1").to(device)
chemberta.eval()

def smiles_to_embedding(smiles):
    with torch.no_grad():
        inputs = tokenizer(smiles, return_tensors="pt", padding=True, truncation=True, max_length=128).to(device)
        outputs = chemberta(**inputs)
        return outputs.last_hidden_state.mean(dim=1).squeeze(0)

# === Define models ===
class CellLineAutoEncoder(nn.Module):
    def __init__(self, input_dim, embed_dim):
        super().__init__()
        self.encoder = nn.Sequential(nn.Linear(input_dim, 128), nn.ReLU(), nn.Linear(128, embed_dim))
        self.decoder = nn.Sequential(nn.Linear(embed_dim, 128), nn.ReLU(), nn.Linear(128, input_dim))

    def forward(self, x):
        z = self.encoder(x)
        return self.decoder(z), z

class Predictor(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 128), nn.ReLU(),
            nn.Linear(128, 64), nn.ReLU(),
            nn.Linear(64, 1)
        )

    def forward(self, x): return self.net(x).squeeze(1)

# === Train AutoEncoder ===
pivot = data.pivot_table(index="cell_line", columns="SMILES", values="score", fill_value=0)
scaler = StandardScaler()
x_scaled = scaler.fit_transform(pivot.values)
cell_input_dim = x_scaled.shape[1]

cell_model = CellLineAutoEncoder(cell_input_dim, EMBED_DIM).to(device)
optim_ae = torch.optim.Adam(cell_model.parameters(), lr=LR)
loss_fn = nn.MSELoss()
x_tensor = torch.tensor(x_scaled, dtype=torch.float32).to(device)

for epoch in range(EPOCHS_AUTOENCODER):
    cell_model.train()
    recon, emb = cell_model(x_tensor)
    loss = loss_fn(recon, x_tensor)
    optim_ae.zero_grad(); loss.backward(); optim_ae.step()

cell_model.eval()
with torch.no_grad():
    _, cell_emb = cell_model(x_tensor)
cell_emb_df = pd.DataFrame(cell_emb.cpu().numpy(), index=pivot.index)

# === Build training dataset ===
def build_dataset(df):
    features, labels = [], []
    for _, row in tqdm(df.iterrows(), total=len(df)):
        smi = row["SMILES"]
        cell = row["cell_line"]
        score = row["score"]
        if cell not in cell_emb_df.index: continue
        smi_emb = smiles_to_embedding(smi)
        cell_vec = torch.tensor(cell_emb_df.loc[cell].values, dtype=torch.float32)
        features.append(torch.cat([smi_emb, cell_vec], dim=0))
        labels.append(score)
    return torch.stack(features), torch.tensor(labels)

X, y = build_dataset(data)
X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.1, random_state=42)

predictor_input_dim = X.shape[1]
model = Predictor(predictor_input_dim).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=LR)

# === Train predictor ===
for epoch in range(EPOCHS_PREDICTOR):
    model.train()
    optimizer.zero_grad()
    pred = model(X_train.to(device))
    loss = loss_fn(pred, y_train.to(device))
    loss.backward(); optimizer.step()

    if epoch % 10 == 0:
        model.eval()
        val_pred = model(X_val.to(device))
        val_loss = mean_squared_error(y_val.cpu().numpy(), val_pred.detach().cpu().numpy())
        print(f"Epoch {epoch}: val MSE = {val_loss:.4f}")

# === Save model and hyperparameters ===
torch.save({
    "model_state_dict": model.state_dict(),
    "cell_model_state_dict": cell_model.state_dict(),
    "cell_embedding": cell_emb_df.to_dict(),
    "cell_input_dim": cell_input_dim,
    "embed_dim": EMBED_DIM,
    "predictor_input_dim": predictor_input_dim,
}, "DDN_final_model.pt")
print("✅ Model training completed and saved to DDN_final_model.pt")

# === Scatter plot: True vs Predicted scores on validation set ===
true_vals = y_val.cpu().numpy()
pred_vals = val_pred.detach().cpu().numpy()
min_val = min(true_vals.min(), pred_vals.min())
max_val = max(true_vals.max(), pred_vals.max())

plt.figure(figsize=(6, 6))
plt.scatter(true_vals, pred_vals, alpha=0.6)
plt.plot([min_val, max_val], [min_val, max_val], 'r--', label='Ideal Fit')
plt.xlabel("True Score")
plt.ylabel("Predicted Score")
plt.title("Validation: True vs Predicted")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("validation_scatter.png", dpi=300)
print("✅ Saved scatter plot to validation_scatter.png")

# === Fragment Contribution Analysis on Training Set ===
def predict_score(smiles, cell_line, model, cell_emb_df):
    try:
        smi_emb = smiles_to_embedding(smiles)
    except:
        return None
    if cell_line not in cell_emb_df.index:
        return None
    cell_vec = torch.tensor(cell_emb_df.loc[cell_line].values, dtype=torch.float32).to(device)
    input_vec = torch.cat([smi_emb, cell_vec])
    with torch.no_grad():
        score = model(input_vec.unsqueeze(0)).item()
    return score

fragment_records = []
for _, row in tqdm(data.iterrows(), total=len(data)):
    smiles = row["SMILES"]
    cell_line = row["cell_line"]
    original_score = predict_score(smiles, cell_line, model, cell_emb_df)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or original_score is None:
        continue
    try:
        fragments = BRICS.BRICSDecompose(mol)
    except:
        continue
    for frag in fragments:
        modified_smiles = ".".join([f for f in fragments if f != frag])
        modified_score = predict_score(modified_smiles, cell_line, model, cell_emb_df)
        if modified_score is None:
            continue
        contribution = original_score - modified_score
        fragment_records.append({
            "compound": smiles,
            "cell_line": cell_line,
            "fragment": frag,
            "original_score": original_score,
            "modified_score": modified_score,
            "contribution": contribution
        })

frag_df = pd.DataFrame(fragment_records)
frag_df.to_csv("fragment_contributions.csv", index=False)
print("✅ Fragment contribution analysis completed and saved to fragment_contributions.csv")

