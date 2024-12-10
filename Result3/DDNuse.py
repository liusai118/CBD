import pandas as pd
import torch
import torch.nn as nn
from rdkit import Chem
from rdkit.Chem import BRICS
from transformers import RobertaModel, RobertaTokenizer
TEST_FILE = "newtest.csv"
MODEL_FILE = "chemberta_final_model.pt"
TARGET_CELL = "HUVEC"
device = "cuda" if torch.cuda.is_available() else "cpu"
tokenizer = RobertaTokenizer.from_pretrained("seyonec/ChemBERTa-zinc-base-v1")
chemberta = RobertaModel.from_pretrained("seyonec/ChemBERTa-zinc-base-v1").to(device)
chemberta.eval()
def smiles_to_embedding(smiles):
    with torch.no_grad():
        inputs = tokenizer(smiles, return_tensors="pt", padding=True, truncation=True, max_length=128).to(device)
        outputs = chemberta(**inputs)
        return outputs.last_hidden_state.mean(dim=1).squeeze(0)
class CellLineAutoEncoder(nn.Module):
    def __init__(self, input_dim, embed_dim):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 128), nn.ReLU(),
            nn.Linear(128, embed_dim)
        )
        self.decoder = nn.Sequential(
            nn.Linear(embed_dim, 128), nn.ReLU(),
            nn.Linear(128, input_dim)
        )
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
    def forward(self, x):
        return self.net(x).squeeze(1)

checkpoint = torch.load(MODEL_FILE)
cell_emb_df = pd.DataFrame.from_dict(checkpoint["cell_embedding"])
input_dim = checkpoint["cell_input_dim"]
embed_dim = checkpoint["embed_dim"]
predictor_input_dim = checkpoint["predictor_input_dim"]

model = Predictor(predictor_input_dim).to(device)
model.load_state_dict(checkpoint["model_state_dict"])
model.eval()

cell_model = CellLineAutoEncoder(input_dim=input_dim, embed_dim=embed_dim).to(device)
cell_model.load_state_dict(checkpoint["cell_model_state_dict"])
cell_model.eval()

test_data = pd.read_csv(TEST_FILE)
if "smiles" not in test_data.columns:
    test_data.columns = [c.lower() for c in test_data.columns]

results = []
for _, row in test_data.iterrows():
    smi = row["smiles"]
    cell = row["cell_line"]
    if cell not in cell_emb_df.index:
        results.append(None)
        continue
    smi_emb = smiles_to_embedding(smi)
    cell_vec = torch.tensor(cell_emb_df.loc[cell].values, dtype=torch.float32).to(device)
    input_vec = torch.cat([smi_emb, cell_vec])
    with torch.no_grad():
        score = model(input_vec.unsqueeze(0)).item()
    results.append(score)

test_data["predicted_score"] = results
test_data = test_data.sort_values("predicted_score")
test_data.to_csv("newtest_scored.csv", index=False)
