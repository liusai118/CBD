import pytest
from torch.utils.data import DataLoader
import myDataSet as ms

@pytest.fixture
def myDataLoader():

    input_data = ['ACGT', 'CGTA', 'GTAC'] 
    label_data = [0, 1, 0]  
    dataset = ms.MyDataSet(input=input_data, label=label_data)
    loader = DataLoader(dataset, batch_size=64, shuffle=True)
    return loader
