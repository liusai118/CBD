from torch.utils.data import Dataset
import torch
import load_data as ld
import numpy as np


class MyDataSet(Dataset):
    def __init__(self, input, label):
        self.input_seq = input
        self.output = label

    def __getitem__(self, index):
        input_seq_origin = self.input_seq[index]  
        input_seq = np.array(ld.k_mer_stride(input_seq_origin, 3, 1)).T
        input_seq = np.array(ld.k_mer_stride(input_seq_origin, 3, 1)).T
        input_seq = torch.from_numpy(input_seq).type(torch.FloatTensor).cuda()
        output_seq = self.output[index]
        output_seq = torch.Tensor([output_seq]).cuda()
        return input_seq, output_seq  

    def __len__(self):
        return len(self.input_seq) 

    def get_original_sequence(self, index):
        """
        返回原始输入序列
        :param index: 序列的索引
        :return: 原始输入序列
        """
        return self.input_seq[index]


