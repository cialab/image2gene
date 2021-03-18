from __future__ import print_function

import numpy as np

import argparse
import torch
import torch.utils.data as data_utils
import torch.optim as optim
from torch.autograd import Variable

from dataloader import HDF5Dataset
from slideModel import Attention

from torch.cuda import memory_cached
from torch.cuda import empty_cache

import pdb

# Training settings
parser = argparse.ArgumentParser()
parser.add_argument('--patch_size', type=int, default=32, metavar='PS',
                    help='patch size')
parser.add_argument('--gene_idx', metavar='gIdx',
                    help='index of gene')
parser.add_argument('--model',
                    help='where is the model')
parser.add_argument('--seed', type=int, default=1, metavar='S',
                    help='random seed (default: 1)')
parser.add_argument('--split', type=float, default=0.8, metavar='Sp',
                    help='random seed (default: 1)')
parser.add_argument('--file_path',
                    help='where to search for dataset .h5')
parser.add_argument('--recursive', type=bool, default=True,
                    help='to search recursively in file_path')
parser.add_argument('--load_data', type=bool, default=False,
                    help='load all data first?')
parser.add_argument('--transform', type=bool, default=True,
                    help='map values?')
parser.add_argument('--data_cache_size', type=int, default=3,
                    help='not sure')                    
parser.add_argument('--no-cuda', action='store_true', default=False,
                    help='disables CUDA training')
parser.add_argument('--multi', type=bool, default=False)

args = parser.parse_args()
args.cuda = not args.no_cuda and torch.cuda.is_available()

torch.manual_seed(args.seed)
if args.cuda:
    torch.cuda.manual_seed(args.seed)

loader_kwargs = {'num_workers': 1, 'pin_memory': True} if args.cuda else {}

if args.multi:
    genes=args.gene_idx.split(',')
    gene_idx=[]
    for i in range(0,len(genes)):
        gene_idx.append(int(genes[i]))
else:
    gene_idx=int(args.gene_idx)
  
test_loader = data_utils.DataLoader(HDF5Dataset(file_path=args.file_path,
                                                load_data=True,
                                                train=False,
                                                seed=args.seed,
                                                split=args.split,
                                                transform=args.transform,
                                                gene_idx=gene_idx,
                                                num_instances=0),
                                     batch_size=1,
                                     shuffle=False,
                                     **loader_kwargs)

def test():
    model.eval()
    with torch.no_grad():
        for batch_idx, (data, label, n) in enumerate(test_loader):
            bag_label = label[0]
            if args.cuda:
                data, bag_label = data.cuda(), bag_label.cuda()
            data, bag_label = Variable(data), Variable(bag_label)
            predicted_value = model.predictFullBag(data,10000)

            n=n.data.cpu().numpy()
            label=label.data.cpu().numpy()
            predicted_value=predicted_value.data.cpu().numpy()
            print(str(n[0])+"\t"+str(np.transpose(label[0])[0])+"\t"+str(predicted_value[0]))


if __name__ == "__main__":
    #print('# of testing bags:\t' + str(len(test_loader)))
    model=torch.load(args.model)
    model.cuda()
    test()
