from __future__ import print_function

import numpy as np

import argparse
import torch
import torch.utils.data as data_utils
import torch.optim as optim
from torch.autograd import Variable

from dataloader import HDF5Dataset
from slideModel import Attention
from slideModel import MultiRegressionAttention
from slideModel import MultiRegressionAttention2

from torch.cuda import memory_cached
from torch.cuda import empty_cache

# Training settings
parser = argparse.ArgumentParser()
parser.add_argument('--epochs', type=int, default=20, metavar='N',
                    help='number of epochs to train (default: 20)')
parser.add_argument('--patch_size', type=int, metavar='PS',
                    help='patch size')
parser.add_argument('--gene_idx', metavar='gIdx',
                    help='index of gene')
parser.add_argument('--lr', type=float, default=0.0003, metavar='LR',
                    help='learning rate (default: 0.01)')
parser.add_argument('--reg', type=float, default=0.0005, metavar='R',
                    help='weight decay')
parser.add_argument('--num_bags_train', type=int, default=200, metavar='NTrain',
                    help='number of bags in training set')
parser.add_argument('--num_bags_test', type=int, default=50, metavar='NTest',
                    help='number of bags in test set')
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
parser.add_argument('--train_only', default=False)
parser.add_argument('--cutoff', type=float, default=0.0)
parser.add_argument('--multi', type=bool, default=False)
parser.add_argument('--patience', type=int, default=1000)
parser.add_argument('--model_prefix', type=str, default='hi')

args = parser.parse_args()
args.cuda = not args.no_cuda and torch.cuda.is_available()

torch.manual_seed(args.seed)
if args.cuda:
    torch.cuda.manual_seed(args.seed)

print('Load Train and Test Set')
loader_kwargs = {'num_workers': 1, 'pin_memory': True} if args.cuda else {}

if args.multi:
    genes=args.gene_idx.split(',')
    gene_idx=[]
    for i in range(0,len(genes)):
        gene_idx.append(int(genes[i]))
else:
    gene_idx=int(args.gene_idx)

train_loader = data_utils.DataLoader(HDF5Dataset(file_path=args.file_path,
                                                 load_data=True,
                                                 train=True,
                                                 seed=args.seed,
                                                 split=args.split,
                                                 transform=args.transform,
                                                 gene_idx=gene_idx,
                                                 num_instances=5000),
                                     batch_size=1,
                                     shuffle=False,
                                     **loader_kwargs)



if not args.train_only:
    test_loader = data_utils.DataLoader(HDF5Dataset(file_path=args.file_path,
                                                    load_data=True,
                                                    train=False,
                                                    seed=args.seed,
                                                    split=args.split,
                                                    transform=args.transform,
                                                    gene_idx=gene_idx,
                                                    num_instances=5000),
                                         batch_size=1,
                                         shuffle=False,
                                         **loader_kwargs)

print('Init Model')
if args.multi:
    model = MultiRegressionAttention2(ps=args.patch_size,nout=len(gene_idx))
else:
    model = Attention(ps=args.patch_size)
    
if args.cuda:
    model.cuda()

optimizer = optim.Adam(model.parameters(), lr=args.lr, betas=(0.9, 0.999), weight_decay=args.reg, amsgrad=True)

patience=0.0
c_loss=1000.0
def train(epoch):
    model.train()
    train_loss = 0.
    train_error = 0.
    for batch_idx, (data, label, n) in enumerate(train_loader):
        bag_label = label[0]
        if args.cuda:
            data, bag_label = data.cuda(), bag_label.cuda()
        data, bag_label = Variable(data), Variable(bag_label)

        # reset gradients
        optimizer.zero_grad()
        # calculate loss and metrics
        loss, _ = model.calculate_objective(data, bag_label)
        train_loss += loss.item()

        # backward pass
        loss.backward()
        # step
        optimizer.step()
        #meep
        empty_cache()

    # calculate loss and error for epoch
    train_loss /= len(train_loader)
    print('Epoch: {}, Loss: {:.4f}'.format(epoch, train_loss))
    if not args.train_only:
        test()
    if train_loss<args.cutoff:
        print([train_loss,args.cutoff])
        torch.save(model,args.model_prefix+'.model')
        exit(0)

patience=0.0
c_loss=1000.0
def test():
    global patience
    global c_loss
    model.eval()
    test_loss = 0.
    test_error = 0.
    for batch_idx, (data, label, n) in enumerate(test_loader):
        bag_label = label[0]
        if args.cuda:
            data, bag_label = data.cuda(), bag_label.cuda()
        data, bag_label = Variable(data), Variable(bag_label)
        loss, attention_weights = model.calculate_objective(data, bag_label)
        test_loss += loss.item()
        
        #n=n.data.cpu().numpy()
        #label=label[0].data.cpu().numpy()
        #predicted_label=predicted_label.data.cpu().numpy()
        #print(str(n[0])+"\t"+str(label)+"\t"+str(int(predicted_label[0][0])))
        
    test_loss /= len(test_loader)
    if c_loss>test_loss:
        c_loss=test_loss
        torch.save(model,args.model_prefix+'.model')
        print(' Test Set, Loss: {:.4f}'.format(test_loss))
        patience=0
    else:
        patience=patience+1
    if patience>args.patience:
        exit(0)


if __name__ == "__main__":
    print('Start Training')
    print('# of training bags:\t' + str(len(train_loader)))
    if not args.train_only:
        print('# of testing bags:\t' + str(len(test_loader)))
    for epoch in range(1, args.epochs + 1):
        train(epoch)
    torch.save(model,args.model_prefix+'.model')
