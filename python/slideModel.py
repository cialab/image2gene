import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import pdb

class Attention(nn.Module):
    def __init__(self,ps):
        super(Attention, self).__init__()
        self.L = 512
        self.D = 128
        self.K = 1
        self.P = 10*(((((ps+2-6)/2)/2+1)-2)/2)*(((((ps+2-6)/2)/2+1)-2)/2) # this is the size of the last three dimensions of the output of the feature extraction step

        self.feature_extractor_part1 = nn.Sequential(
            nn.Conv2d(3, 48, kernel_size=7, stride=1, padding=2),
            nn.ReLU(),
            nn.MaxPool2d(2, stride=2),
            nn.Conv2d(48, 48, kernel_size=1, stride=1, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(2, stride=2),
            nn.Conv2d(48, 10, kernel_size=3, stride=1, padding=0),
            nn.ReLU(),
            nn.MaxPool2d(2, stride=2)
        )

        self.feature_extractor_part2 = nn.Sequential(
            nn.Linear(self.P, self.L),
            nn.ReLU(),
            nn.Dropout(),
            nn.Linear(self.L,self.L),
            nn.ReLU(),
            nn.Dropout()
        )

        self.attention = nn.Sequential(
            nn.Linear(self.L, self.D),
            nn.Tanh(),
            nn.Linear(self.D, self.K)
        )

        self.regressor = nn.Sequential(
            nn.Linear(self.L*self.K, 1),
        )

    def forward(self, x):
        x = x.squeeze(0)

        H = self.feature_extractor_part1(x)
        H = H.view(-1, self.P) # flatten each instance
        H = self.feature_extractor_part2(H)  # NxL; does both fc-512 layers
        
        A = self.attention(H)  # NxK
        A = torch.transpose(A, 1, 0)  # KxN
        A = F.softmax(A, dim=1)  # softmax over N

        M = torch.mm(A, H)  # KxL

        Y_hat = self.regressor(M)

        return Y_hat, A

    # AUXILIARY METHODS
    def calculate_objective(self, X, Y):
        Y = Y.float()
        Y_hat, A = self.forward(X)
        error = torch.pow(Y-Y_hat,2)

        return error, A
        
    def attentionWeights(self, x):
        x = x.squeeze(0)

        H = self.feature_extractor_part1(x)
        H = H.view(-1, self.P) # flatten each instance
        H = self.feature_extractor_part2(H)  # NxL; does both fc-512 layers
        
        A = self.attention(H)  # NxK
        A = torch.transpose(A, 1, 0)  # KxN

        return A
        
class MultiRegressionAttention(nn.Module):
    def __init__(self,ps,nout):
        super(MultiRegressionAttention, self).__init__()
        self.L = 90
        self.D = 32
        self.K = 1
        self.P = 10*(((((ps+2-6)/2)/2+1)-2)/2)*(((((ps+2-6)/2)/2+1)-2)/2) # this is the size of the last three dimensions of the output of the feature extraction step
        self.nout=nout
        self.cos = nn.CosineSimilarity(dim=1, eps=1e-6)

        self.feature_extractor_part1 = nn.Sequential(
            nn.Conv2d(3, 48, kernel_size=7, stride=1, padding=2),
            nn.ReLU(),
            nn.MaxPool2d(2, stride=2),
            nn.Conv2d(48, 48, kernel_size=1, stride=1, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(2, stride=2),
            nn.Conv2d(48, 10, kernel_size=3, stride=1, padding=0),
            nn.ReLU(),
            nn.MaxPool2d(2, stride=2)
        )

        self.feature_extractor_part2 = nn.Sequential(
            nn.Linear(self.P, self.L),
            nn.ReLU(),
            nn.Dropout(),
            nn.Linear(self.L,self.L),
            nn.ReLU(),
            nn.Dropout()
        )

        self.attention = nn.Sequential(
            nn.Linear(self.L, self.D),
            nn.Tanh(),
            nn.Linear(self.D, self.K)
        )

        self.regressor = nn.Sequential(
            nn.Linear(self.L*self.K, self.nout),
        )

    def forward(self, x):
        x = x.squeeze(0)

        H = self.feature_extractor_part1(x)
        H = H.view(-1, self.P) # flatten each instance
        H = self.feature_extractor_part2(H)  # NxL; does both fc-512 layers
        
        A = self.attention(H)  # NxK
        A = torch.transpose(A, 1, 0)  # KxN
        A = F.softmax(A, dim=1)  # softmax over N

        M = torch.mm(A, H)  # KxL
        
        # Legendre polynomial n=2
        # M = torch.sub(torch.mul(1.5,torch.pow(M,2)),0.5)
        
        # Legendre polynomial n=3
        # M = torch.sub(torch.mul(2.5,torch.pow(M,3)),torch.mul(3,M))

        Y_hat = self.regressor(M)

        return Y_hat, A

    # AUXILIARY METHODS
    def calculate_objective(self, X, Y):
        Y = torch.transpose(Y.float(),1,0)
        Y_hat, A = self.forward(X)
        #error = 1.-self.cos(Y,Y_hat)
        error = torch.mean(torch.pow(Y-Y_hat,2))

        return error, A
        
    def attentionWeights(self, x):
        x = x.squeeze(0)

        H = self.feature_extractor_part1(x)
        H = H.view(-1, self.P) # flatten each instance
        H = self.feature_extractor_part2(H)  # NxL; does both fc-512 layers
        
        A = self.attention(H)  # NxK
        A = torch.transpose(A, 1, 0)  # KxN

        return A
        
class MultiRegressionAttention2(nn.Module):
    def __init__(self,ps,nout):
        super(MultiRegressionAttention2, self).__init__()
        self.L = 90
        self.D = 32
        self.K = 1
        self.P = 10*(((((ps+2-6)/2)/2+1)-2)/2)*(((((ps+2-6)/2)/2+1)-2)/2) # this is the size of the last three dimensions of the output of the feature extraction step
        self.nout=nout
        self.cos = nn.CosineSimilarity(dim=1, eps=1e-6)

        self.feature_extractor_part1 = nn.Sequential(
            nn.Conv2d(3, 48, kernel_size=7, stride=1, padding=2),
            nn.ReLU(),
            nn.MaxPool2d(2, stride=2),
            nn.Conv2d(48, 48, kernel_size=1, stride=1, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(2, stride=2),
            nn.Conv2d(48, 10, kernel_size=3, stride=1, padding=0),
            nn.ReLU(),
            nn.MaxPool2d(2, stride=2)
        )

        self.feature_extractor_part2 = nn.Sequential(
            nn.Linear(self.P, self.L),
            nn.ReLU(),
            nn.Dropout(),
            nn.Linear(self.L,self.L),
            nn.ReLU(),
            nn.Dropout()
        )

        self.attention = nn.Sequential(
            nn.Linear(self.L, self.D),
            nn.Tanh(),
            nn.Linear(self.D, self.K)
        )

        self.regressor = nn.Sequential(
            nn.Linear(self.L*self.K, self.nout),
        )

    def forward(self, x):
        x = x.squeeze(0)

        H = self.feature_extractor_part1(x)
        H = H.view(-1, self.P) # flatten each instance
        H = self.feature_extractor_part2(H)  # NxL; does both fc-512 layers
        
        A = self.attention(H)  # NxK
        A = torch.transpose(A, 1, 0)  # KxN
        A = F.softmax(A, dim=1)  # softmax over N

        M = torch.mm(A, H)  # KxL
        
        # Legendre polynomial n=2
        # M = torch.cat((M,torch.sub(torch.mul(1.5,torch.pow(M,2)),0.5)),0)
        
        # Legendre polynomial n=3
        # M = torch.sub(torch.mul(2.5,torch.pow(M,3)),torch.mul(3,M))

        Y_hat = self.regressor(M)

        return Y_hat, A

    # AUXILIARY METHODS
    def calculate_objective(self, X, Y):
        Y = torch.transpose(Y.float(),1,0)
        Y_hat, A = self.forward(X)
        #error = 1.-self.cos(Y,Y_hat)
        error = torch.mean(torch.pow(Y-Y_hat,2))

        return error, A
        
    def attentionWeights(self, x):
        x = x.squeeze(0)

        H = self.feature_extractor_part1(x)
        H = H.view(-1, self.P) # flatten each instance
        H = self.feature_extractor_part2(H)  # NxL; does both fc-512 layers
        
        A = self.attention(H)  # NxK
        A = torch.transpose(A, 1, 0)  # KxN

        return A

    def predictFullBag(self,x,n):
        x=x.squeeze(0)
        l=np.round(np.linspace(0,x.shape[0],np.ceil(float(x.shape[0])/float(n))+1)).astype(np.uint64)
        Hs=[]
        As=[]
        
        for i in range(0,len(l)-1):
            H = self.feature_extractor_part1(x[l[i]:l[i+1]])
            H = H.view(-1, self.P) # flatten each instance
            H = self.feature_extractor_part2(H)
            Hs.append(H)

            A = self.attention(H)  # NxK
            A = torch.transpose(A, 1, 0)  # KxN
            As.append(A)

        Hs=torch.cat(Hs,0)
        As=torch.cat(As,1)
        As = F.softmax(As, dim=1)
        M = torch.mm(As, Hs)

        Y_hat = self.regressor(M)

        return Y_hat

    def patchPredict(self, x):
        x = x.squeeze(0)

        H = self.feature_extractor_part1(x)
        H = H.view(-1, self.P) # flatten each instance
        
        # Legendre polynomial n=2
        # M = torch.cat((H,torch.sub(torch.mul(1.5,torch.pow(H,2)),0.5)),0)
        
        # Legendre polynomial n=3
        # M = torch.sub(torch.mul(2.5,torch.pow(M,3)),torch.mul(3,M))

        Y_hat = self.regressor(H)

        return Y_ha