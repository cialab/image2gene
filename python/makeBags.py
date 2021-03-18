from __future__ import print_function

import numpy as np
import argparse
from dataloader import SlideDataset
from slideModel import Attention
import openslide as opsl
import scipy.io
import h5py as h5
import random
import math
import pdb

parser = argparse.ArgumentParser()
parser.add_argument('--patch_size', type=int, metavar='PS',
                    help='patch size')               
parser.add_argument('--slide_path',
                    help='which slide to process')
parser.add_argument('--mask_path',
                    help='where is associated mask')
parser.add_argument('--model_path', default='saved.model',
                    help='which model you want to use')
parser.add_argument('--num_query', type=int, default=0,
                    help='override default query points calculatation')
parser.add_argument('--num_instances', type=int, default=1000)
parser.add_argument('--method', type=int, default=0)
parser.add_argument('--centers', type=float, default=3)
parser.add_argument('--rd')
parser.add_argument('--wd')

args = parser.parse_args()

def loadmat():
    #fn=args.mask_path.split("/")[-1].split(".")[0] # train
    fn=args.mask_path.split("/")[-1] # test
    gn=args.model_path.split("/")[-1].split('.')[0].split('_')[0]
    m=scipy.io.loadmat(args.rd+gn+'/'+fn+'.mat')
    weights=m['weights2']
    weights=np.swapaxes(weights,0,1)
    bag_points=m['bag_points']
    if bag_points.shape[0]==1:
        bag_points=np.concatenate(np.concatenate(bag_points,axis=0),axis=0)
    else:
        bag_points=np.concatenate(bag_points,axis=0)
    return bag_points,weights
    
def readimg(minx,miny,maxx,maxy,bs=4096):
    f=opsl.OpenSlide(args.slide_path)
    im=np.zeros((maxy-miny,maxx-minx,3),dtype=np.uint8)
    xs=np.append(np.arange(minx,maxx,bs),[maxx],axis=0)
    ys=np.append(np.arange(miny,maxy,bs),[maxy],axis=0)
    for i in range(0,len(xs)-1):
        for j in range(0,len(ys)-1):
            img=f.read_region((xs[i],ys[j]),0,(xs[i+1]-xs[i],ys[j+1]-ys[j]))
            im[ys[j]-miny:ys[j+1]-miny,xs[i]-minx:xs[i+1]-minx,:]=np.array(img)[:,:,0:3]
    return im
    
def makeBag(im,p,minx,miny,maxx,maxy,n):
    x=np.zeros((args.patch_size,args.patch_size,3,p.shape[0]),dtype=np.uint8)
    ii=0
    iii=0
    for i in range(0,p.shape[0]):
        img=im[p[i,1]-miny:p[i,1]-miny+args.patch_size,p[i,0]-minx:p[i,0]-minx+args.patch_size,:]
        x[:,:,:,i]=img
    return x
    
def rgb2gray(rgb):
    return np.dot(rgb[...,:3], [0.2989, 0.5870, 0.1140])

def makeBags_segments(bp,w,im,minx,miny,maxx,maxy):

    # Process by segments
    m=100.0/args.centers
    bags=[]
    for i in range(0,int(args.centers)+1):
        center=m*i
        idx=np.argsort(np.abs(np.subtract(np.percentile(w,center),w)),axis=0)[0:int(args.num_instances/args.centers)]
        bag=makeBag(im,np.squeeze(bp[idx]),minx,miny,maxx,maxy,int(args.num_instances/args.centers))
        bags.append(bag)

    bag=np.concatenate(bags,axis=3)
    bag=np.transpose(bag)

    fn=args.mask_path.split("/")[-1].split(".")[0]
    #gn=args.model_path.split("/")[-1].split('.')[0].split('_')[0]
    f=h5.File(args.wd+"/"+fn+".h5","w")
    dset=f.create_dataset('bag0',data=bag,dtype='u1')
    f.close()

def makeBags_random(bp,w,im,minx,miny,maxx,maxy):

    # Number of weights
    n=w.shape[0]
    
    # Random
    if args.num_instances>n:
        idx=range(0,n)
    else:
        idx=random.sample(range(0,n),args.num_instances)
    bag=makeBag(im,np.squeeze(bp[idx]),minx,miny,maxx,maxy,args.num_instances)
    bag=np.transpose(bag)

    # Write bags
    fn=args.mask_path.split("/")[-1].split(".")[0]
    #gn=args.model_path.split("/")[-1].split('.')[0].split('_')[0]
    f=h5.File(args.wd+"/"+fn+".h5","w")
    dset=f.create_dataset('bag0',data=bag,dtype='u1')
    f.close()

def makeBags_uniform(bp,w,im,minx,miny,maxx,maxy):

    # Number of weights
    n=w.shape[0]
    
    # Uniform
    idx=np.round(np.linspace(0,n-1,args.num_instances))
    bag=makeBag(im,np.squeeze(bp[idx.astype(np.uint64)]),minx,miny,maxx,maxy,args.num_instances)
    bag=np.transpose(bag)
    
    # Make and write bags
    fn=args.mask_path.split("/")[-1].split(".")[0]
    #gn=args.model_path.split("/")[-1].split('.')[0].split('_')[0]
    f=h5.File(args.wd+"/"+fn+".h5","w")
    dset=f.create_dataset('bag0',data=bag,dtype='u1')
    f.close()

def preProcess():
    bp,w=loadmat()

    # Read portion of slide
    minx=np.min(bp[:,0])
    miny=np.min(bp[:,1])
    maxx=np.max(bp[:,0])
    maxy=np.max(bp[:,1])
    im=readimg(minx,miny,maxx,maxy)

    # Throw out whitespace
    keep=np.zeros((bp.shape[0]),dtype=np.bool)
    for i in range(0,bp.shape[0]):
        img=im[bp[i,1]-miny:bp[i,1]-miny+args.patch_size,bp[i,0]-minx:bp[i,0]-minx+args.patch_size,:]
        if np.sum(rgb2gray(img)<230)/(args.patch_size*args.patch_size)>0.50:
            keep[i]=True
    bp=bp[keep,]
    w=w[keep]

    return bp,w,im,minx,miny,maxx,maxy

if __name__ == "__main__":
    bp,w,im,minx,miny,maxx,maxy=preProcess()

    if args.method==0:
        makeBags_random(bp,w,im,minx,miny,maxx,maxy)
    elif args.method==1:
        makeBags_uniform(bp,w,im,minx,miny,maxx,maxy)
    elif args.method==2:
        makeBags_segments(bp,w,im,minx,miny,maxx,maxy)
    else:
        print('Must choose method')
        exit(0)
