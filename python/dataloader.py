import h5py
#import helpers
import numpy as np
from pathlib import Path
import torch
from torch.utils import data
from datetime import datetime
import random
import openslide as opsl
import pdb
        
class SlideDataset(data.Dataset):
    """.
    
    Input params:
        slide_path: where the slide is located
        bag_points: list of sets of points for each bag
        load_data: If True, loads all the data immediately into RAM. Use this if
            the dataset is fits into memory. Otherwise, leave this at false and 
            the data will load lazily.
        transform: PyTorch transform to apply to every data instance (default=None).
    """
    def __init__(self, slide_path, bag_points, load_data, transform, patch_size, minx, miny, maxx, maxy):
        super(SlideDataset,self).__init__()
        self.bag_points=bag_points
        self.load_data=load_data
        self.transform=transform
        self.patch_size=patch_size
        self.minx=minx
        self.miny=miny
        self.maxx=maxx
        self.maxy=maxy
        
        f=opsl.OpenSlide(slide_path)
        print('Reading image...')
        self.im=self.readit(f,maxx,minx,maxy,miny)
        print('Done')
        
    def __getitem__(self, index):
        p=self.bag_points[index]
        x=np.zeros((p.shape[0],3,self.patch_size,self.patch_size),dtype=np.uint8)
        for i in range(0,p.shape[0]):
            im=self.im[p[i,1]-self.miny:p[i,1]-self.miny+self.patch_size,p[i,0]-self.minx:p[i,0]-self.minx+self.patch_size,:]
            if im.shape[0]!=32 | im.shape[1]!=32:
                print(im.shape)
                print([self.maxx-self.minx,self.maxy-self.miny,p[i,1]-self.miny,p[i,0]-self.minx,self.im.shape])
                print('fuck')
            x[i,:,:,:]=np.swapaxes(im,0,2)
        
        if self.transform:
            x = np.float32(x)
            x = torch.from_numpy(x/127.5-1.)
        else:
            x = torch.from_numpy(np.float32(x))

        return x

    def __len__(self):
        return len(self.bag_points)
        
    def readit(self,f,maxx,minx,maxy,miny,bs=1600):
        im=np.zeros((maxy-miny,maxx-minx,3),dtype=np.uint8)
        xs=np.append(np.arange(minx,maxx,bs),[maxx],axis=0)
        ys=np.append(np.arange(miny,maxy,bs),[maxy],axis=0)
        for i in range(0,len(xs)-1):
            for j in range(0,len(ys)-1):
                img=f.read_region((xs[i],ys[j]),0,(xs[i+1]-xs[i],ys[j+1]-ys[j]))
                im[ys[j]-miny:ys[j+1]-miny,xs[i]-minx:xs[i+1]-minx,:]=np.array(img)[:,:,0:3]
        return im

class HDF5Dataset(data.Dataset):
    """Represents an abstract HDF5 dataset.
    
    Input params:
        file_path: Path to the folder containing the dataset (one or multiple HDF5 files).
        recursive: If True, searches for h5 files in subdirectories.
        load_data: If True, loads all the data immediately into RAM. Use this if
            the dataset is fits into memory. Otherwise, leave this at false and 
            the data will load lazily.
        data_cache_size: Number of HDF5 files that can be cached in the cache (default=3).
        transform: PyTorch transform to apply to every data instance (default=None).
    """
    def __init__(self, file_path, load_data, train, seed, split, gene_idx, num_instances, data_cache_size=3, transform=False):
        super(HDF5Dataset,self).__init__()
        self.data_info = []
        self.data_cache = {}
        self.data_cache_size = data_cache_size
        self.transform = transform
        self.idx=0
        self.train=train
        self.gene_idx=gene_idx
        self.num_instances=num_instances
        self.gene_vals=[]
        self.onegene=isinstance(self.gene_idx,int)

        # Search for all h5 files
        p = Path(file_path)
        bags = sorted(p.glob('*.h5'))
        
        # For consistent folds
        random.seed(seed)
        random.shuffle(bags)
        random.seed(datetime.now()) # assuming this makes things non-deterministic again

        # split is the number testing slides per class
        split=int(split)
        l=np.round(np.linspace(0,len(bags),11)).astype(np.uint64)
        if split==-1:
            files=bags
        elif train:
            bags=bags[l[0]:l[split]]+bags[l[split+1]:]
            files=bags
        else:
            bags=bags[l[split]:l[split+1]]
            files=bags

        for h5dataset_fp in files:
            a=1
            while a==1:
                try:
                    self._add_data_infos(str(h5dataset_fp.resolve()), load_data)
                    a=0
                except:
                    print('collision. trying again.')
                    a=1
        
        # Sampling if one gene
        if self.onegene:
            ming=np.min(np.floor(np.multiply(self.gene_vals,2.0)))/2
            maxg=np.max(np.ceil(np.multiply(self.gene_vals,2.0)))/2
            n=np.linspace(ming,maxg,(maxg-ming)/0.5+1)
            self.bins=[ [] for _ in range(len(n)-1) ]
            for i in range(0,len(self.data_info)):
                gene_val=self.data_info[i]['gene_exp'][self.gene_idx][0]
                for j in range(0,len(n)-1):
                    if (gene_val>=n[j]) & (gene_val<n[j+1]):
                        self.bins[j].append(i)
                        break
        
    def __getitem__(self, index):
        # get random data; ignore index
        if self.train:
            if self.onegene:
                i=random.randint(0,len(self.bins)-1)
                while not self.bins[i]:
                    i=random.randint(0,len(self.bins)-1)
                i=random.choice(self.bins[i])
            else:
                i=random.choice(range(0,len(self)))
        # for testing
        else:
            i=index

        x=self.get_data('bag',i)
        if self.transform:
            x = np.float32(x)
            x = torch.from_numpy(x/127.5 - 1.)
        else:
            x = torch.from_numpy(np.float32(x))
        
        # sample it
        if self.num_instances>0:
            x=x[np.random.choice(x.shape[0],self.num_instances),:,:,:]
        
        # get gene expression value
        if self.gene_idx==-1:
            # no gene expression value
            y = self.data_info[i]['label']
            y = torch.from_numpy(np.array(y))
        else:
            y = self.data_info[i]['gene_exp'][self.gene_idx]
            y = torch.from_numpy(np.array(y))

        # get mouse number
        n = self.data_info[i]['mousenum']
        
        return (x, y, n)

    def __len__(self):
        return len(self.get_data_infos('bag'))
    
    def _add_data_infos(self, file_path, load_data):
        with h5py.File(file_path) as h5:
            if len(h5.keys())==2:
                # This is a slide without gene expression data (classification)
                label=h5['label'].value
                mousenum=int(file_path.split('/')[-1].split('.')[0])
                idx = -1
                if load_data:
                        # add data to the data cache
                    idx = self._add_to_cache(h5['bag'].value, file_path)
                    
                # Append the bag
                self.data_info.append({'file_path': file_path, 'type': 'bag', 'shape': h5['bag'].shape, 'cache_idx': idx, 'label': label, 'mousenum': mousenum})
            elif len(h5.keys())==3:
                # One big bag
                label=h5['label'].value
                gene_exp=h5['gene_exp'].value
                self.gene_vals.append(h5['gene_exp'].value[self.gene_idx][0])
                mousenum=int(file_path.split('/')[-1].split('.')[0])
                for k in h5.keys():
                    # It's a bag!
                    if "bag" in str(k):
                        shape = h5[k].shape
                        
                        # if data is not loaded its cache index is -1
                        idx = -1
                        if load_data:
                            # add data to the data cache
                            idx = self._add_to_cache(h5[k].value, file_path)
                        
                        # Append the bag
                        self.data_info.append({'file_path': file_path, 'type': 'bag', 'shape': shape, 'cache_idx': idx, 'label': label, 'mousenum': mousenum, 'gene_exp': gene_exp})
            else:
                # Get the mouse number, gene info, and label first
                label=h5['label'].value
                gene_exp=h5['gene_exp'].value
                self.gene_vals.append(h5['gene_exp'].value[self.gene_idx][0])
                mousenum=int(file_path.split('/')[-1].split('.')[0])
                for k in h5.keys():
                    # It's a bag!
                    #if "bag" in str(k):
                    if ("2" in str(k)) | ("3" in str(k)) | ("5" in str(k)) | ("6" in str(k))| ("7" in str(k)):   # this is temporary
                        i=int(k.split("bag")[1])
                        shape = h5[k].shape
                        
                        # if data is not loaded its cache index is -1
                        idx = -1
                        if load_data:
                            # add data to the data cache
                            idx = self._add_to_cache(h5[k].value, file_path)
                        
                        # Append the bag
                        self.data_info.append({'file_path': file_path, 'type': 'bag', 'num': i, 'shape': shape, 'cache_idx': idx, 'label': label, 'mousenum': mousenum, 'gene_exp': gene_exp})

    def _load_data(self, file_path):
        """Load data to the cache given the file
        path and update the cache index in the
        data_info structure.
        """
        with h5py.File(file_path) as h5_file:
            for dname, ds in h5_file.items():
                # add data to the data cache and retrieve
                # the cache index
                idx = self._add_to_cache(ds.value, file_path)

                # find the beginning index of the hdf5 file we are looking for
                file_idx = next(i for i,v in enumerate(self.data_info) if v['file_path'] == file_path)

                # the data info should have the same index since we loaded it in the same way
                self.data_info[file_idx + idx]['cache_idx'] = idx

        # remove an element from data cache if size was exceeded
        if len(self.data_cache) > self.data_cache_size:
            # remove one item from the cache at random
            removal_keys = list(self.data_cache)
            removal_keys.remove(file_path)
            self.data_cache.pop(removal_keys[0])
            # remove invalid cache_idx
            self.data_info = [{'file_path': di['file_path'], 'type': di['type'], 'shape': di['shape'], 'cache_idx': -1} if di['file_path'] == removal_keys[0] else di for di in self.data_info]

    def _add_to_cache(self, data, file_path):
        """Adds data to the cache and returns its index. There is one cache
        list for every file_path, containing all datasets in that file.
        """
        if file_path not in self.data_cache:
            self.data_cache[file_path] = [data]
        else:
            self.data_cache[file_path].append(data)
        return len(self.data_cache[file_path]) - 1

    def get_data_infos(self, type):
        """Get data infos belonging to a certain type of data.
        """
        data_info_type = [di for di in self.data_info if di['type'] == type]
        return data_info_type

    def get_data(self, type, i):
        """Call this function anytime you want to access a chunk of data from the
            dataset. This will make sure that the data is loaded in case it is
            not part of the data cache.
        """
        fp = self.get_data_infos(type)[i]['file_path']
        if fp not in self.data_cache:
            self._load_data(fp)
        
        # get new cache_idx assigned by _load_data_info
        cache_idx = self.get_data_infos(type)[i]['cache_idx']
        return self.data_cache[fp][cache_idx]