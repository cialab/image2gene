import glob
import numpy
import sys
from subprocess import Popen
import subprocess
import os
import numpy as np
import random

d='/isilon/datalake/cialab/scratch/cialab/tet/python/AttentionDeepMIL/image2gene/revision/run_makeBags_test/*.sh'
p=glob.glob(d)
hey=np.linspace(0,len(p),10+1)
hey=list(np.round(hey).astype(np.int))

# Figure hostname
n=int(sys.argv[1])#int(sys.argv[1])

a=[]
limit=5
q=0
i=hey[n]
iend=hey[n+1]

print(str(i)+":"+str(iend))
print('Total jobs: '+str(iend-i))
while True:
    if q<limit:
        pr=Popen(['sh',p[i]])
        a.append(pr)
        q=q+1
        i=i+1
        print(q,i)
    # Check if any process has ended
    for ii in range(len(a)-1,-1,-1):
        if a[ii].poll() is not None:
            a.pop(ii)
            q=q-1
    if i==iend:
        while len(a)>0:
            for p in range(len(a)-1,-1,-1):
                if a[p].poll() is not None:
                    print('Popping '+str(p))
                    a.pop(p)
                    q=q-1
        break
