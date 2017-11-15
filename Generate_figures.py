# author = 'Helio WANG'
# Last edited : 06/2016

import numpy as np
import matplotlib.pyplot as plt
import sys
name = sys.argv[1]+'.txt'
size = int(name.split('size')[1][0:2])
print("Generating images from "+ name)

t=1
with open('History/'+name) as f: 
    while 1:
        mat = np.array(map(float, f.readline().split()))
        try :
            mat = mat.reshape((size,size))
        except :
            break
        fig = plt.gcf()
        img = plt.imshow(mat, clim=[-1, 2], interpolation='nearest')
        #plt.colorbar(img, ticks=[-1, 0, 1, 2])
        if t<10 :	
            nb = '0'+str(t)
        else :
            nb = str(t)
        plt.savefig('Figures/'+name+'_fig_'+nb+'.png')
        # plt.show()
        fig.clear()
        t+=1
print(str(t-1)+ " images were generated.")