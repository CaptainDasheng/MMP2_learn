#!/usr/bin/env python 
from matplotlib import pyplot as plt

data={"Cs_4ZnSnAg_2Cl_12": -62.81720875,
"ZnCl_2": -16.79261651/2.0,
"SnCl_4": -63.67686922/4.0,
"Cs_2ZnCl_4": -89.21373208/4,
"Cs_2SnCl_6": -30.373410446,
"CsAgCl_2": -24.62158239/2.0,
"Cs_2AgCl_3": -69.08917582/4.0,
"CsCl": -6.63994848,
"AgCl": -22.13092452/4.0} 

path1= {'CsCl':4, "ZnCl_2": 1, "SnCl_4":1, "AgCl":1}
path2= {"Cs_2AgCl_3":2,"ZnCl_2":1, "SnCl_4":1}
path3= {"CsAgCl_2":2,"CsCl":2,"ZnCl_2":1,"SnCl_4":1}
path4= {"CsAgCl_2":2,"Cs_2ZnCl_4":1,"SnCl_4":1}
path5= {"Cs_2SnCl_6":1,"Cs_2AgCl_3":1,"AgCl":1,"ZnCl_2":1}
path6= {"Cs_2SnCl_6":1,"CsAgCl_2":2,"ZnCl_2":1}


fig = plt.figure(figsize=(8,7))
ax = plt.subplot(111)
n = 0 
for path in [path1, path2, path3, path4, path5,path6]:
    energy = 0.0
    for p in path:
	energy  += data[p]*path[p]
    if n ==3 or n ==5 : 
        cls = 'r'
        cl = 'b'
    else:
        cls = 'g'
        cl = 'k'
    fe = (energy -  data["Cs_4ZnSnAg_2Cl_12"])*1000.0/(4+1+1+2+12)
    dpath = ' + '.join(c for c in path)
    ax.bar(n, fe, 0.8, color=cls, alpha=0.6)
    ax.text(n, 220, '$\mathregular{'+dpath+'}$',color=cl,rotation=90, fontsize=14).set_weight("bold") 
    n += 1
ax.axes.set_linewidth=2.0
ax.axhline(linestyle='--',color='r')
ax.set_xticks([])
#ax.set_xtickslabel(fontsize=12)
#ax.set_ytickslabel(fontsize=12)
plt.title("$\mathregular{Cs_{4}[Ag_{2}ZnSn]Cl_{12}}$",fontsize=16)
plt.ylabel("$\mathregular{\Delta{E} \ (meV/atom)}$",fontsize=16)
plt.xlabel("Decomposition pathways", fontsize=16)
plt.savefig('decomp.png', dpi=400)
#plt.show()

