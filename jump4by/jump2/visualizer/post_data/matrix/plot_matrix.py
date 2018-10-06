#!/usr/bin/env python 

from matplotlib import pyplot as plt 


def get_data(inf):

    data = {}
    with open(inf, 'r') as f:

        while True:
            line = f.readline()
            if not line:
                break
            else:
                line = f.readline()
                nk = (line.split()[-1]).strip(':')
                line = f.readline()
                if nk not in data:
                    data[nk] = {}
                    data[nk]['x'] = []
                    data[nk]['y'] = []
                energy  = float(line.split()[2]) 
                ampmx = float(line.split()[5])
                data[nk]['x'].append(energy)
                data[nk]['y'].append(ampmx)
                
    return data 
 
fig=plt.figure(figsize=(6,4))
ax1=plt.subplot(121)
ax2=plt.subplot(122)
data1 = get_data('AgIn.matrix')
data2 = get_data('ZnSn.matrix')

for k in data1:
    ax1.bar(data1[k]['x'], data1[k]['y'], color='k', width=0.01)
for k in data2:
    ax2.bar(data2[k]['x'], data2[k]['y'], color='k', width=0.01)

ax1.arrow(1.03,0.08, 0.0, -0.06, fc='r', ec='r',
	head_width=0.05, head_length=0.02)
ax1.text(1.03-0.3, 0.085, "$\mathregular{E_{g}="+str(1.03)+"}$")
ax2.arrow(0.81,0.08, 0.0, -0.06, fc='r', ec='r',
	head_width=0.05, head_length=0.02)
ax2.text(0.81-0.3, 0.085, "$\mathregular{E_{g}="+str(0.81)+"}$")
ax1.set_title('$\mathregular{Cs_{4}[AgIn]_{2}Cl_{12}}$')
ax2.set_title('$\mathregular{Cs_{4}[Ag_{2}ZnSn]Cl_{12}}$')
ax1.set_xlim(0.5,3.0)
ax1.set_ylim(0.,0.6)
ax2.set_ylim(0.,0.6)
ax2.set_xlim(0.5,3.0)
plt.xlabel("Energy (eV)")
ax1.set_ylabel("Ampitude")
plt.savefig('matrix.png', dpi=300)
#plt.show()



