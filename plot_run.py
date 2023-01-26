from asyncore import read
from ase.io import read
import matplotlib.pyplot as plt

atoms = read("structures_-1.0.traj@:")
lenergy = 1000.0
num = 0
nums = []
engs = []
stc = []
for i in atoms:
    eng = -i.info['key_value_pairs']['raw_score']

    stc.append(i.info['key_value_pairs']['var_stc'])
    lenergy = min(eng,lenergy)
    engs.append(lenergy)
    nums.append(num)
    num+=1
   
plt.xlabel("Steps")
plt.ylabel("Most stable structure Energy")
plt.plot(nums,engs,color='green')

plt.ylabel("Stoichiometry ID")
plt.plot(nums,stc,color='red',linestyle='dashed')
plt.show()
