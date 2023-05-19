
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from pytools import natsorted

list1 = [i.strip() for i in open("Chr1.ordered.list") if i.strip()]
list2 = natsorted(list1)

real_db = dict(zip(list2, range(len(list2))))

list1 = list(map(lambda x: real_db[x], list1))
list2 = list(map(lambda x: real_db[x], list2))

fig, ax = plt.subplots(figsize=(7, 7))
plt.scatter(list2, list1, color='black', s=1)
plt.xlabel('C-Phasing', fontsize=18)
plt.ylabel('Real', fontsize=18)
plt.savefig("Chr1.ordered.png", dpi=600)