import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('results.csv')

N = 8
tau = 0.0125
x = 1.047198
ma = 4.5
env = 'delta'
sig = 3.0
aT = 10.0
lam = 0.0
aD_0 = 1.0
l_0_init = 0.0
cutoff = 1e-18
l_0_list = [0.4] # [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
cutoff_list = [1e-18] # [1e-18, 1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4]

for cutoff in cutoff_list: 
    for l_0 in l_0_list:

        df_tmp = df[(df.N == N) & (df.x == x) & (df.ma == ma) & (df.l_0_init == l_0_init) & (df.l_0 == l_0) & (df.tau == tau) & (df.env == env) & (df.sig == sig) & (df.aT == aT) & (df.lam == lam) & (df.aD_0 == aD_0)]

        print(df_tmp)

        steps = df_tmp.step
        pn = df_tmp.pn
        plt.plot(steps[:300], pn[:300], label = f'l_0={l_0},c={cutoff}')
    
plt.title(f'N={N},x={x},ma={ma},tau={tau},env={env},aT={aT},aD_0={aD_0}')
plt.legend()
plt.savefig('specific_plotter.png')