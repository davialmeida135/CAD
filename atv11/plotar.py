import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Ler arquivo CSV
df = pd.read_csv("resultado.csv")

# Normalizar cor de acordo com u
u_min, u_max = df['u'].min(), df['u'].max()
colors = (df['u'] - u_min) / (u_max - u_min)
df['deviation'] = df['u'] - 1.0
colors = df['deviation']
# Plot 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(df['i'], df['j'], df['k'], c=colors, cmap='coolwarm')

# Rótulos e cor
ax.set_xlabel("i")
ax.set_ylabel("j")
ax.set_zlabel("k")
fig.colorbar(sc, label='u (velocidade)')
plt.title("Distribuição do campo de velocidades (u) em 3D")
plt.show()
