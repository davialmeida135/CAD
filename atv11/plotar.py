import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Ler arquivo CSV
df = pd.read_csv("simulation_results.csv")

sns.set(style="whitegrid")
# Criar gráfico de dispersão
plt.figure(figsize=(10, 6))
sns.lineplot(data=df, x="TimeStep", y="Field_U_Center", palette="deep", alpha=0.7)
plt.title("Perturbação no centro vs Tempo")
plt.xlabel("Tempo")
plt.ylabel("Perturbação no centro")
# Set ticks for x and y axes to show every number
plt.yticks(ticks=np.sort(range(0,12)))

plt.grid(True)
plt.savefig("plot.png", dpi=300)
plt.show()