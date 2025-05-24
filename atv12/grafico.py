import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the CSV file
try:
    df = pd.read_csv("efficiency_results.csv")
except FileNotFoundError:
    print("Error: 'efficiency_results.csv' not found. Make sure the file is in the same directory.")
    exit()

# Pivot the DataFrame to create a matrix for the heatmap
# 'LoadScale' will be columns (x-axis), 'Threads' will be index (y-axis), and 'Efficiency' will be the values
heatmap_data = df.pivot(index="Threads", columns="LoadScale", values="Efficiency")

# Create the heatmap
plt.figure(figsize=(12, 8))
sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap="viridis", linewidths=.5)

# Set the title and labels
plt.title("Efficiency Heatmap (Load Scale vs. Number of Threads)", fontsize=16)
plt.xlabel("Load Scale", fontsize=12)
plt.ylabel("Number of Threads", fontsize=12)
plt.xticks(rotation=45)
plt.yticks(rotation=0)

# Ensure layout is tight
plt.tight_layout()

# Save the figure
plt.savefig("efficiency_heatmap.png", dpi=300)
print("Heatmap saved as efficiency_heatmap.png")

# Show the plot
plt.show()