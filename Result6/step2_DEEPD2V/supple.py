import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
data = pd.read_csv("data2.csv")

# Check if the data needs to be transposed for horizontal orientation
# Transpose the DataFrame to swap rows and columns
data_transposed = data.T

# Set figure size for horizontal layout
plt.figure(figsize=(10, 2))  # Adjust the height to make it more horizontal
sns.set(font_scale=0.7)

# Create the heatmap with adjusted parameters
ax = sns.heatmap(data_transposed, square=False, cbar_kws={'orientation': 'horizontal', "shrink": 0.3},
                 annot=True, annot_kws={"size": 8}, vmax=1.0, vmin=-1.0, cmap='coolwarm')

# Set title
ax.set(title="House SalePrice Correlation Heatmap")

# Rotate x-axis labels if necessary
ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha='center')  # Horizontal and centered
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha='right')

# Adjust the tick parameters
ax.tick_params(axis='x', labelsize=7)
ax.tick_params(axis='y', labelsize=7)

# Show the plot
plt.show()
