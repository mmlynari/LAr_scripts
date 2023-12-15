import matplotlib.pyplot as plt
from matplotlib.text import Text

# Sample data
x = [1, 2, 3, 4]
y1 = [2, 4, 6, 8]
y2 = [1, 2, 1, 2]

# Create a figure and axis
fig, ax = plt.subplots()

# Plot data
line1, = ax.plot(x, y1, label='Line 1')
line2, = ax.plot(x, y2, label='Line 2')

# Create custom legend
legend_text_before = "Text before"
legend_text_after = "Text after"
legend_col1 = "Label 1"
legend_col2 = "Label 2"

# Add custom legend text
text_before = Text(0, 0, legend_text_before, ha='center', va='center')
text_col1 = Text(0, 0, f"{legend_col1}  {legend_col2}", ha='center', va='center')
text_after = Text(0, 0, legend_text_after, ha='center', va='center')

# Add legend items and text to the axis
legend = ax.legend(handles=[line1, line2], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False)
ax.add_artist(text_before)
ax.add_artist(text_col1)
ax.add_artist(text_after)

# Hide default legend
ax.legend().set_visible(False)

# Adjust layout to make room for the legend
plt.subplots_adjust(bottom=0.3)

# Show the plot
plt.show()
plt.savefig("test_legend.png")