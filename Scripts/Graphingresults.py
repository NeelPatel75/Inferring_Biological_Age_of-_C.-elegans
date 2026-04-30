# This script creates a bar plot comparing corrected biological ages across sample groups
# It reads the final biological age results CSV created by the main pipeline
# It separates samples into WT embryo, WT L1, xol-1 embryo, and xol-1 L1 groups.
# It plots the mean corrected biological age for each group with SEM error bars
# It also runs Welch's t-tests only to decide whether to add significance asterisks to the graph.
# The final plot is saved as a PNG file in the Results folder

import os  # Imports os so the script can create file paths that work on different computers.
import pandas as pd  # Imports pandas to read and work with the CSV data.
import matplotlib.pyplot as plt  # Imports matplotlib to create and save the bar plot.
from scipy.stats import ttest_ind  # Imports Welch's t-test function for significance annotations.

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))  # Gets the folder where this graphing script is located
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)  # Gets the main project folder by going one level above Scripts.
RESULTS_DIR = os.path.join(PROJECT_ROOT, "Results")  # Creates the path to the Results folder

INPUT_FILE = os.path.join(RESULTS_DIR, "Final_Biological_Age_Results.csv")  # Sets the input CSV file created by the main analysis script.
PLOT_FILE = os.path.join(RESULTS_DIR, "Corrected_Biological_Age_Comparisons.png")  # Sets the output path for the saved plot.

df = pd.read_csv(INPUT_FILE)  # Reads the final biological age results into a pandas DataFrame.

sample_col = "Sample"  # Defines the column containing the sample names/SRR IDs.
age_col = "Corrected_Biological_Age"  # Defines the column containing corrected biological age values.

wt_embryo = ["SRR28479540", "SRR28479539", "SRR28479538", "SRR28479525", "SRR28479526", "SRR28479527"]  # Defines WT embryo sample IDs.
wt_l1 = ["SRR28479534", "SRR28479533", "SRR28479532", "SRR28479531"]  # Defines WT L1 sample IDs.
xol1_embryo = ["SRR28479537", "SRR28479536", "SRR28479535"]  # Defines xol-1 embryo sample IDs.
xol1_l1 = ["SRR28479530", "SRR28479529", "SRR28479528"]  # Defines xol-1 L1 sample IDs.

all_embryo = wt_embryo + xol1_embryo  # Combines WT embryo and xol-1 embryo samples into one embryo group
all_l1 = wt_l1 + xol1_l1  # Combines WT L1 and xol-1 L1 samples into one L1 group.
wt_all = wt_embryo + wt_l1  # Combines all WT samples into one WT group.
xol1_all = xol1_embryo + xol1_l1  # Combines all xol-1 samples into one xol-1 group

def get_group_ages(df, srr_list):  # Defines a function that extracts corrected ages for a specific list of SRR IDs.
    mask = df[sample_col].apply(lambda x: any(str(x).startswith(srr) for srr in srr_list))  # Finds rows where the Sample column starts with one of the SRR IDs
    return df[mask][age_col].dropna()  # Returns only the corrected biological ages for matching samples and removes missing values

ages_wt_embryo = get_group_ages(df, wt_embryo)  # Extracts corrected biological ages for WT embryo samples.
ages_wt_l1 = get_group_ages(df, wt_l1)  # Extracts corrected biological ages for WT L1 samples.
ages_xol1_embryo = get_group_ages(df, xol1_embryo)  # Extracts corrected biological ages for xol-1 embryo samples.
ages_xol1_l1 = get_group_ages(df, xol1_l1)  # Extracts corrected biological ages for xol-1 L1 samples.

ages_all_embryo = get_group_ages(df, all_embryo)  # Extracts corrected biological ages for all embryo samples.
ages_all_l1 = get_group_ages(df, all_l1)  # Extracts corrected biological ages for all L1 samples.
ages_wt_all = get_group_ages(df, wt_all)  # Extracts corrected biological ages for all WT samples.
ages_xol1_all = get_group_ages(df, xol1_all)  # Extracts corrected biological ages for all xol-1 samples.

plot_groups = {  # Creates a dictionary where each graph label is connected to its age values.
    "WT Embryo": ages_wt_embryo,  # Stores WT embryo ages.
    "WT L1": ages_wt_l1,  # Stores WT L1 ages.
    "All Embryos": ages_all_embryo,  # Stores all embryo ages.
    "All L1 Larvae": ages_all_l1,  # Stores all L1 ages.
    "xol-1 Embryo": ages_xol1_embryo,  # Stores xol-1 embryo ages.
    "xol-1 L1": ages_xol1_l1,  # Stores xol-1 L1 ages.
    "WT All": ages_wt_all,  # Stores all WT ages.
    "xol-1 All": ages_xol1_all  # Stores all xol-1 ages.
}  # Ends the plot group dictionary.

group_names = list(plot_groups.keys())  # Creates a list of group names for the x-axis labels.
means = [plot_groups[g].mean() for g in group_names]  # Calculates the mean corrected biological age for each group.
errors = [plot_groups[g].sem() for g in group_names]  # Calculates SEM for each group to use as error bars.

plt.figure(figsize=(13, 6))  # Creates a new figure and sets its size.

x = range(len(group_names))  # Creates x-axis positions for each group.
plt.bar(x, means, yerr=errors, capsize=5)  # Creates the bar plot with SEM error bars

plt.ylabel("Corrected Biological Age")  # Labels the y-axis.
plt.xlabel("Group")  # Labels the x-axis.
plt.title("Corrected Biological Age Across Statistical Comparison Groups")  # Adds the plot title.
plt.xticks(x, group_names, rotation=45, ha="right")  # Adds group names to the x-axis and rotates them for readability.

def add_sig_bar(x1, x2, y, text):  # Defines a helper function to draw significance brackets.
    plt.plot([x1, x1, x2, x2], [y, y + 2, y + 2, y], color="black")  # Draws the bracket line between two bars.
    plt.text((x1 + x2) / 2, y + 2.5, text, ha="center", fontsize=14)  # Places the asterisk text above the bracket.

comparisons = [  # Defines the exact group comparisons used for significance annotations.
    ("WT L1", "WT Embryo", ages_wt_l1, ages_wt_embryo),  # Compares WT L1 against WT embryo.
    ("All L1 Larvae", "All Embryos", ages_all_l1, ages_all_embryo),  # Compares all L1 samples against all embryo samples.
    ("xol-1 Embryo", "WT Embryo", ages_xol1_embryo, ages_wt_embryo),  # Compares xol-1 embryo against WT embryo.
    ("xol-1 L1", "WT L1", ages_xol1_l1, ages_wt_l1),  # Compares xol-1 L1 against WT L1.
    ("xol-1 All", "WT All", ages_xol1_all, ages_wt_all)  # Compares all xol-1 samples against all WT samples.
]  # Ends the comparisons list

max_y = max([m + e for m, e in zip(means, errors)])  # Finds the highest bar plus its error bar to determine where to start placing significance brackets
current_y = max_y + 5  # Sets the first significance bracket slightly above the tallest bar

for label1, label2, group1, group2 in comparisons:  # Loops through each comparison.
    t_stat, p_val = ttest_ind(group1, group2, equal_var=False)  # Runs Welch's t-test for graph annotation 

    if p_val < 0.05:  # Checks if the comparison is statistically significant.
        x1 = group_names.index(label1)  # Finds the x-position of the first group
        x2 = group_names.index(label2)  # Finds the x-position of the second group

        if p_val < 0.001:  # Checks if the p-value is below 0.001
            stars = "***"  # Uses three stars for p < 0.001.
        elif p_val < 0.01:  # Checks if the p-value is below 0.01
            stars = "**"  # Uses two stars for p < 0.01.
        else:  # Runs if p-value is below 0.05 but not below 0.01
            stars = "*"  # Uses one star for p < 0.05.

        add_sig_bar(x1, x2, current_y, stars)  # Adds the significance bracket and asterisk to the plot.
        current_y += 10  # Moves the next bracket higher so brackets do not overlap

plt.tight_layout()  # Adjusts spacing so labels and title fit cleanly.
plt.savefig(PLOT_FILE, dpi=300)  # Saves the plot as a  PNG file
plt.show()  # Displays the plot if the system supports graphical display.

print(f"Plot saved to: {PLOT_FILE}")  #prints the saved plot location in the terminal.
