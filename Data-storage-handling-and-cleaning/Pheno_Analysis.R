## Preliminary analyses (checking for signficant differences between males and females)

# Define the data for Eye Color frequencies from the user's table
data_eye_color = {
    "Blue": [1986, 2157],
    "Brown": [1315, 1575],
    "Green": [904, 1547],
    "Hazel": [926, 1676]
}

# Create a DataFrame from the data
df_eye_color = pd.DataFrame(data_eye_color, index=["Males", "Females"]).T
df_eye_color['Total'] = df_eye_color['Males'] + df_eye_color['Females']

# Calculate the chi-square test of independence
chi2_stat_eye, p_value_eye, dof_eye, expected_eye = chi2_contingency(df_eye_color[['Males', 'Females']])

# Convert the expected frequencies into a DataFrame for easy comparison
df_expected_eye = pd.DataFrame(expected_eye, index=df_eye_color.index, columns=['Expected Males', 'Expected Females'])

# Calculate odds ratios and 95% confidence intervals for eye color
odds_ratios_eye = {}
for color in data_eye_color.keys():
    count_eye = [df_eye_color.loc[color, 'Males'], df_eye_color.loc[color, 'Females']]
    nobs_eye = [df_eye_color['Males'].sum(), df_eye_color['Females'].sum()]
    stat_eye, pval_eye = proportions_ztest(count_eye, nobs_eye)
    odds_ratio_eye = (count_eye[0] / (nobs_eye[0] - count_eye[0])) / (count_eye[1] / (nobs_eye[1] - count_eye[1]))
    # CI for odds ratio
    se_eye = ((1/count_eye[0]) + (1/(nobs_eye[0]-count_eye[0])) + (1/count_eye[1]) + (1/(nobs_eye[1]-count_eye[1])))**0.5
    ci_low_eye = odds_ratio_eye / np.exp(1.96*se_eye)
    ci_high_eye = odds_ratio_eye * np.exp(1.96*se_eye)
    odds_ratios_eye[color] = {'Odds Ratio': odds_ratio_eye, '95% CI': (ci_low_eye, ci_high_eye), 'p-value': pval_eye}

# Present the results
(df_eye_color, df_expected_eye, odds_ratios_eye, chi2_stat_eye, p_value_eye, dof_eye)

## Repeat for other phenotypes

## For genotype-stratefied analysis:

import pandas as pd
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt

# Input data
data = {
    "Eye Colour": ["Blue", "Green", "Hazel", "Brown"],
    "GG (Blue background)": [2135, 860, 256, 231],
    "AA+AG (Non-blue background)": [411, 960, 1537, 1804],
}
df = pd.DataFrame(data)

# Total counts for GG and AA/AG
total_gg = sum(df["GG (Blue background)"])
total_aa_ag = sum(df["AA+AG (Non-blue background)"])

# Expected proportions based on 23andMe data
gg_expected_proportions = {
    "Blue": 0.72,  # Blue + greenish-blue
    "Green": 0.17,
    "Hazel": 0.11,
    "Brown": 0.01,
}

aa_ag_expected_proportions = {
    "Blue": 0.02,
    "Green": 0.11,
    "Hazel": 0.42,
    "Brown": 0.45,
}

# Calculate expected counts for GG and AA/AG
gg_expected = {k: gg_expected_proportions[k] * total_gg for k in gg_expected_proportions.keys()}
aa_ag_expected = {k: aa_ag_expected_proportions[k] * total_aa_ag for k in aa_ag_expected_proportions.keys()}

# Combine expected counts
expected_total = {
    eye_color: gg_expected.get(eye_color, 0) + aa_ag_expected.get(eye_color, 0)
    for eye_color in ["Blue", "Green", "Hazel", "Brown"]
}

# Combine observed counts
observed_total = {
    "Blue": df.loc[df["Eye Colour"] == "Blue", ["GG (Blue background)", "AA+AG (Non-blue background)"]].sum(axis=1).values[0],
    "Green": df.loc[df["Eye Colour"] == "Green", ["GG (Blue background)", "AA+AG (Non-blue background)"]].sum(axis=1).values[0],
    "Hazel": df.loc[df["Eye Colour"] == "Hazel", ["GG (Blue background)", "AA+AG (Non-blue background)"]].sum(axis=1).values[0],
    "Brown": df.loc[df["Eye Colour"] == "Brown", ["GG (Blue background)", "AA+AG (Non-blue background)"]].sum(axis=1).values[0],
}

# Perform goodness-of-fit test
observed_values = list(observed_total.values())
expected_values = list(expected_total.values())
chi2, p, dof, expected_table = chi2_contingency([observed_values, expected_values])

# Summarize results
goodness_of_fit_results = pd.DataFrame({
    "Eye Colour": expected_total.keys(),
    "Observed": observed_values,
    "Expected": expected_values,
})

# Visualize the observed vs expected frequencies for GG and AA/AG
expected_gg = [gg_expected[ec] for ec in ["Blue", "Green", "Hazel", "Brown"]]
expected_aa_ag = [aa_ag_expected[ec] for ec in ["Blue", "Green", "Hazel", "Brown"]]
observed_gg = df["GG (Blue background)"].tolist()
observed_aa_ag = df["AA+AG (Non-blue background)"].tolist()

# Plotting observed vs expected, broken down by genetic background
fig, ax = plt.subplots(figsize=(10, 6))
bar_width = 0.2
x = range(len(["Blue", "Green", "Hazel", "Brown"]))

# Observed bars
ax.bar(x, observed_gg, bar_width, label="Observed GG", align="center", alpha=0.7)
ax.bar([i + bar_width for i in x], observed_aa_ag, bar_width, label="Observed AA/AG", align="center", alpha=0.7)

# Expected bars
ax.bar([i + 2 * bar_width for i in x], expected_gg, bar_width, label="Expected GG", align="center", alpha=0.7)
ax.bar([i + 3 * bar_width for i in x], expected_aa_ag, bar_width, label="Expected AA/AG", align="center", alpha=0.7)

# Adding labels, title, and legend
ax.set_xticks([i + 1.5 * bar_width for i in x])
ax.set_xticklabels(["Blue", "Green", "Hazel", "Brown"])
ax.set_xlabel("Eye Colour")
ax.set_ylabel("Number of Individuals")
ax.set_title("Observed vs. Expected Frequencies by Genetic Background")
ax.legend()

# Display the plot
plt.tight_layout()
plt.show()

# Print results
print(f"Chi-squared statistic: {chi2}")
print(f"P-value: {p}")
print(f"Degrees of Freedom: {dof}")

# Display results table
print(goodness_of_fit_results)
