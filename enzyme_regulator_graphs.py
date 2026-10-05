import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_rel, wilcoxon, spearmanr
from numpy.polynomial.polynomial import Polynomial
import numpy as np
from sklearn.metrics import r2_score

# Dosyayı oku
df = pd.read_csv("Enzyme_Regulator_summary.csv")



# Stil
sns.set(style="whitegrid")

# Veri yükle
df['label'] = df['filename'].str.split('_').str[2]
df['diff'] = df['main_diversity'] - df['reg_diversity']

# Genel istatistiksel testler
ttest_p = ttest_rel(df['main_diversity'], df['reg_diversity']).pvalue
wilcoxon_p = wilcoxon(df['main_diversity'], df['reg_diversity']).pvalue

print(f"Paired t-test p-value: {ttest_p:.4g}")
print(f"Wilcoxon test p-value: {wilcoxon_p:.4g}")

# 1. Scatter plot (main vs reg)
plt.figure(figsize=(10, 6))
sns.scatterplot(data=df, x='main_diversity', y='reg_diversity', s=100)

x = df['main_diversity']
y = df['reg_diversity']
coeffs = Polynomial.fit(x, y, deg=1).convert().coef  # Get actual coefficients

# Tahmin çizgisi
x_fit = np.linspace(x.min(), x.max(), 200)
y_fit = coeffs[0] + coeffs[1] * x_fit  # Linear: y = a + bx

# Plot çiz
plt.plot(x_fit, y_fit, color='black', linestyle='--', label='Linear Fit')
# R² hesapla
y_pred = coeffs[0] + coeffs[1] * x
r2 = r2_score(y, y_pred)

# Sapmalar (residuals)
residuals = np.abs(y - y_pred)

# En çok sapan 10 noktanın indeksleri
top10_idx = residuals.nlargest(5).index

# Bu noktalara label ekle
for idx in top10_idx:
    plt.text(x[idx], y[idx], str(df.loc[idx, 'filename']).split('_')[2], fontsize=9, color='black')

# Opsiyonel: R² yaz
plt.text(0.05, 0.95, f"$R^2$ = {r2:.2f}", transform=plt.gca().transAxes, fontsize=12)


plt.plot([0, 1], [0, 1], 'k--', alpha=0.5)
plt.xlabel("Main Gene Diversity")
plt.ylabel("Regulator Gene Diversity")
plt.title(f"Main vs Reg Diversity\n(paired t-test p = {ttest_p:.4g}, wilcoxon p = {wilcoxon_p:.4g}, R^2:{r2})")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("scatter_main_vs_reg.pdf")
plt.close()

# 2. Boxplot (main vs reg by label)
df_melted = df.melt(id_vars='label', value_vars=['main_diversity', 'reg_diversity'],
                    var_name='Gene Type', value_name='Diversity')

plt.figure(figsize=(12, 6))
sns.boxplot(data=df_melted, x='label', y='Diversity', hue='Gene Type')
plt.title("Diversity Distribution by Label and Gene Type")
plt.ylabel("Diversity (p-distance)")
plt.xlabel("Gene Label")
plt.tight_layout()
plt.savefig("boxplot_diversity_by_label.pdf")
plt.close()

# 3. Violin plot of diversity difference
plt.figure(figsize=(10, 6))
sns.violinplot(data=df, x='label', y='diff', inner='point', linewidth=1)
plt.axhline(0, linestyle='--', color='gray')
plt.ylabel("Main - Reg Diversity")
plt.title("Difference in Diversity (Main - Reg) by Label")
plt.tight_layout()
plt.savefig("violin_diff_by_label.pdf")
plt.close()

# 4. Pair count vs diversity diff
spearman_corr, spearman_p = spearmanr(df['pair_count'], df['diff'])

plt.figure(figsize=(10, 6))
sns.scatterplot(data=df, x='pair_count', y='diff', hue='label', s=100)
plt.axhline(0, color='gray', linestyle='--')
plt.xlabel("Pair Count")
plt.ylabel("Main - Reg Diversity")
plt.title(f"Pair Count vs Diversity Difference\n(Spearman ρ = {spearman_corr:.2f}, p = {spearman_p:.4g})")
plt.tight_layout()
plt.savefig("scatter_paircount_vs_diff.pdf")
plt.close()

# 5. Her bir label için Wilcoxon testi
print("\nLabel bazında Wilcoxon testleri:")
for label in sorted(df['label'].unique()):
    sub = df[df['label'] == label]
    if len(sub) >= 2:  # En az 2 örnek gerek
        try:
            p_val = wilcoxon(sub['main_diversity'], sub['reg_diversity']).pvalue
            print(f"{label}: Wilcoxon p = {p_val:.4g} (n = {len(sub)})")
        except:
            print(f"{label}: Test yapılamadı (veri yetersiz)")
    else:
        print(f"{label}: Veri yetersiz (n = {len(sub)})")

print("\nGrafikler PDF olarak kaydedildi.")
