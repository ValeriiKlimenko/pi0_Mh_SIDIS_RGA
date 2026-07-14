import matplotlib.pyplot as plt

iterations = [1, 2, 3, 4, 5, 6, 7, 8, 9]
chi2_change = [3.13172e+08, 1.06791e+07, 3.62924e+06, 1.87379e+06,
               1.39052e+06, 1.34152e+06, 1.50438e+06, 1.74404e+06, 1.89136e+06]

plt.figure(figsize=(7, 4.5))
plt.plot(iterations, chi2_change, marker="o", linewidth=2)
plt.xlabel("Iteration")
plt.ylabel(r"$\Delta \chi^2$")
plt.title(r"RooUnfold: $\Delta \chi^2$ vs Iteration")
plt.yscale("log")
plt.grid(True, which="both", linestyle="--", alpha=0.4)
plt.xticks(iterations)
plt.tight_layout()

plt.savefig("roounfold_chi2_change_vs_iteration.png", dpi=300, bbox_inches="tight")
