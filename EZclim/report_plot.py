import matplotlib.pyplot as plt
import seaborn as sns

# Serial
# Ens  | Time   | Parallel
# 1      | 2.786   |  3.008
# 5      |  10.38   | 8.139
# 10    |  18.086 | 12.567
# 15    | 25.37   |  17.45
# 20    | 37.313 | 21.345
# 25    |  49.255 | 23.789
# 30    |  54.265 | 25.879
# 35    | 65.035 | 26.80
# 40    | 74.544 |27.913

ens = [1, 5, 10, 15, 20, 25, 30, 35, 40]
s = [2.786, 10.38, 18.086, 25.37, 37.313, 49.255, 54.265, 65.035, 74.544]
p = [3.008, 8.139, 12.567, 15.45, 17.46, 18.34, 18.881, 18.77, 18.76]
sns.set(rc={'figure.figsize': (11, 4)})
fig, ax = plt.subplots()
sns.lineplot(ens, s, ax=ax, label="Serial")
sns.lineplot(ens, p, color='g', ax=ax, label="Parallel")
ax.set_title("Measuring time taken to run a range of ensembles using serialised and parallelised version of code.")
ax.set_xlabel("Number of ensembles")
ax.set_ylabel("Time taken to run (seconds)")
plt.legend()
plt.show()