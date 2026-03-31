import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_pkfg_results():
    # ログファイルの読み込み
    py_log = "../Logs/worm_log_python.txt"
    ft_log = "../Logs/worm_log_fortran.txt"
    
    if not (os.path.exists(py_log) and os.path.exists(ft_log)):
        print("Logs not found. Please run the simulations first.")
        return

    df_py = pd.read_csv(py_log)
    df_ft = pd.read_csv(ft_log)
    
    # 共通のステップ数
    steps = min(len(df_py), len(df_ft))
    df_py = df_py.iloc[:steps]
    df_ft = df_ft.iloc[:steps]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)

    # 1. 餌への距離 (1.0, 1.0)
    target = np.array([1.0, 1.0])
    dist_py = np.sqrt((df_py['HeadX'] - target[0])**2 + (df_py['HeadY'] - target[1])**2)
    dist_ft = np.sqrt((df_ft['HeadX'] - target[0])**2 + (df_ft['HeadY'] - target[1])**2)
    
    ax1.plot(df_py['Step'], dist_py, label='Python (Distance to Food)', color='blue', alpha=0.7)
    ax1.plot(df_ft['Step'], dist_ft, label='Fortran (Distance to Food)', color='red', linestyle='--', alpha=0.7)
    ax1.set_ylabel("Distance to Food (1.0, 1.0)")
    ax1.set_title("Worm Navigation: Approaching and Passing the Food")
    ax1.legend()
    ax1.grid(True)

    # 2. 蛇行 (HeadX の変動)
    ax2.plot(df_py['Step'], df_py['HeadX'], label='Python HeadX (Undulation)', color='cyan')
    ax2.plot(df_ft['Step'], df_ft['HeadX'], label='Fortran HeadX (Undulation)', color='orange', linestyle=':')
    ax2.set_ylabel("Head X Position")
    ax2.set_xlabel("Steps")
    ax2.set_title("Undulation Pattern: Sinusoidal Locomotion Emergence")
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig("../Logs/worm_movement_analysis.png")
    print("Graph saved to Logs/worm_movement_analysis.png")

if __name__ == "__main__":
    plot_pkfg_results()
