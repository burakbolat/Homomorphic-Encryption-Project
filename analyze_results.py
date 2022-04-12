import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def main():
    results = pd.read_csv('results.csv')
    group = results.groupby('NodeCount')

    x_axis = results['NodeCount'].unique()
    fig, ax = plt.subplots(1)
    for metric in results.columns[2:-1]:
        ax.plot(x_axis, group[metric].mean(), lw=2, label=metric)
        ax.fill_between(x_axis, group[metric].max(), group[metric].min(), alpha=0.4)
        
    ax.legend(loc='best')
    ax.set_xlabel('num of nodes')
    ax.set_xticks(x_axis)
    ax.set_ylabel('time (ms)')
    ax.set_title("Time analysis with respect to num of nodes")
    ax.grid()
    plt.savefig("results/first.png")
    
    fig.clf()
    fig, ax = plt.subplots(1)
    ax.plot(x_axis, group["Mse"].mean(), lw=2, label=metric)
    ax.fill_between(x_axis, group["Mse"].max(), group["Mse"].min(), alpha=0.4)
    ax.set_xlabel('num of nodes')
    ax.set_xticks(x_axis)
    ax.set_ylabel('Mse')
    ax.set_title("Mse with respect to num of nodes")
    ax.grid()
    plt.savefig("results/second.png")
    
if __name__=="__main__": main()