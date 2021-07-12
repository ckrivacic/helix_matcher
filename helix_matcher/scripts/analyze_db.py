import sys
import pandas as pd
import matplotlib.pyplot as plt


def main():
    dfpath = 'nr_dataframes/final.pkl'
    df = pd.read_pickle(dfpath)
    df.hist(column='length', bins=100)
    df = df[df[show] > 400]
    plt.show()


if __name__=="__main__":
    show = sys.argv[1]
    main()
