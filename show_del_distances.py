from cdg import DataSet
import matplotlib.pyplot as plt
import numpy as np


def main():
    dataset = DataSet.load_dataset('./datasets/delaunay', precomputed_distance=True)
    print(dataset.prec_distance_mat)

    cl = list(dataset.elements.keys())[1::2] + [0]
    n_cl = len(cl)

    xl = [dataset.prec_distance_mat.min(), dataset.prec_distance_mat.max()]

    import seaborn as sns; #sns.set(color_codes=True)
    sns.set_palette("Set2", 8)
    
    fig = plt.figure(figsize=(5, 5))
    
    for i, c in enumerate(cl[:-1]):
        plt.subplot(n_cl-1, 1, i+1)
        # plt.plot()
        ut0 = np.triu_indices(len(dataset.elements[0]), 1)
        ut = np.triu_indices(len(dataset.elements[c]), 1)
        ax = sns.kdeplot(dataset.distance_measure(dataset.elements[0], dataset.elements[0])[ut0].ravel(), label="d(0, 0)", shade=True)
        ax = sns.kdeplot(dataset.distance_measure(dataset.elements[c], dataset.elements[c])[ut].ravel(),  label="d(c, c)", dashes=(6, 3), shade=True)
        ax = sns.kdeplot(dataset.distance_measure(dataset.elements[0], dataset.elements[c]).ravel(),      label="d(c, 0)", dashes=(3, 1.5), shade=True)
    
        plt.xlim(xl)
        plt.ylabel("Cl. 0 vs {}".format(c))
    plt.legend()
    plt.xlabel('Graph distance')
    plt.tight_layout(pad=0)
    fig.savefig('delaunay_distances.pdf')

if __name__ == "__main__":
    main()