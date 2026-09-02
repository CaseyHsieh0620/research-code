import numpy as np
import matplotlib.pyplot as plt
rng = np.random.default_rng()
import math

#pauli Matrices
sigmaZ = np.array([[1.,0.], [0.,-1.]])
sigmaX = np.array([[0.,1.], [1.,0.]])
I = np.eye(2)
# this describes coupling THi si sthe J
K = [.2, .5, 1]

def generateHamiltonian(K, n):
    list = []
    out = []
    couplings = []
    # build the randomized array
    for i in range(n):
        # this describes the actual Hamiltonians. This is h
        Jn = rng.uniform(1,1)
        list.append(Jn)
    # build hamiltonian terms individually 

    for i in range(n):
        temp = [[1.]]
        for j in range(n):
            if (i == j):
                temp = np.kron(temp, sigmaX)
            else:
                temp = np.kron(temp, I)
        out.append(temp)
    # build coupling terms

    for i in range(n):
        for j in range(i + 1, n):
            temp = [[1.]]
            for k in range(n):
                if (k == j):
                    temp = np.kron(temp, sigmaZ)
                elif (k == i):
                    temp = np.kron(temp, sigmaZ)
                else:
                    temp = np.kron(temp, I)
            couplings.append(temp)
    H = np.zeros((2**n, 2**n))
    for j in range(n):
        H = H + list[j] * out[j]
    for coupling in couplings:
        H = H + K * coupling
    return H
                

def generate(n, number):
    for k, coupling in enumerate(K):
        x = 0
        while(x < number):
            final = generateHamiltonian(coupling, n)
            eigenvalues = np.linalg.eigvalsh(final)
            s = np.diff(eigenvalues)
            #s = s/s.mean()
            for a in s:
                position = int(a / .1)
                if position >= 40 or position < 0:
                    continue
                counts[k][position] = counts[k][position] + 1
            x = x + 1
    probs = counts / counts.sum(axis=1, keepdims=True)

    ncols = math.ceil(math.sqrt(len(K)))
    nrows = math.ceil(len(K) / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5*ncols, 3*nrows),
                             sharex=True, sharey=True, squeeze=False)
    axes = axes.ravel()

    centers = np.arange(40) * 0.1 + 0.05
    x = np.linspace(0.01, 4, 200)

    for i in range(len(K)):
        ax = axes[i]
        ax.bar(centers, probs[i] / 0.1, width=0.1, alpha=0.8)
        ax.plot(x, np.exp(-x), 'k--', lw=1.5)
        ax.plot(x, (np.pi/2)*x*np.exp(-np.pi*x**2/4), 'r--', lw=1.5)
        ax.set_title(f'J = {K[i]:g}')

    for j in range(len(K), len(axes)):
        axes[j].axis('off')          # hide unused panels

    fig.supxlabel('spacing'); fig.supylabel('P(s)')
    plt.tight_layout()
    plt.show()



while True:
    inp = input("Number of spins (q to quit):")
    if inp.lower() in ('q'):
        break
    try:
        n = int(inp)
    except ValueError:
        print("Print a Number or Q")

    trials = input("Input the amoutn of randomized Hamiltonains:")
    if trials.lower() in ('q'):
            break
    try:
        trial = int(trials)
    except ValueError:
        print("Print a Number or Q")



    counts = np.zeros((len(K), 40))
    generate(n, trial)


