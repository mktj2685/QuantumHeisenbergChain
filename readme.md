# Quantum Heisenberg Chain

## Description

Exact diagonalization of the following types of Hamiltonians :

$$
    \hat{H} = \sum_{i=1}^N (\hat{S}^x_i\hat{S}^x_{i+1} + \hat{S}^y_i\hat{S}^y_{i+1}+\hat{S}^z_{i}\hat{S}^z_{i+1})
$$

where $N$ is number of spins.

## Install

```
git clone https://github.com/mktj2685/QuantumHeisenbergChain.git
cd QuantumHeisenbergChain
pip install -r requirements.txt
```

## Usage

### Case : S=1/2
```
python heisenberg.py --S 0.5 --N 10
```

### Case : S=1
```
python heisenberg.py --S 1.0 --N 10
```

### Case : S=3/2
```
python heisenberg.py --S 1.5 --N 10
```