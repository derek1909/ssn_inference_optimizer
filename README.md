# Training of the stochastic SSN

This is a basic engine to train an SSN to exhibit desired patterns of first and
second-order moments as presented in:

Echeveste, R., Aitchison, L., Hennequin, G., & Lengyel, M. (2019). Cortical-like dynamics in recurrent circuits optimized for sampling-based probabilistic inference. bioRxiv, 696088.

To source the Ocaml environment (you'll have to do that everytime you want to
compile an OCaml program):

```sh
eval `opam config env`
```

To compile the executable:

```sh
ocamlbuild -use-ocamlfind train.native
```

To use the executable:

```sh
./train.native -d result_dir -m 50 -n_targ 10
```

There are other command-line parameters with good default values; you can check these
out directly in `train.ml` (look for code of the form `Cmdargs.(...)`).

You must first save all your desired moments ("targets") in `result_dir`.
Format should be plain text files, called `h0, h1, ...`, `mu0, mu1, ...`,
`sigma0, sigma1, ...`. Optimization progress will be displayed online, and
results will be saved every 100 iterations.

The input baseline is learnt.
The saved "h_true" include that baseline. Target patterns need to already
include their own baseline and therefore, the saved mu_aux include it too.
 

# In case of segfault

1. Talk to your operator
2. You can reuse the state vector from a previous round of optimization, using
```sh
./train.native [...] -reuse /path/to/state_vector.bin
```

# Reading parameters from the state_vector.bin file

The code `read_out_state_vector.ml` allows one to read out all relevant parameters 
out of the binary state vector file, and can compute the corresponding moments.

# Using the 'train_mult_init.sh' script

The script runs the train optimizer several times starting from random initial conditions.
As it is, the targets need to be placed in a subdir called `targets` and the results will be 
stored in a subdir `results`.

# Terms of use

This repository is provided for reproducibility purposes only. Please contact the authors if you wish to use any parts of this code for other purposes.

