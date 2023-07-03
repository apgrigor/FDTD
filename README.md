## FDTD simulations for Electromagnetism

This repository contains FDTD simulation scripts in Python and C for the graduate EM course at the National and Kapodistrian University of Athens.

### Installation (Linux)

You need to download the repo, the required python libraries and compile the C extention module.

```sh
git clone https://github.com/ShadyStokes/FDTD
cd FDTD
```

For the 3D simulation:

```sh
cd 3d
make
```
This will compile the module ``FDTD/3d/fd3d.so`` which is then imported by ``FDTD/3d/*.py``. You may then run 

```sh
python fd3d_0.py
```

Make sure the libraries in ``FDTD/requirements.txt`` are installed. This can be done via:

```sh
pip install -r requirements.txt
```

(note that this command will only work as-is if you're in the ``FDTD`` directory, otherwise you'll have to change the path to ``requirements.txt``)

The output will be in ``FDTD/3d/img``