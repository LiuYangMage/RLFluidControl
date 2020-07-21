# RLFluidControl
Reinforcement learning for active fluid control.

See our paper: Reinforcement Learning for Active Flow Control in Experiments (https://arxiv.org/abs/2003.03419)

## Dependencies

Python Packages: TensorFlow 1.x (1.13+ suggested), numpy

## How to run

Start the server (RL agent) by

```shell
$ cd server
$ python server.py -env CFD -fil None
```

`-fil None` means that we don't use any filter, in consistent to our practice for the CFD environment in the paper.

Start the client:
### 2D BDIM:

We provide a simple example using Lilypad, a code based on a 2D BDIM method developed by Dr. Gab Weymouth. See https://github.com/weymouth/lily-pad.

Run **`clientLilypad/clientCFD.pde`**  in Processing (download from https://processing.org/download/). 

### 3D LES:
See details in **`clientNektar/Readme`**
