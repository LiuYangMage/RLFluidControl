# RLFluidControl
Reinforcement learning for active fluid control.

See our paper: Reinforcement Learning for Active Flow Control in Experiments (https://arxiv.org/abs/2003.03419)

## Dependencies
Processing from https://processing.org/download/

Python Packages: TensorFlow 1.x (1.13+ suggested), numpy

## How to run

Start the server (RL agent) by

```shell
$ cd server
$ python server.py -env CFD -fil None
```

`-fil None` means that we don't use any filter, in consistent to our practice for the CFD environment in the paper.

Start the client (CFD environment):

Run **`clientCFD/clientCFD.pde`**  in Processing. 

## External Links

The CFD processing code is based on a 2D BDIM code ''lilypad'' developed by Dr. Gab Weymouth. Here is the Github link: https://github.com/weymouth/lily-pad 