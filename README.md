# TweakAESKS

This repo contains the implementations used for the paper "Tweaking the AES Key-Schedule for Better Truncated Differential Bounds"

Required :
- [Minizinc](http://www.minizinc.org/index.html)
- [Gecode](http://www.gecode.org/)
- [Choco](http://www.choco-solver.org/)

The folder organization is the following :
- `search_cycle` contains the code used to prove Theorem 2 in Section 4.1
- `meta_search` contains the code for Algorithm 1 (Tweaked Simulated Annealing) in Section 4.2
- `searchAllPath`contains the Minizinc model to search for all truncated paths of a given length (end of Section 4.2)
- `search_instanciate`contains the Choco model to search for an instantiation with probability > 2^{-128}

