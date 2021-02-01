# Inferring network connectivity from event timing patterns



Model-free method for inferring synaptic interactions from spike train recordings.

By mapping spike timing data in **event spaces** _(spanned by inter-spike and cross-spike intervals)_,  
we identify synaptic interactions in networks of spiking neurons through **Event Space Linearisations (ESL)** .  
Here, we provide implementations of network simulations and reconstructions as described in:  
**Casadiego\*, Jose, Maoutsa\*, Dimitra, Timme, Marc,
 _Inferring network connectivity from event timing patterns_, Physical Review Letters 2018**  
For further information refer to the [article](https://github.com/dimitra-maoutsa/Connectivity_from_event_timing_patterns/blob/di.ma-master-patch-52805/PhysRevLett.121.054101.pdf) and the [supplementary info](https://github.com/dimitra-maoutsa/Connectivity_from_event_timing_patterns/blob/di.ma-master-patch-52805/Supplementary_Inferring_network_connectivity_from_event_timing_patterns.pdf). (can be found [here](https://github.com/dimitra-maoutsa/Connectivity_from_event_timing_patterns/blob/di.ma-master-patch-52805/PhysRevLett.121.054101.pdf) and [here](https://github.com/dimitra-maoutsa/Connectivity_from_event_timing_patterns/blob/di.ma-master-patch-52805/Supplementary_Inferring_network_connectivity_from_event_timing_patterns.pdf) as pdf) .


<br>

### Running the code:
1. Generate input data
    - Either extract provided data
    
        ```
        cd Connectivity_from_event_timing_patterns/simulate_network
        tar -xzvf Data.tar.gz
        ```
    - Or simulate network (requires [NEST simulator] (http://www.nest-simulator.org/) )
    
        ```
        python Connectivity_from_event_timing_patterns/simulate_network/simulate_network.py
        ```
2. Reconstruct
    ```bash
    python Connectivity_from_event_timing_patterns/reconstruct_network/inferring_connections_from_spikes.py
    
    ```
    **Caution:** Input data files should be stored in folder `simulate_network/Data/`



<br>

### Support:
For questions please contact: Dimitra Maoutsa [ dimitra.maoutsa <-at-> tu-berlin.de ] 

### Cite:
```
@article{ESL18,
  title = {Inferring Network Connectivity from Event Timing Patterns},
  author = {Casadiego, Jose and Maoutsa, Dimitra and Timme, Marc},
  journal = {Phys. Rev. Lett.},
  volume = {121},
  issue = {5},
  pages = {054101},
  numpages = {6},
  year = {2018},
  month = {Aug},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevLett.121.054101},
  url = {https://link.aps.org/doi/10.1103/PhysRevLett.121.054101}
}

```
