# Inferring network connectivity from event timing patterns




<p style="text-align:justify;">By mapping spike timing data in **event spaces** _(spanned by inter-spike and cross-spike intervals)_,  
we identify synaptic interactions in networks of spiking neurons through **Event Space Linearisations (ELS)** without assuming any particular neuronal dynamical model. 

Here, we provide implementations of network simulations and reconstructions as described in:  
**Casadiego*, Maoutsa*, Timme, _Inferring network connectivity from event timing patterns_, Physical Review Letters 2018**  
For further information refer to the article.
</p>

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
@article{ELS18,
    Author = {Jose Casadiego and Dimitra Maoutsa and Marc Timme},
    Title = {Inferring network connectivity from event timing patterns},
    Year = {2018},
    Journal = {Physical Review Letters},
    Volume={},
    Number={},
    Pages={},
    Publisher={APS}
}
```