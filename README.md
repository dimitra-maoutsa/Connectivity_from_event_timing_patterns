# Inferring network connectivity from event timing patterns

Implementations of network simulations and reconstructions as described in:


**Casadiego*, Maoutsa*, Timme, _Inferring network connectivity from event timing patterns_, Physical Review Letters 2018**  

For further information refer to the article.


<br>
### Running the code:
1. Generate input data
    - Either extract provided data
    
        ```
        cd Connectivity_from_event_timing_patterns/simulate_network
        tar -xzvf Data.tar.gz
        ```
    - Or simulate network 
    
        ```
        python Connectivity_from_event_timing_patterns/simulate_network/simulate_network.py
        ```
2. Reconstruct
    ```
    python Connectivity_from_event_timing_patterns/reconstruct_network/inferring_connections_from_spikes.py
    
    ```


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