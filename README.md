# Inferring network connectivity from event timing patterns

Implementations of network simulations and reconstructions as described in:


**Casadiego*, Maoutsa*, Timme, _Inferring network connectivity from event timing patterns_, Physical Review Letters 2018**  

For further information refer to the article.


<br>
### Running the code:
1. Generate input data
    - Either extract provided data
    
        ```
        tar -xzvf simulate_network/Data.tar.gz
        ```
    - Or simulate network 
    
        ```
        python simulate_network.py
        ```
2. Reconstruct
    ```
    python reconstruct_network/inferring_connections_from_spikes.py
    
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