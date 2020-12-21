# PhylUp - updating of phylogenetic alignments with custom sampling
![](./workflow_Dez2020.png)   
    
## Short introduction to the program:

PhylUp is a command-line program written in python3 to automatically update alignments and phylogenies with a focus on different sampling strategies.
As input it needs a alignment or a single sequence (and if available a phylogeny) and 
a file with the information about the sequence names and the corresponding species names. 
PhylUp will take every input sequence and blasts it against the ncbi GenBank database. 
Sequences that are similar to the input sequence will be added to the alignment, 
if they are a different taxon and/or they are longer than existing sequences.
They are then filtered to user settings provided in the configuration file.
Newly found and filtered sequences will be blasted again until no new sequences were found.
Finally, it will place the newly found sequences into the alignment and if enabled to update phylogenies.

After the single-gene datasets are updated, the data can be concatenated. 
The tool decides randomly which sequences to combine if there are more than a single sequence for a taxon in one of the alignments.


## Tutorial:

To get started please view the Wiki for more details.

## Examples:

To update an alignment I provide different example files in `examples`:
  * The file `example_analysis.py` updates a small subset of an alignment provided in `data`.
 * Filtering with different number of sequences for different taxonomic ranks is shown in ` example_different_levels_wrapper.py`.
 * How to concatenate different alignments and calculate the phylogeny is shown in `example_concat.py`.
   




  
