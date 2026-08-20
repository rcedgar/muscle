![Muscle5](http://drive5.com/images/muscle5_header.jpg)

Muscle is widely-used software for making multiple alignments of biological sequences. 

Muscle achieves highest scores on Balibase, Bralibase and Balifam benchmark tests and scales to thousands of sequences or structures on a commodity desktop computer.

Muscle supports generating an ensemble of alternative alignments with the same high accuracy obtained with default parameters. By comparing downstream predictions from different alignments, such as trees, a biologist can evaluation the robustness of conclusions against alignment variation caused by ambiguities and errors.

### Multiple structure alignment

Structure alignment ("Muscle-3D") is supported as well as conventional amino acid sequence alignment. Muscle accepts the same **STRUCTS** inputs as [reseek](https://github.com/rcedgar/reseek): `.pdb` / `.cif`, `.cal`, `.bca` / `.bcb`, `.files` lists, or a directory. Feature alphabets default to reseek `-stats sf`.

Legacy text `.mega` files (from older `reseek -pdb2mega`) are still accepted but **deprecated**.


[<img src="https://drive5.com/reseek/youtube_snip_muscle3d.gif" width="150">](https://www.youtube.com/watch?v=BzIgqdm9xDs)

<pre>
# for up to ~100 structures
muscle -align STRUCTS -output structs.afa

# for up to ~10,000 structures (guide tree / distance matrix from reseek)
reseek -convert STRUCTS -bca structs.bca
reseek -distmx structs.bca -output structs.distmx
muscle -super7 structs.bca -distmxin structs.distmx -reseek -output structs.afa

# deprecated: text mega from pdb2mega
# muscle -align structs.mega -output structs.afa
</pre>

### Downloads and installation

Binary files are self-contained, no dependencies. To install, download the binary and make sure the execute bit is set.

https://github.com/rcedgar/muscle/releases

### Documentation

[Muscle v5 home page](https://drive5.com/muscle5)   
[Manual](https://drive5.com/muscle5/manual)   

### Building MUSCLE from source

[https://github.com/rcedgar/muscle/wiki/Building-MUSCLE](https://github.com/rcedgar/muscle/wiki/Building-MUSCLE)


### References
Edgar RC., Muscle5: High-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny. <i>Nature Communications</i> 13.1 (2022): 6968.    
[https://www.nature.com/articles/s41467-022-34630-w.pdf](https://www.nature.com/articles/s41467-022-34630-w.pdf)

Edgar RC. and Tolstoy I., Muscle-3D: scalable multiple protein structure alignment (2024) <i>BioRxiv</i>.