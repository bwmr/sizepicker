# sizepicker

Density-based cryo-ET particle picker. 

Fork of `cet_pick_size` by the [Bartesaghi Lab](https://gitlab.cs.duke.edu/bartesaghilab/cet_pick_size) using more recent Python and some convenience features. 

### How it works

### Design changes vs. reference implementation:

### Installation

Install via:
`pip install "git+https://github.com/bwmr/sizepicker.git"`

Requires `newstack` from the IMOD package to be on `$PATH` (this is only used for lowpass-filtering and will be removed in future versions).

### Reference

The paper below describes the original implementation:

```@article{JIN2024100104,title = {Accurate size-based protein localization from cryo-ET tomograms},  
journal = {Journal of Structural Biology: X},  
volume = {10},  
pages = {100104},  
year = {2024},  
issn = {2590-1524},  
doi = {https://doi.org/10.1016/j.yjsbx.2024.100104},  
url = {https://www.sciencedirect.com/science/article/pii/S2590152424000096},  
author = {Weisheng Jin and Ye Zhou and Alberto Bartesaghi}
