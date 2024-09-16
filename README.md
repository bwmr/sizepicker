# sizepicker

Density-based cryo-ET particle picker. 

Fork of `cet_pick_size` by the [Bartesaghi Lab](https://gitlab.cs.duke.edu/bartesaghilab/cet_pick_size) using more recent Python and some convenience features. 

### How it works

Run as `sizepicker input_files_*.mrc output_folder`.  
Full help via `sizepicker --help`.

1. The tomogram will be low-pass filtered to half the particle radius, controlled by `--radius`. The pixel size is automatically determined from the mrc header, or given via `--angpix`. The filtered tomogram is written to the output folder.
2. A contamination mask is created for each tomogram:
    - The tomogram is binned to speed up processing (set by `--contamination-binning`).
    - If `--gaussian` is passed, a low-pass filtered version of the tomogram is used to flatten the background (controlled by `--sigma_cont`).
    - Everything that is more dense than n*SD below mean (set by `--stdtimes_cont`) is considered contamination.
    - This is smoothed, small patches are removed (set by `--minsize_cont`), holes are removed and the mask is dilated (set by `--dilate_cont`).
    - The contamination mask is also written to the output folder. 
3. Local minima outside the contamination mask are picked, in a local neighborhood of `--radiustimes` * particle radius. If `--inhibit` is passed, a slightly more refined algorithm is used, usually giving more picks. Only coordinates within the `--detection_z` voxels around the center are considered.
4. Subvolumes are extracted from the lowpass-filtered tomogram and z-projected. Only picks, where the SD within the particle is `--stdtimes_pick` * SD above the mean particle density SD are retained. Additionally, only particles with a higher foreground SD than background SD are kept, unless `--remove_edge` is passed.
5. XYZ coordinates of the retained picks are written to the output folder. 

### Design changes vs. reference implementation:
- added ability to loop over several tomograms
- using `mrcfile` to parse tomograms, removing requirement for Python 3.6.
- using `click` for CLI
- removed `--binning` parameter, instead read pixel size out of header
- `--sigma_cont` is now given in A rather than px
- choice of a few of the default values changed

### Installation

Install via:
`pip install "git+https://github.com/bwmr/sizepicker.git"`

Requires `newstack` from the IMOD package to be on `$PATH` (it is only used for lowpass-filtering, this requirement will be removed in future versions).

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
