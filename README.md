> **Important**: This repository is *silent*. No active development is currently being
carried out, but it is still maintained. Please open an issue to request help
or features.

# kym
###### Wavelet transform-based analyzer of biological signals

### Introduction

*kym* is a software tool for the quantitative analysis of non-stationary
oscillatory signals. It has been developed within the field of neurosciences to
study neuronal activity in terms of cytosolic calcium signals, but it can be
considered as a general-purpose set of MATLAB functions that implement a
wavelet-based time/frequency analysis.

*kym* has been developed starting from 2010 and released for the first time in
2011 as supplementary material for a scientific methodological paper ([1]). From
then on, *kym* has been modified several times to extend its analysis
capabilities, and different methods to handle the cone of influence as well as
many filters in the transformed domain are now available and already used in
many other works (see e.g., [2], [3], [4], [5], [6], [7]).

*kym* ver.0.6 is the last and most up-to-date release of the software. It
consists of a collection of 30 .m files, each of them containing a single MATLAB
function, but **only the "capital" ones (`VX`, `WT`, `WTX`, `IWT`, `KYM`,
`FILT*`, `FEAT`) need to be directly managed by the end user**, while the other
ones are to be considered as auxiliary functions, invoked on the fly by the
previous ones.

In order to use *kym* (`VX`, `WT` and `WTX` functions in particular), data need
to be stored in a .csv (comma separated values) file, organized as it follows:
**the first column must contain the vector of time samples, while the actual
signals should fill all the other columns**. In the case of calcium signal
analysis, columns are the time courses of the fluorescence intensity emitted by
the calcium indicator loaded inside the cells. Each column may represent a
different cell or a different region of interest (ROI) over the same cell. The
file `data.csv` provides an example of this and it can be used to test the
software. It contains the same set of traces showed in [3] to present the
analytical method: 56 time traces of calcium intracellular concentration
recorded from ROIs drawn over the same cell, with a sampling time of *dt = 2 s*.

`WTX` function is specifically tailored for the spatial wavelet analysis of
neural calcium signals as described in [3]. It can take two data arguments when
using a ratiometric fluorescent indicator for calcium (such as Fura-2). Further
details about this function are given in [3], in particular see "Surface to
Volume Ratio Assessment" subsection.

Many other parameters, beyond the raw data, can be passed as argument to the
main functions, in order to specify the unit of measurement of time, the
changing points of *in acuto* treatments, the threshold levels, and so on. For
further information about the meaning and the syntax of each *kym* function you
can use the `help` command followed by the name of the function, or you can
follow the examples below to get started. Nevertheless an exhaustive and
complete user guide for *kym* has never been written, so a future effort will be
to produce such a documentation.

From a technical point of view, wavelet transform computation has been
implemented as a product in the Fourier transformed domain. A standard code for
this algorithm can be found, for instance, in WaveLab850
(http://www-stat.stanford.edu/~wavelab/). Peak detection is based on image
dilation (see e.g., `localMaximum.m` m-file by Yonathan Nativ,
http://www.mathworks.com/matlabcentral/fileexchange/authors/). The rest of the
code has been written and developed *ad hoc* to perform the analyses presented
in [1] and [3].

**NOTE:** `imdilate` function is required by `peak`, `paths` (and hence by `KYM`
and `WTX`) *kym* functions. Octave users can find `imdilate.m` within the
`image` package, while MATLAB users need to have Image Processing Toolbox
installed.

### References

1. Ruffinatti F.A., Lovisolo D., Distasi C., Ariano P., Erriquez J., Ferraro M. Calcium signals: analysis in time and frequency domains. *J Neurosci Methods.* 2011 Aug; 199(2):310-320.

2. Zamburlin P., Ruffinatti F.A., Gilardino A., Farcito S., Parrini M., Lovisolo D. Calcium signals and FGF-2 induced neurite growth in cultured parasympathetic neurons: spatial localization and mechanisms of activation. *Pflugers Arch.* 2013 Sep; 465(9):1355-1370.

3. Ruffinatti F.A., Gilardino A., Lovisolo D., Ferraro M. Spatial wavelet analysis of calcium oscillations in developing neurons. *PLoS One.* 2013 Oct 14; 8(10):e75986.

4. Moccia F., Ruffinatti F.A., Zuccolo E. Intracellular Ca2+ Signals to Reconstruct a Broken Heart: Still a Theoretical Approach? *Curr Drug Targets.* 2015 Jul; 16(8):793-815.

5. Gilardino A., Catalano F., Ruffinatti F.A., Alberto G., Nilius B., Antoniotti S., Martra G., Lovisolo D. Interaction of SiO2 nanoparticles with neuronal cells: Ionic mechanisms involved in the perturbation of calcium homeostasis. *Int J Biochem Cell Biol.* 2015 Sep; 66:101-111.

6. Dragoni S., Reforgiato M., Zuccolo E., Poletto V., Lodola F., Ruffinatti F.A., Bonetti E., Guerra G., Barosi G., Rosti V., Moccia F. Dysregulation of VEGF-induced proangiogenic Ca2+ oscillations in primary myelofibrosis-derived endothelial colony-forming cells. *Exp Hematol.* 2015 Dec; 43(12):1019-1030.

7. Lodola F., Laforenza U., Cattaneo F., Ruffinatti F.A., Poletto V., Massa M., Tancredi R., Zuccolo E., Kheder D., Riccardi A., Biggiogera M., Rosti V., Guerra G., Moccia F. VEGF-induced intracellular Ca2+ oscillations are down-regulated and do not stimulate angiogenesis in breast cancer-derived endothelial colony forming cells. *Oncotarget. 2017.* IN PRESS

### Syntax to get started

#### Trace Visualization

General Syntax:

`VX(filename,sub,unit,mark,roi,output)`

> Example with comments on function arguments:
>
> `VX('data.csv','C','s',[110],[]);`
>
> `VX` allows a preliminary inspection of time courses of the raw data
  contained in `data.csv` and returns a graph of sampling quality control.
>
> * `'C'` stands for *Color map mode*: it displays the two-variable function
    \[Ca<sup>2+</sup>\]\(*t*, *x*\) as a color map (*t* = time, *x* = space).
    Different visualization modes are available:
>  * *`n`* (being an integer greater or equal to one) shows the selected traces
     arranging them in subplot groups of *n* traces per window.
>  * `'S'` stands for "Superimposition" and shows all the traces selected
     superimposed according to a color gradient from blue (first ROIs) to red
     (last ROIs).
>  * `'C'` is the already mentioned Color map mode.
>  * `'M'` shows a single trace representing the point-by-point mean and
     standard error of the mean (SEM) computed from all the selected trace.
>
> * Second (`'s'`) is the unit of measurement of the first column entries.
>
> * A time marker is placed at step number 110, that corresponds in this case to
    *t* = 220s being *dt* = 2s the sampling time.
>
> * The last void argument (`[]`) means that all traces included in the file are
    selected and then displayed. As an alternative, you can pass an array of
    ROIs according to MATLAB syntax. E.g. `[3 5 9:12 15:17]` means ROIs 3, 5, 9,
    10, 11, 12, 15, 16, 17.
 
#### Wavelet Transform

General syntax:

`[wt,par,sig] = WT(filename,unit,mark,roi,lft,nvoice,wavelet,cone,output)`

> Example with comments on function arguments:
>
> `[wt,par,sig]=WT('data.csv','s',[110],[2:7],3,48,'M1','PAD');`
>
> Computes the wavelet transform starting from the raw data arranged in the file
  `data.csv`.
>
> * Second (`'s'`) is the unit of measurement of the first column entries.
>
> * A time marker is placed at step number 110, that corresponds in this case to
    *t* = 220s being *dt* = 2s the sampling time.
>
> * The wavelet transform and related scaleogram are computed just for ROIs from
    2 to 7.
>
> * The first 3 lowest frequencies are discarded (low frequency threshold,
    recommended).
>
> * Scaleograms are drawn using 48 inter-octave voices.
>
> * `'M1'` (Morlet version 1) mother wavelet is used to compute wavelet
    transform.
>
> * `'PAD'` technique is use to handle the cone of influence. Three different
    ways to handle the cone of influence (COI) issue have been implemented:
>  * `'NONE'`: this option does nothing on COI artifacts, but it is useful to
     preserve the whole information carried by the original signal within the
     scaleogram.
>  * `'COI'`: this option computes the e-folding cone of influence and masks it
     in order to completely discard those regions out of the scaleogram.
>  * `'PAD'`: this is a less drastic algorithm that provides a linear detrending
     ensuring a matching between the starting and the ending point of the signal
     and, at the same time, pads the original signal with a number of zeros
     equal to the signal length, finally removing them after wavelet
     transformation: this is the procedure we recommend because it allows to
     avoids artificial edge discontinuities and low-frequencies cross talk
     between the two edges of the scaleogram, with a very poor loss of
     information.
>
> * Output variables (`wt`, `par`, `sig`) contain scaleograms, generic
    descriptive parameters and the original signal respectively, concerning last
    ROI analyzed (ROI 7 in this case). They serves as input argument of `PD`
    function.

#### Wavelet Analysis and Peak Detection

General Syntax:

`ratio = KYM(wt,par,sig,thr,output)`

> Example with comments on function arguments:
>
> `KYM(wt,par,sig,50);`
>
> * The first three variables (`wt`, `par`, and `sig`) are the outputs of `WT`
    functions concerning the last ROI analyzed by `WT` (ROI 7 in this case).
>
> * Peak detection and frequency path detection are performed after a 35%
    thresholding of the scaleogram.
>
> * The output consists in several plots deriving the power and the energy
    density of the signal starting from the related wavelet scaleograms: further
    details about them can be found in [1]. Notice that the time-averaged
    activity values and related *r* values are calculate respect to the time
    markers set in `WT`.
