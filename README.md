# KYM
###### Wavelet transform-based analyzer of biological signals

### Introduction

KYM is a software tool that aims at the quantitative analysis of non-stationary
oscillatory signals. It has been developed within the field of neurosciences to
evaluate the neuronal activity in terms of cytosolic calcium signals, but it can
be considered as a general-purpose set of functions to implement a wavelet-based
time-frequency analysis.

KYM has been developed starting from 2010 and released for the first time in
2011 as supplementary material for the scientific methodological paper [1]. From
then on, KYM has been modified several times in order to upgrade its features
(e.g. new handling methods for the cone of influence, transformed domain
filtering, inverse wavelet transform computation, spatial analysis of the
oscillatory activity) and it has been already employed in many other works (see
[2], [3], [4], [5], [6], [7]).

KYM has been originally developed under the freely redistributable GNU Octave
environment [8], but starting from version 0.5 it has also been made fully
compatible with MATLAB environment and, in particular, it has been tested with
the following releases of these two environments: Octave 3.6.4 and MATLAB 8
(R2012b).

Note for Octave users: Octave has two official plotting backends: Gnuplot (a
classic Linux plotting utility) and FLTK Backend (a new experimental OpenGL
backend based on the FLTK GUI toolkit). We strongly recommend using the first
one by typing

`graphics_toolkit gnuplot`

at octave prompt. Otherwise, you could encounter many serious problems,
especially in the subplot management. To make this change permanent add or
modify the line

`graphics_toolkit('gnuplot');`

in your `<your octave dir>\share\octave\site\m\startup\octaverc` file.

KYM ver.0.5 is the last and most up-to-date release of the software. It consists
of a coordinated set of 19 .m files, each of them containing a single function,
but **only the "upper-case" ones (`VX`, `WT`, `PD`, `FEAT`, `FILT`, `IWT`,
`WTX`) need to be directly managed by the end user**, while the others are just
auxiliary functions invoked on the fly by the first ones.

KYM has a command line user interface that has not been developed taking into
account the end users and, therefore, a future effort will be aimed at
developing a more user-friendly interface.

In order to use KYM (`VX`, `WT` and `WTX` functions in particular), data need to
be stored in a .csv (comma separated values) file, structured as it follows:
**the first column must contain the vector of time samples, while fluorescence
signals must fill all the subsequent columns**. In the case of calcium
signalling analysis, columns are the time courses of the fluorescence intensity
emitted by the calcium indicator loaded inside the cells. Each column may
represent a different cell or a different region of interest (ROI) over the
same cell. The file `data.csv` provides an example of this and you can use it to
test the software. It contains the same set of traces showed in [3] to better
describe the analytical method: 56 calcium intracellular concentration time
courses recorded from ROIs drawn over the same cell, with a sampling time of 2s.

WTX function is specifically tailored for the spatial wavelet analysis of neural
calcium signals described in [3]. It can load two data arguments, but, in this
case, it is assumed that you are using a ratiometric fluorescent indicator for
calcium (such as Fura-2). Further details about this function are given in [3],
in particular see "Surface to Volume Ratio Assessment" subsection.

Many other parameters, beyond the raw data, can be passed as argument to the
main functions, in order to specify the unit of measurement of time, the
changing points of in-acuto treatments, the threshold levels, and so on. For
further information about the meaning and the syntax of each KYM function you
can use help command, followed by the name of a function, to display the
documentation header inside every .m file, or you can follow the examples below
to get started. Nevertheless an exhaustive and complete user guide for KYM has
never been written, so another future effort will be to produce such a
documentation. At present, inspecting the code line by line is the only actual
way to have a full understanding of KYM algorithms.

Actual wavelet transform computation has been implemented as a product in the
Fourier transformed domain. A standard code for this algorithm can be found, for
instance, in WaveLab850 (http://www-stat.stanford.edu/~wavelab/). Peak detection
uses a technique that is based on image dilation (see, for instance,
`localMaximum.m` m-file by Yonathan Nativ,
http://www.mathworks.com/matlabcentral/fileexchange/authors/). The rest of the
code has been written and developed ad hoc to perform the analyses presented
in [1] and [3].

**NOTE:** `imdilate` function is required by `peak`, `paths` (and hence by `PD`
and `WTX`) KYM functions. Octave users can find `imdilate.m` inside the `image`
package. In order to have no problems with function dependences we recommend to
install ALL Octave packages at once, simply by following OcatveForge guidelines
about this. MATLAB users need to have Image Processing Toolbox installed to
enjoy imdilate function.

Since KYM is an open source code you can inspect it to see exactly what
algorithms have been used, and then modify the source to produce a better code
or to satisfy other particular needs. Users may redistribute it and/or modify it
under the terms of the GNU General Public License (GPL) as published by the Free
Software Foundation.

To our knowledge this is the first open source tool specifically dedicated to
the analysis of the time course of cellular calcium signals and, more generally,
of oscillatory signals recorded by means of fluorescent dyes from biological
systems.

### References

1. Ruffinatti F.A., Lovisolo D., Distasi C., Ariano P., Erriquez J., Ferraro M. Calcium signals: analysis in time and frequency domains. *J Neurosci Methods.* 2011 Aug; 199(2):310-320.

2. Zamburlin P., Ruffinatti F.A., Gilardino A., Farcito S., Parrini M., Lovisolo D. Calcium signals and FGF-2 induced neurite growth in cultured parasympathetic neurons: spatial localization and mechanisms of activation. *Pflugers Arch.* 2013 Sep; 465(9):1355-1370.

3. Ruffinatti F.A., Gilardino A., Lovisolo D., Ferraro M. Spatial wavelet analysis of calcium oscillations in developing neurons. *PLoS One.* 2013 Oct 14; 8(10):e75986.

4. Moccia F., Ruffinatti F.A., Zuccolo E. Intracellular Ca2+ Signals to Reconstruct a Broken Heart: Still a Theoretical Approach? *Curr Drug Targets.* 2015 Jul; 16(8):793-815.

5. Gilardino A., Catalano F., Ruffinatti F.A., Alberto G., Nilius B., Antoniotti S., Martra G., Lovisolo D. Interaction of SiO2 nanoparticles with neuronal cells: Ionic mechanisms involved in the perturbation of calcium homeostasis. *Int J Biochem Cell Biol.* 2015 Sep; 66:101-111.

6. Dragoni S., Reforgiato M., Zuccolo E., Poletto V., Lodola F., Ruffinatti F.A., Bonetti E., Guerra G., Barosi G., Rosti V., Moccia F. Dysregulation of VEGF-induced proangiogenic Ca2+ oscillations in primary myelofibrosis-derived endothelial colony-forming cells. *Exp Hematol.* 2015 Dec; 43(12):1019-1030.

7. Lodola F., Laforenza U., Cattaneo F., Ruffinatti F.A., Poletto V., Massa M., Tancredi R., Zuccolo E., Kheder D., Riccardi A., Biggiogera M., Rosti V., Guerra G., Moccia F. VEGF-induced intracellular Ca2+ oscillations are down-regulated and do not stimulate angiogenesis in breast cancer-derived endothelial colony forming cells. *Oncotarget. 2017.* IN PRESS

8. Eaton JW, Bateman D, Hauberg S (2008) GNU Octave Manual Version 3. UK: Network Theory Ltd.


### Syntax to get started


#### Trace Visualization

 General Syntax:

 `VX(filename,sub,unit,mark,roi,output)`

> Example with argument comments:
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

> Example with argument comments:
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

#### Wavelet Analysis and Peak Detection:

General Syntax:

`ratio = PD(wt,par,sig,thr1,thr2,output)`

> Example with argument comments:
>
> `PD(wt,par,sig,50,35);`
>
> * The first three variables (`wt`, `par`, and `sig`) are the outputs of `WT`
    functions concerning the last ROI analyzed by `WT` (ROI 7 in this case).
>
> * Peak detection is performed after a 50% thresholding procedure.
>
> * Frequency paths are detected on scaleograms after a 35% thresholding
    procedure.
>
> * The output consists in several plots deriving the power and the energy
    density of the signal starting from the related wavelet scaleograms: further
    details about them can be found in [1]. Notice that the time-averaged
    activity values and related *r* values are calculate respect to the time
    markers set in `WT`.
