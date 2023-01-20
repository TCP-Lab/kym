# KYM
###### Wavelet transform-based analyzer of biological signals

KYM is a GNU Octave set of functions that implement an automated procedure of
quantification and characterization of oscillatory biological signals. In
particular it implements a wavelet-based analysis for the oscillatory component
of the cytosolic calcium concentration time course.

KYM is tested on both Microsoft and Linux Debian systems, under Octave 3.2.3
version.

At the present KYM is made up of 14 .m files each containing a single function,
for a total of about 2000 code lines, but only 4 of them (`VX`, `WT`, `PD`,
`FEAT`) need to be directly managed by the end user. The other ones are
auxiliary functions invoked on the fly by the mains.

KYM has a command-line user-interaction that has not been developed taking into
account the end users and, therefore, a future effort will be to make a
user-friendly interface.

In order to use KYM, data need to be stored in a .csv (comma separated values)
file. It must contain the vector of time samples as first column, while the time
courses of fluorescence value for each cell (or region of interest - ROI) must
fill the next columns in the matrix. Some parameters can be passed as argument
to the main functions, in order to specify the unit of measurement of time, the
changing points of *in acuto* treatments, the threshold levels, and so on.

Wavelet actual computation consists in the usual method that implements a time
convolution as a product in the Fourier transformed domain; a standard code for
this algorithm can be downloaded at http://www-stat.stanford.edu/~wavelab/.

Peaks detection uses a technique that is based on images dilation (see for
instance http://www.mathworks.com/matlabcentral/fileexchange/authors/26510/).

The rest of code has been written and developed *ad hoc* to perform the analysis
presented in "Ruffinatti F.A. et al., Calcium signals: analysis in time and
frequency domains. *Journal of Neuroscience Methods* 2011; 199:310-320".

At the beginning we made KYM available as supplementary material for this early
article, but the up-to-date version of KYM (ver. 0.3) offers some new features,
among which 3 different ways to handle the cone of influence (COI) issue: `NONE`
(this option does nothing on COI artefacts, but it is useful to preserve in the
scaleogram the whole information carried by the original signal), `COI` (this
option computes the e-folding cone of influence and masks it in order to
completely discard those regions of the scaleogram) and `PAD` (this is a less
drastic algorithm that provides a linear detrending ensuring a matching between
the starting and the ending point of the signal and, at the same time, pads the
original signal with a number of zeros equal to the signal length, finally
removing them after wavelet transformation: this is the procedure we recommend
because it allows to avoids artificial edge discontinuities and low-frequencies
cross talk between the two edges of the scaleogram, with a very poor loss of
information).

As it is an open source code one can inspect it to see exactly what algorithms
have been used, and then modify the source to produce a better code or to
satisfy other particular needs. Users may redistribute it and/or modify it under
the terms of the GNU General Public License (GPL) as published by the Free
Software Foundation. Because KYM is a free software users are encouraged to help
make Octave more useful by writing and contributing additional functions for it
and by reporting any problems they may have.

Each .m file comes with an extended documentation in the heading (displayable
using 'help' function), explaining the syntax and the meaning of each KYM
function and its related arguments.
