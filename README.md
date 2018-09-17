These codes are used for fusing the concurrent EEG and fMRI. This repository contains source codes for the fusion of concurrent EEG-fMRI recording.

Here, there are 3 approaches for fusing the EEG and fMRI features, two of them is based on the source-space EEG and one in the sensor-space EEG.

In one of the source-space based approaches, the EEG features are extracted using the source-space ICA as demonstrated in this paper:

https://scholar.google.com.au/scholar?cluster=697165797441660148&hl=en&as_sdt=0,5&sciodt=0,5

For the second source-space based approach, the EEG features are extracted using the enveloping and convolving the source-space projected EEG data as described here:

https://scholar.google.com.au/scholar?cluster=6999230742770998304&hl=en&as_sdt=2005&sciodt=0,5

In the sensor-space based approach, the EEG features are obtained using the temporal ICA on the EEG data.

In the case of fMRI, the date was first preprocessed (brain extraction, movement correction, ...) and then spatial ICA was applied to remove the artifacts. In the second step, the fRMI time-courses were temporally band-pass filtered using the multibanding approach as described here:

https://github.com/YJonmo/Multi-band-ICA-for-EEG-and-MEG

Then temporal ICA was applied to the fMRI data to extract the features.

The fusion of the EEG and fMRI features are based SOBI or PCA (SVD). As SOBI is based on finding the correlation of maximization of the time-courses using a predefined time-shift, we have the option of choosing the time shift to be examined by the SOBI. For small time-shifts (e.g., less than 10 s), the SVD and SOBI behave similarly.
