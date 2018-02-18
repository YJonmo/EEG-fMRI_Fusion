# EEG-fMRI_Fusion
These codes are used for fusing the cuncurrent EEG and fMRI
This repository contains source codes for the fusion of concurrent EEG-fMRI recording.

Here, there are 3 approaches in fusing the EEG and fMRI features, two of them is based on the source-sapce EEG and one in the sensor-sapce EEG.

In one of the source-space based approaches, the EEG features are extraced using the source-space ICA as demonstraited in this paper: 

https://scholar.google.com.au/scholar?cluster=697165797441660148&hl=en&as_sdt=0,5&sciodt=0,5

For the second source-space based approach, the EEG features are extracted using the envelopping and convolving the soure-space projected EEG data as described here:

https://scholar.google.com.au/scholar?cluster=6999230742770998304&hl=en&as_sdt=2005&sciodt=0,5

In the sensor-space based approach, the EEG features are obtained using the temporal ICA on the EEG data.

In the case of fMRI, the date was first preporcessed (brain extraction, movement correction, ...) and then spatial ICA was applied to remove the artifacts. In the second step, the fRMI time-courses were temporally band-pass filtered using multibanding approach as described here:

https://github.com/YJonmo/Multi-band-ICA-for-EEG-and-MEG

Then temporal ICA was applied on the fMRI data to extract the features.

The fusion of the EEG and fMRI features are based SOBI or PCA (SVD). As SOBI is based on finding the correlation of maximization of the time-courses using a predefined time-shift, we have the option of choosing the time shift to be examined by the SOBI. For small time-shifts (e.g., less than 10 s), the SVD and SOBI behave similarly. 


