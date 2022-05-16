# AI-SRS-multimets
Work developed from my PhD project

The directory TESTS contains the statistical tests and the supervised and unsupervised classification algorithms: Decision tree, K-means and PCA.

Yhe directory SOFTWARE contains the "heart" of the project. The software was built in Matlab v2019.b. It is based on a GUI called VisorIPRO.

The other files complements the GUI and they are specific to perform rotations, translations of the structures of the DICOM files. Other files are based on genetic algorithms.

The SOFTWARE was validated with the DICOM files generated on the treatment planning system Brainlab Elements TM v2.0 and v3.0.

The SOFTWARE works with DICOM files prepared in the following form:

1. The DICOM files have to have RS, RT, RP, RD or CT as head.
2. Create a directory called Patient.
3. Inside the directory, create two directories.
4. The first one has to call CT. Put the CT files on this directory.
5. The second one has to call Dose. Put the RT, RP, RD and RS files on this directory.
