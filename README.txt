Things you should know:

sa.rb contains all the info you need to use the algorithm. The license info is in the file.

noyce.csv, Sample_3_750-800.csv and Sample_3_1000-1050.csv are annotated direct injection lipidomics data sets. They are free to use for any purpose, provided the paper associated with this work (see below) is cited in any publications where the results of said work are discussed. Additionally, if the Noyce data set is used, please also cite the paper (below) that discusses how the data was generated. The first paper below gives info on the Sample_3 datasets. The format of the datasets are:

<m/z>, <intensity>, <rt>, <peak ID>

where peak ID is the unique identifier of the peak the point was assigned to. Peaks with ID = 0 are considered noise.

<paper citations will go here>

bin-point is just an organizational helpfile for sa.rb.

mzML_to_csv will convert an mzML file to the csv format used by sa.rb. You have to have mspire installed for it to work. That is a freely available gem: "gem install mspire". If you use it, please cite the upcoming paper on it.

<mspire paper citation will go here>

Enjoy! Send any questions to Rob 2robsmith@gmail.com.
