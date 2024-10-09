[My version of the IEEE 118 bus system for dynamic studies]
IEEE118 - Power System + Dynamic data

Set-up for PSSe v.33:
1- Create a v.33 sav case from a raw file "IEEE118_v33.raw", save it as "IEEE118.sav".
2- run dynFLAT.PY from within PSSe:
   It creates a converted *.cnv file using the "conl.idv" and
   a snap *.snp file, using the "ieee118.dyr" file
   It performs a no-disturbance test, creating an *.out file
   with channels (monitored vars) defined from within.
3- run GETDEVNS.IDV from within PSSPLT, the plotting tool in PSSe
  It output to a file the max deviations for selected channels
 
The original IEEE 118-bus can be translated from its cdf version in "IEEE118cdf.txt"
or from "ieee118.doc", "ieee118.xlsx".  A close version is in "ieee118_v32.raw"
