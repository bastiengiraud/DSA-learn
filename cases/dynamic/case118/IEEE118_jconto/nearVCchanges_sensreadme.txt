This file describes in detail how the standard IEEE 118 test data was
slightly modified in order to perform the studies described in detail
in reference [1].  This paper also used the Point of Collapse power
flow [2], which requires the creation of a "k" file.  The information
in this writeup should assist users in the creation of a proper
k-file.  We have also included a sample k-file in this directory.
The intent of this file is to permit anyone seriouly interested
in replicating the results in [1] to do so.

The studies performed in [1] ignore all line flow limits and all generation
active power limits, and concentrate instead on the determination of
maximum loadability based on voltage collapse considerations.


Reactive power limits.

The Q-limits under three operating conditions are described
below.  The first column refers to those Q-limits that are enforced
during an ordinary base case solution of the 118 bus test case "as given."
The second column corresponds to a list of those Q-limits that are enforced
after the initial load expansion described in [1].  The third set corresponds
to those limits in effect at the maximum loadability point.

Case as given       Loads expanded      Point of Collapse
Bus Limit           Bus Limit           Bus Limit
 19 Qmin              1 Qmax      	  1 Qmax
 32 Qmin             12 Qmax      	  6 Qmax
 34 Qmin             15 Qmax      	 12 Qmax
 92 Qmin             18 Qmax      	 15 Qmax
103 Qmax             34 Qmax      	 18 Qmax
105 Qmin             36 Qmax      	 19 Qmax
		     55 Qmax      	 32 Qmax
		     56 Qmax      	 34 Qmax
		     62 Qmax      	 36 Qmax
		     70 Qmax      	 49 Qmax
		     74 Qmax      	 55 Qmax
		     76 Qmax      	 56 Qmax
		     77 Qmax      	 59 Qmax
		     85 Qmax      	 62 Qmax
		     92 Qmax      	 70 Qmax
		    103 Qmax      	 74 Qmax
		    104 Qmax      	 76 Qmax
		    110 Qmax      	 77 Qmax
		      6 Qmax      	 85 Qmax
		     19 Qmax      	 92 Qmax
		     49 Qmax      	103 Qmax
		    100 Qmax      	104 Qmax
		    105 Qmax      	105 Qmax
				      	110 Qmax
				      	  8 Qmax
				      	100 Qmax
				      	 80 Qmax
				      	 65 Qmax
				      	 10 Qmax
				      	 46 Qmax

The Qmax limit of Generator 4 has been increased from the given value
to avoid an immediate instability.
	
Area interchanges:

The modified case considers inter-area power dispatch, as follows:

  Main Area to Area 2 -> 102.90 MW
  Main Area to Area 3 -> 143.80 MW
  Main Area Total     -> 246.70 MW

Assumption about load variaton:
Bus   1, delta P-> 0.0930946  delta Q-> 0.0492854
Bus   2, delta P-> 0.0365077  delta Q-> 0.0164285
Bus   3, delta P-> 0.07119    delta Q-> 0.0182538
Bus   4, delta P-> 0.0547615  delta Q-> 0.0219046
Bus   6, delta P-> 0.09492    delta Q-> 0.0401584
Bus   7, delta P-> 0.0346823  delta Q-> 0.00365077
Bus  11, delta P-> 0.127777   delta Q-> 0.0419838
Bus  12, delta P-> 0.085793   delta Q-> 0.0182538
Bus  13, delta P-> 0.062063   delta Q-> 0.0292061
Bus  14, delta P-> 0.0255554  delta Q-> 0.00182538
Bus  15, delta P-> 0.164285   delta Q-> 0.0547615
Bus  16, delta P-> 0.0456346  delta Q-> 0.0182538
Bus  17, delta P-> 0.0200792  delta Q-> 0.00547615
Bus  18, delta P-> 0.109523   delta Q-> 0.062063
Bus  19, delta P-> 0.0821423  delta Q-> 0.0456346
Bus  20, delta P-> 0.0328569  delta Q-> 0.00547615
Bus  21, delta P-> 0.0255554  delta Q-> 0.0146031
Bus  22, delta P-> 0.0182538  delta Q-> 0.00912692
Bus  23, delta P-> 0.0127777  delta Q-> 0.00547615
Bus  27, delta P-> 0.113174   delta Q-> 0.02373
Bus  28, delta P-> 0.0310315  delta Q-> 0.0127777
Bus  29, delta P-> 0.0438092  delta Q-> 0.00730154
Bus  31, delta P-> 0.0784915  delta Q-> 0.0492854
Bus  32, delta P-> 0.107698   delta Q-> 0.0419838
Bus  33, delta P-> 0.0419838  delta Q-> 0.0164285
Bus  34, delta P-> 0.107698   delta Q-> 0.04746
Bus  35, delta P-> 0.0602377  delta Q-> 0.0164285
Bus  36, delta P-> 0.0565869  delta Q-> 0.0310315
Bus  39, delta P-> 0.0492854  delta Q-> 0.0200792
Bus  40, delta P-> 0.0365077  delta Q-> 0.0419838
Bus  41, delta P-> 0.0675392  delta Q-> 0.0182538
Bus  42, delta P-> 0.0675392  delta Q-> 0.0419838
Bus  43, delta P-> 0.0328569  delta Q-> 0.0127777
Bus  44, delta P-> 0.0292061  delta Q-> 0.0146031
Bus  45, delta P-> 0.0967453  delta Q-> 0.0401584
Bus  46, delta P-> 0.0511107  delta Q-> 0.0182538
Bus  48, delta P-> 0.0365077  delta Q-> 0.0200792
Bus  49, delta P-> 0.158808   delta Q-> 0.0547615
Bus  50, delta P-> 0.0310315  delta Q-> 0.00730154
Bus  51, delta P-> 0.0310315  delta Q-> 0.0146031
Bus  52, delta P-> 0.0328569  delta Q-> 0.00912692
Bus  53, delta P-> 0.0419838  delta Q-> 0.0200792
Bus  54, delta P-> 0.206268   delta Q-> 0.0584123
Bus  55, delta P-> 0.114999   delta Q-> 0.0401584
Bus  56, delta P-> 0.153332   delta Q-> 0.0328569
Bus  57, delta P-> 0.0219046  delta Q-> 0.00547615
Bus  58, delta P-> 0.0219046  delta Q-> 0.00547615
Bus  59, delta P-> 0.505631   delta Q-> 0.206268
Bus  60, delta P-> 0.14238    delta Q-> 0.00547615
Bus  62, delta P-> 0.140555   delta Q-> 0.0255554
Bus  66, delta P-> 0.07119    delta Q-> 0.0328569
Bus  67, delta P-> 0.0511107  delta Q-> 0.0127777
Bus  70, delta P-> 0.120475   delta Q-> 0.0365077
Bus  74, delta P-> 0.124126   delta Q-> 0.0492854
Bus  75, delta P-> 0.085793   delta Q-> 0.0200792
Bus  76, delta P-> 0.124126   delta Q-> 0.0657138
Bus  77, delta P-> 0.111348   delta Q-> 0.0511107
Bus  78, delta P-> 0.129602   delta Q-> 0.04746
Bus  79, delta P-> 0.07119    delta Q-> 0.0584123
Bus  80, delta P-> 0.2373     delta Q-> 0.04746
Bus  82, delta P-> 0.0985707  delta Q-> 0.0492854
Bus  83, delta P-> 0.0365077  delta Q-> 0.0182538
Bus  84, delta P-> 0.0200792  delta Q-> 0.0127777
Bus  85, delta P-> 0.0438092  delta Q-> 0.0273808
Bus  86, delta P-> 0.0383331  delta Q-> 0.0182538
Bus  88, delta P-> 0.0876184  delta Q-> 0.0182538
Bus  90, delta P-> 0.14238    delta Q-> 0.0766661
Bus  92, delta P-> 0.11865    delta Q-> 0.0182538
Bus  93, delta P-> 0.0219046  delta Q-> 0.0127777
Bus  94, delta P-> 0.0547615  delta Q-> 0.0292061
Bus  95, delta P-> 0.0766661  delta Q-> 0.0565869
Bus  96, delta P-> 0.0693646  delta Q-> 0.0273808
Bus  97, delta P-> 0.0273808  delta Q-> 0.0164285
Bus  98, delta P-> 0.062063   delta Q-> 0.0146031
Bus 100, delta P-> 0.0675392  delta Q-> 0.0328569
Bus 101, delta P-> 0.0401584  delta Q-> 0.0273808
Bus 102, delta P-> 0.00912692 delta Q-> 0.00547615
Bus 103, delta P-> 0.0419838  delta Q-> 0.0292061
Bus 104, delta P-> 0.0693646  delta Q-> 0.0456346
Bus 105, delta P-> 0.0565869  delta Q-> 0.04746
Bus 106, delta P-> 0.0784915  delta Q-> 0.0292061
Bus 107, delta P-> 0.0511107  delta Q-> 0.0219046
Bus 108, delta P-> 0.00365077 delta Q-> 0.00182538
Bus 109, delta P-> 0.0146031  delta Q-> 0.00547615
Bus 110, delta P-> 0.07119    delta Q-> 0.0547615
Bus 112, delta P-> 0.0456346  delta Q-> 0.02373
Bus 114, delta P-> 0.0146031  delta Q-> 0.00547615
Bus 115, delta P-> 0.0401584  delta Q-> 0.0127777
Bus 117, delta P-> 0.0365077  delta Q-> 0.0146031
Bus 118, delta P-> 0.0602377  delta Q-> 0.0273808


Generator Participation by Area:

As demands increase, generators particiapate in the following
proportions:

  Area 2, Bus  10 -> 0.420954
  Area 2, Bus  12 -> 0.0795136
  Area 2, Bus  25 -> 0.2058
  Area 2, Bus  26 -> 0.293732

  Area 1, Bus  46 -> 0.0100903
  Area 1, Bus  49 -> 0.108338
  Area 1, Bus  54 -> 0.0254913
  Area 1, Bus  59 -> 0.0823154
  Area 1, Bus  61 -> 0.0849712
  Area 1, Bus  65 -> 0.207648
  Area 1, Bus  66 -> 0.208179
  Area 1, Bus  69 -> 0.272969

  Area 3, Bus  80 -> 0.337819
  Area 3, Bus  89 -> 0.429887
  Area 3, Bus 100 -> 0.17847
  Area 3, Bus 103 -> 0.0283287
  Area 3, Bus 111 -> 0.0254957


Modification necessary to the standard 118 bus data include:

(1) Increase the loading from the initial loading in the 118 bus
    test system data file to a total load of 5677 MW in direct
    proportion to load size.  While doing this, enforce all reactive
    power limits.  At this loading, there will be no generators
    at their minimum limit.

(2) Remove the Q max limit at bus 4.  If this is not done, an immediate
    instability results just before voltage collapse.

(3) Establish an area interchange control as specified above.

Other changes are necessary for some of the sensitivity simulations
as a result of the need for model changes required by each sensitivity
study.  For example, to study the effect of load type the load model
themselves change.  The paper details these changes.

References:

[1] S. Greene, I. Dobson and F. Alvarado, "Sensitivity of the loading
    margin to voltage collapse with respect to arbitrary parameters,"
    ECE Report 95-8, Department of Electrical and Computer Engineering,
    the University of Wisconsin, Madison, Wisconsin, July 1995.

[2] C. Canizares, F. Alvarado et al, "Point of Collapse and Continuation
    Methods for Large AC/DC Systems," IEEE Transactions on Power Systems, 
    Vol. 7, No. 1, Feb 1993, pp. 1-8.
