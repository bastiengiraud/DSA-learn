%% MATPOWER Case Format : Version 2
function mpc = pglib_opf_case39_epri
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100.0;
%% bus data
%    bus_i    type    Pd    Qd    Gs    Bs    area    Vm    Va    baseKV    zone    Vmax    Vmin
mpc.bus = [
	1	1	78.08	35.36	0.0	0.0	2	1.0118964863737825	-3.9709910327478166	345.0	1	1.06	0.94				
	2	1	0.0	0.0	0.0	0.0	2	1.0361640849635887	-1.8690313829285197	345.0	1	1.06	0.94				
	3	1	257.6	1.9200000000000002	0.0	0.0	2	1.0261883116764676	-5.740113468215608	345.0	1	1.0554	0.94				
	4	1	400.0	147.20000000000002	0.0	0.0	1	1.0173098456700833	-8.317386674088182	345.0	1	1.0479	0.94				
	5	1	0.0	0.0	0.0	0.0	1	1.023388123160793	-7.867646658862885	345.0	1	1.0567	0.9439				
	6	1	0.0	0.0	0.0	0.0	1	1.0265247442640153	-7.505190023067316	345.0	1	1.06	0.9469				
	7	1	187.04000000000002	67.2	0.0	0.0	1	1.0159580876766972	-8.845044138976156	345.0	1	1.0506	0.94				
	8	1	417.6	141.28	0.0	0.0	1	1.0138873991176334	-9.052110051014683	345.0	1	1.0489	0.94				
	9	1	5.2	-53.279999999999994	0.0	0.0	1	1.018875200017187	-6.298306665615597	345.0	1	1.06	0.9419				
	10	1	0.0	0.0	0.0	0.0	1	1.035150963907695	-7.315457852876776	345.0	1	1.06	0.9552				
	11	1	0.0	0.0	0.0	0.0	1	1.0314968525055337	-7.381589176719547	345.0	1	1.0593	0.951				
	12	1	6.824	70.4	0.0	0.0	1	1.0224641818801348	-7.434477238431447	345.0	1	1.0497	0.94				
	13	1	0.0	0.0	0.0	0.0	1	1.0314742219185242	-7.3866934063364225	345.0	1	1.0573	0.9498				
	14	1	0.0	0.0	0.0	0.0	1	1.0251068268717614	-7.550328996052732	345.0	1	1.0531	0.9423				
	15	1	256.00000000000006	122.40000000000002	0.0	0.0	3	1.0178612794190642	-6.600260592798639	345.0	1	1.0466	0.94				
	16	1	263.2	25.839999999999996	0.0	0.0	3	1.0263171867993004	-4.913205073313745	345.0	1	1.0546	0.9449				
	17	1	0.0	0.0	0.0	0.0	2	1.0253125823706668	-5.12227876817776	345.0	1	1.0554	0.9406				
	18	1	126.40000000000002	24.0	0.0	0.0	2	1.0247895344668265	-5.706263805135811	345.0	1	1.0546	0.94				
	19	1	0.0	0.0	0.0	0.0	3	1.0526323513288294	0.5188860566726127	345.0	1	1.06	0.9893				
	20	1	544.0	82.4	0.0	0.0	3	0.9847772584224658	0.17344821825115456	345.0	1	1.0067	0.94				
	21	1	219.20000000000002	92.0	0.0	0.0	3	1.0184527595767052	-3.516889881271142	345.0	1	1.0519	0.94				
	22	1	0.0	0.0	0.0	0.0	3	1.0227914249451646	-0.4044390104528665	345.0	1	1.06	0.9447				
	23	1	198.00000000000003	67.68	0.0	0.0	3	1.0200030762435746	-1.8055404736010432	345.0	1	1.0572	0.94				
	24	1	246.88000000000002	-73.76	0.0	0.0	3	1.0291617756101732	-5.159155384055236	345.0	1	1.0587	0.9459				
	25	1	179.20000000000002	37.760000000000005	0.0	0.0	2	1.0287605565249105	-2.149164091535543	345.0	1	1.06	0.94				
	26	1	111.19999999999999	13.600000000000001	0.0	0.0	2	1.0231063317721694	-1.747769185132527	345.0	1	1.06	0.94				
	27	1	224.8	60.400000000000006	0.0	0.0	2	1.0192902566696194	-4.271388921456328	345.0	1	1.054	0.94				
	28	1	164.8	22.080000000000002	0.0	0.0	3	1.006957790163245	3.752589612276562	345.0	1	1.06	0.9414				
	29	1	226.8	21.519999999999996	0.0	0.0	3	1.002739683500088	6.976076503910422	345.0	1	1.06	0.9457				
	30	2	0.0	0.0	0.0	0.0	2	1.028880704422357	4.701501970986763	345.0	1	1.06	0.94				
	31	3	7.359999999999999	3.6799999999999997	0.0	0.0	1	1.0025327813323415	2.143014216208063e-26	345.0	1	1.06	0.94				
	32	2	0.0	0.0	0.0	0.0	1	1.0009800204386916	-6.424357213349509	345.0	1	1.048	0.94				
	33	2	0.0	0.0	0.0	0.0	3	1.0148726231500715	5.186387339895661	345.0	1	1.0266	0.94				
	34	2	0.0	0.0	0.0	0.0	3	0.9896056325043593	5.489216057087287	345.0	1	1.0282	0.94				
	35	2	0.0	0.0	0.0	0.0	3	1.0042964039851696	5.09685939606773	345.0	1	1.06	0.94				
	36	2	0.0	0.0	0.0	0.0	3	1.0201330945141516	-0.23063966197464406	345.0	1	1.06	0.94				
	37	2	0.0	0.0	0.0	0.0	2	1.0039744002582713	-1.1113776033926894	345.0	1	1.06	0.94				
	38	2	0.0	0.0	0.0	0.0	3	0.9738323929806453	14.980347571816104	345.0	1	1.06	0.94				
	39	2	883.1999999999999	200.0	0.0	0.0	1	0.9898419823544102	-4.202379989555321	345.0	1	1.06	0.94				
];

%% generator data
%    bus    Pg    Qg    Qmax    Qmin    Vg    mBase    status    Pmax    Pmin    Pc1    Pc2    Qc1min    Qc1max    Qc2min    Qc2max    ramp_agc    ramp_10    ramp_30    ramp_q    apf
mpc.gen = [
	30	649.6817660355186	140.00001921824492	400.0	140.0	1.0	100.0	1	1040.0	0.0															
	31	509.8652427796324	209.73153507442294	300.0	-100.0	1.0	100.0	1	646.0	0.0															
	32	82.61639313273726	168.49603790769677	300.0	150.0	1.0	100.0	1	725.0	0.0															
	33	591.4766136850048	216.88650255400455	250.0	0.0	1.0	100.0	1	652.0	0.0															
	34	508.0	72.87709803315818	167.0	0.0	1.0	100.0	1	508.0	0.0															
	35	682.4436007623773	77.58478174213114	300.0	-100.0	1.0	100.0	1	687.0	0.0															
	36	159.6265184729694	7.876733466710098e-6	240.0	0.0	1.0	100.0	1	580.0	0.0															
	37	77.20890603498184	0.00030871742395357967	250.0	0.0	1.0	100.0	1	564.0	0.0															
	38	852.0612485879291	-11.862547007201407	300.0	-150.0	1.0	100.0	1	865.0	0.0															
	39	1039.264932233541	-99.99999084090418	300.0	-100.0	1.0	100.0	1	1100.0	0.0															
];

%% branch data
%    f_bus    t_bus    r    x    b    rateA    rateB    rateC    ratio    angle    status    angmin    angmax
mpc.branch = [
	1	2	0.0035	0.0411	0.6987	600.0	600.0	600.0	0	0	1	-16.157409822689214	4.365938398896874	-97.8005377340316	-85.47380694801593	98.21192681157365	17.026064345256714				
	1	39	0.001	0.025	0.75	1000.0	1000.0	1000.0	0	0	1	-4.073729923380153	8.250592249883853	19.720537734031588	50.11380694801593	-19.640228311787283	-123.24563323566571				
	2	3	0.0013	0.0151	0.2572	500.0	500.0	500.0	0	0	1	-2.538203032429547	4.938896194027696	479.1239315711298	29.46367636019593	-476.32166233359607	-24.263583591565922				
	2	25	0.007	0.0086	0.146	500.0	500.0	500.0	0	0	1	-3.409098881028398	3.4320171928336314	80.19599475261852	16.235612573212617	-79.73889010977301	-31.237512212630985				
	2	30	0.0	0.0181	0.0	900.0	900.0	2500.0	1.025	0.0	1	-10.679933301238545	0.1604281826366305	-657.531853135322	-62.72535327866525	657.531853135322	140.00001921824492				
	3	4	0.0013	0.0213	0.2214	500.0	500.0	500.0	0	0	1	-4.285724307578558	6.348372370049521	222.47460326994695	22.496500152110695	-221.84919050890133	-35.36331083416861				
	3	18	0.0011	0.0133	0.2138	500.0	500.0	500.0	0	0	1	-3.1512678732195276	3.179915762976069	-3.7529409363509654	-0.15291656054476846	3.7543760797223915	-22.313538329429672				
	4	5	0.0008	0.0128	0.1342	600.0	600.0	600.0	0	0	1	-3.311696055856158	2.8705185536054243	-66.5872219543504	-50.840564818338926	66.63639092275353	37.65542064981719				
	4	14	0.0008	0.0129	0.1382	500.0	500.0	500.0	0	0	1	-3.40336930307709	0.26929016371148695	-111.56358753674829	-60.9961243474925	111.68221088951115	48.4963013101941				
	5	6	0.0002	0.0026	0.0434	1200.0	1200.0	1200.0	0	0	1	-1.1172677005051053	0.2578310078088704	-263.4793480235146	-104.65736389780893	263.63193508378595	102.08166010770256				
	5	8	0.0008	0.0112	0.1476	900.0	900.0	900.0	0	0	1	0.6130648407899808	3.7929806037660496	196.84295710076108	67.00194324799175	-196.50432697625166	-77.57676841117829				
	6	7	0.0006	0.0092	0.113	900.0	900.0	900.0	0	0	1	0.6474423084978302	3.477853816444097	271.80049239535595	97.32103573426922	-271.31912052380704	-101.72547032924935				
	6	11	0.0007	0.0082	0.1389	480.0	480.0	480.0	0	0	1	-2.5553917662834715	0.8479775367936184	-32.92718469950949	-66.72121800837515	32.95782786650218	52.37248480492123				
	6	31	0.0	0.025	0.0	1800.0	1800.0	1800.0	1.07	0.0	1	-11.092462913732737	0.29220847551671986	-502.5052427796324	-132.68147783359663	502.5052427796324	206.05153507442293				
	7	8	0.0004	0.0046	0.078	900.0	900.0	900.0	0	0	1	-0.24064227395494575	1.2776958831417358	84.279120523807	34.525470329249345	-84.245834774708	-42.17722434122036				
	8	9	0.0023	0.0363	0.3804	900.0	900.0	900.0	0	0	1	-6.600473799907084	11.676879864766178	-136.84983824904035	-21.526007247601356	137.26894757126522	-11.156088100333752				
	9	39	0.001	0.025	1.2	900.0	900.0	900.0	0	0	1	-4.772738433439757	7.849521793292279	-142.46894757126523	64.43608810033375	142.81916196684816	-176.75435760523848				
	10	11	0.0004	0.0043	0.0729	600.0	600.0	600.0	0	0	1	-0.4526366581533504	1.3349916626548182	36.529120168175794	80.67916878291203	-36.497431176677644	-88.12249446097898				
	10	13	0.0004	0.0043	0.0729	600.0	600.0	600.0	0	0	1	-0.6302535746439055	1.5126085791453734	38.772087461558804	81.01796445398287	-38.73955360206427	-88.45203762972058				
	10	32	0.0	0.02	0.0	900.0	900.0	2500.0	1.07	0.0	1	-9.963736057325017	0.13178029288008933	-75.3012076297346	-161.69713323689493	75.3012076297346	168.49603790769677				
	12	11	0.0016	0.0435	0.0	500.0	500.0	500.0	1.006	0.0	1	-0.7677634454753032	0.888084582452776	-3.520195652148627	-35.222363953453275	3.5396033101754694	35.75000965605775				
	12	13	0.0016	0.0435	0.0	500.0	500.0	500.0	1.006	0.0	1	-1.008405719430249	0.6531718864491385	-3.303804347851373	-35.17763604654672	3.323140362126325	35.703333934646196				
	13	14	0.0009	0.0101	0.1723	600.0	600.0	600.0	0	0	1	-1.5126085791453734	3.552338329811104	35.41641323993795	52.7487036950744	-35.3733754446183	-70.48458020023254				
	14	15	0.0018	0.0217	0.366	600.0	600.0	600.0	0	0	1	-4.927437038125079	6.508800552686152	-76.30883544489284	21.988278890038448	76.43768081853946	-58.62498263865934				
	15	16	0.0009	0.0094	0.171	600.0	600.0	600.0	0	0	1	-3.684118622691193	1.4495832216809827	-332.4376808185395	-63.775017361340694	333.4239115603337	56.211548433301395				
	16	17	0.0007	0.0089	0.1342	600.0	600.0	600.0	0	0	1	-2.7960340402384176	3.50077212824933	43.79078292950795	1.151425057140391	-43.77759013994458	-15.105506248420056				
	16	19	0.0016	0.0195	0.304	600.0	600.0	2500.0	0	0	1	-6.90414143132642	6.778090716397639	-530.2190345526326	-86.12608151577398	534.5640999430259	106.22886606945121				
	16	21	0.0008	0.0135	0.2548	600.0	600.0	600.0	0	0	1	-4.004974987964455	2.9450030669724314	-184.34501334030716	59.59199344741434	184.64360055287713	-81.18723414583701				
	16	24	0.0003	0.0059	0.068	600.0	600.0	600.0	0	0	1	-0.7219268218648373	1.5412564689019146	74.14935340309813	-56.66888542208214	-74.1256672619054	49.952209709330006				
	17	18	0.0007	0.0082	0.1319	600.0	600.0	600.0	0	0	1	-1.3808282862652839	2.635605857601787	130.2674726589629	-10.84781796752093	-130.1543760797224	-1.686461670570326				
	17	27	0.0013	0.0173	0.3216	600.0	600.0	600.0	0	0	1	-3.328884789710083	5.9415723355066365	-86.48988251901832	25.95332421594098	86.60510041713407	-58.03075366353638				
	19	20	0.0007	0.0138	0.0	900.0	900.0	2500.0	1.06	0.0	1	0.23491269600363754	4.950355349930313	45.63281860694396	57.3404130068616	-45.59469859405256	-56.58890418128806				
	19	33	0.0007	0.0142	0.0	900.0	900.0	2500.0	1.07	0.0	1	-6.090541362240652	0.17188733853924698	-580.1969185499698	-163.56927907631282	582.8252323833773	216.88650255400455				
	20	34	0.0009	0.018	0.0	900.0	900.0	2500.0	1.009	0.0	1	-5.9874089591171025	0.1546986046853223	-498.40530140594745	-25.811095818711948	500.7586015166698	72.87709803315818				
	21	22	0.0008	0.014	0.2565	900.0	900.0	900.0	0	0	1	-6.119189251997192	1.105808544602489	-403.8436005528771	-10.812765854162983	405.1015150664066	6.107327605760707				
	22	23	0.0006	0.0096	0.1846	600.0	600.0	600.0	0	0	1	-1.6615776058793874	2.2574537128154435	266.73584835715906	6.629841944333685	-266.3262523203922	-19.33477949857958				
	22	35	0.0	0.0143	0.0	900.0	900.0	2500.0	1.025	0.0	1	-6.543178020394001	0.13750987083139757	-671.8373634235656	-12.737169550094393	671.8373634235656	77.58478174213114				
	23	24	0.0022	0.035	0.361	600.0	600.0	600.0	0	0	1	-3.0710537819012127	9.33348248268111	173.41249392877282	-51.23446733374204	-172.7543327380946	23.80779029066999				
	23	36	0.0005	0.0272	0.0	900.0	900.0	2500.0	0	0	1	-10.284592422598276	0.18334649444186343	-105.08624160838065	2.8892468323216023	105.13935290817946	7.876733466710098e-6				
	25	26	0.0032	0.0323	0.531	600.0	600.0	600.0	0	0	1	-6.193673765364199	11.390400967200765	-20.832000390502163	-7.946527501192096	20.857401489425794	-47.68734839145745				
	25	37	0.0006	0.0232	0.0	900.0	900.0	2500.0	1.025	0.0	1	-8.703228908037206	0.23491269600363754	-78.62910949972486	1.4240397138230783	78.66594609708463	0.00030871742395357967				
	26	27	0.0014	0.0147	0.2396	600.0	600.0	600.0	0	0	1	-3.019487580339438	4.910248304271155	312.713192464597	-8.882463502754161	-311.4051004171341	-2.3692463364636245				
	26	28	0.0043	0.0474	0.7802	600.0	600.0	600.0	0	0	1	-6.640580845566241	7.041651302157817	-202.59427741542288	22.408607382403712	204.44467057519097	-82.39956906546493				
	26	29	0.0057	0.0625	1.029	600.0	600.0	600.0	0	0	1	-10.324699468257435	7.64898656499649	-242.1763165385999	20.561204511807894	245.67159230668926	-87.82322275321918				
	28	29	0.0014	0.0151	0.249	600.0	600.0	600.0	0	0	1	-3.7643327140095084	0.6990085100596044	-369.244670575191	60.319569065464925	371.20063073605616	-64.36530487336498				
	29	38	0.0008	0.0156	0.0	1200.0	1200.0	2500.0	1.025	0.0	1	-8.94960075994346	0.24064227395494575	-843.6722230427455	130.66852762658417	849.764837433483	-11.862547007201407				
];

%%-----  OPF Data  -----%%
%% cost data
%    1    startup    shutdown    n    x1    y1    ...    xn    yn
%    2    startup    shutdown    n    c(n-1)    ...    c0
mpc.gencost = [
	2	0.0	0.0	2	6.724778000000001	0.0
	2	0.0	0.0	2	14.707625	0.0
	2	0.0	0.0	2	24.804734	0.0
	2	0.0	0.0	2	34.844643	0.0
	2	0.0	0.0	2	24.652994	0.0
	2	0.0	0.0	2	32.306483	0.0
	2	0.0	0.0	2	18.157477	0.0
	2	0.0	0.0	2	31.550181	0.0
	2	0.0	0.0	2	22.503168000000002	0.0
	2	0.0	0.0	2	27.434444	0.0
];

%column_names% pf 
mpc.load_data = {
	0.9109409447828044
	0.999972224451119
	0.9384712337524265
	0.9411027374711993
	0.9472583511480259
	0.09713606937401495
	0.09647962731984765
	0.9021819276640016
	0.9952152687424872
	0.9824472495793353
	0.9887220293322959
	0.9220781944803617
	0.9462471496400356
	0.9581503989055199
	0.9785126902171042
	0.9926039598508953
	0.965748365291303
	0.9911436587590456
	0.9955285546245434
	0.894427190999916
	0.9753061185369931
};

