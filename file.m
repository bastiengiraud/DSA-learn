%% MATPOWER Case Format : Version 2
function mpc = pglib_opf_case39_epri
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100.0;
%% bus data
%    bus_i    type    Pd    Qd    Gs    Bs    area    Vm    Va    baseKV    zone    Vmax    Vmin
mpc.bus = [
	1	1	78.08	35.36	0.0	0.0	2	1.0058894633016264	-14.488993987175416	345.0	1	1.06	0.94				
	2	1	0.0	0.0	0.0	0.0	2	1.043281035341058	-8.59938685595842	345.0	1	1.06	0.94				
	3	1	257.6	1.9200000000000002	0.0	0.0	2	1.0376096766562346	-11.34694724555093	345.0	1	1.0554	0.94				
	4	1	400.0	147.20000000000002	0.0	0.0	1	1.027309508860541	-11.460797670206427	345.0	1	1.0479	0.94				
	5	1	0.0	0.0	0.0	0.0	1	1.0313703817057018	-10.09837888538134	345.0	1	1.0567	0.9439				
	6	1	0.0	0.0	0.0	0.0	1	1.034841810333174	-9.37107081284452	345.0	1	1.06	0.9469				
	7	1	187.04000000000002	67.2	0.0	0.0	1	1.0226325712362054	-11.47321801416334	345.0	1	1.0506	0.94				
	8	1	417.6	141.28	0.0	0.0	1	1.019818384295193	-12.073938871350201	345.0	1	1.0489	0.94				
	9	1	5.2	-53.279999999999994	0.0	0.0	1	1.0155097071729087	-15.135173263959318	345.0	1	1.06	0.9419				
	10	1	0.0	0.0	0.0	0.0	1	1.0488006657007456	-7.292746536797321	345.0	1	1.06	0.9552				
	11	1	0.0	0.0	0.0	0.0	1	1.0431373374755202	-8.006611167105458	345.0	1	1.0593	0.951				
	12	1	6.824	70.4	0.0	0.0	1	1.0347717433689512	-8.12328818466921	345.0	1	1.0497	0.94				
	13	1	0.0	0.0	0.0	0.0	1	1.0439438140388158	-8.141587554641134	345.0	1	1.0573	0.9498				
	14	1	0.0	0.0	0.0	0.0	1	1.0358739574397509	-10.17042262385697	345.0	1	1.0531	0.9423				
	15	1	256.00000000000006	122.40000000000002	0.0	0.0	3	1.029292409864174	-12.407166877997666	345.0	1	1.0466	0.94				
	16	1	263.2	25.839999999999996	0.0	0.0	3	1.037389197840436	-12.127799957154824	345.0	1	1.0546	0.9449				
	17	1	0.0	0.0	0.0	0.0	2	1.0403946974646303	-11.647480248637226	345.0	1	1.0554	0.9406				
	18	1	126.40000000000002	24.0	0.0	0.0	2	1.0385317334153081	-11.872070220791777	345.0	1	1.0546	0.94				
	19	1	0.0	0.0	0.0	0.0	3	1.0563390320691342	-10.569178979109868	345.0	1	1.06	0.9893				
	20	1	544.0	82.4	0.0	0.0	3	0.9903985188132509	-14.431855230472916	345.0	1	1.0067	0.94				
	21	1	219.20000000000002	92.0	0.0	0.0	3	1.025271459936635	-11.396353113265647	345.0	1	1.0519	0.94				
	22	1	0.0	0.0	0.0	0.0	3	1.0242769435451415	-8.983100276574842	345.0	1	1.06	0.9447				
	23	1	198.00000000000003	67.68	0.0	0.0	3	1.0301653812315785	-8.144454179470895	345.0	1	1.0572	0.94				
	24	1	246.88000000000002	-73.76	0.0	0.0	3	1.039939470979061	-12.234629172229276	345.0	1	1.0587	0.9459				
	25	1	179.20000000000002	37.760000000000005	0.0	0.0	2	1.0599999697748832	-6.786705484681473	345.0	1	1.06	0.94				
	26	1	111.19999999999999	13.600000000000001	0.0	0.0	2	1.0584122304182162	-7.349208170726318	345.0	1	1.06	0.94				
	27	1	224.8	60.400000000000006	0.0	0.0	2	1.0451368042294853	-10.225325004323533	345.0	1	1.054	0.94				
	28	1	164.8	22.080000000000002	0.0	0.0	3	1.0596768364227631	-2.1735915359535274	345.0	1	1.06	0.9414				
	29	1	226.8	21.519999999999996	0.0	0.0	3	1.0599998556294137	0.7667054939369136	345.0	1	1.06	0.9457				
	30	2	0.0	0.0	0.0	0.0	2	1.0410631455780794	-5.920261071939391	345.0	1	1.06	0.94				
	31	3	7.359999999999999	3.6799999999999997	0.0	0.0	1	1.013858514904257	1.3357164091245737e-22	345.0	1	1.06	0.94				
	32	2	0.0	0.0	0.0	0.0	1	1.0286803680949648	0.830241933523871	345.0	1	1.048	0.94				
	33	2	0.0	0.0	0.0	0.0	3	1.011901814602239	-5.361570847713273	345.0	1	1.0266	0.94				
	34	2	0.0	0.0	0.0	0.0	3	0.9966722968720466	-13.8195681962058	345.0	1	1.0282	0.94				
	35	2	0.0	0.0	0.0	0.0	3	0.9862013688287783	-7.73013719970755	345.0	1	1.06	0.94				
	36	2	0.0	0.0	0.0	0.0	3	1.0593826896458347	0.1286977130218837	345.0	1	1.06	0.94				
	37	2	0.0	0.0	0.0	0.0	2	1.0295944616779054	0.27225529833532663	345.0	1	1.06	0.94				
	38	2	0.0	0.0	0.0	0.0	3	1.0375631613015914	7.977273217711297	345.0	1	1.06	0.94				
	39	2	883.1999999999999	200.0	0.0	0.0	1	0.9819134982454834	-17.171539551451694	345.0	1	1.06	0.94				
];

%% generator data
%    bus    Pg    Qg    Qmax    Qmin    Vg    mBase    status    Pmax    Pmin    Pc1    Pc2    Qc1min    Qc1max    Qc2min    Qc2max    ramp_agc    ramp_10    ramp_30    ramp_q    apf
mpc.gen = [
	30	278.932834509921	140.00005993442485	400.0	140.0	1.0	100.0	1	1040.0	0.0															
	31	646.0000057752618	245.4796883368223	300.0	-100.0	1.0	100.0	1	646.0	0.0															
	32	713.954241565489	299.99923780289487	300.0	150.0	1.0	100.0	1	725.0	0.0															
	33	648.5516988021433	172.93438132497795	250.0	0.0	1.0	100.0	1	652.0	0.0															
	34	63.49544011179571	80.8574573648014	167.0	0.0	1.0	100.0	1	508.0	0.0															
	35	152.95952288591928	-88.64968580077489	300.0	-100.0	1.0	100.0	1	687.0	0.0															
	36	580.0	144.88796361563834	240.0	0.0	1.0	100.0	1	580.0	0.0															
	37	564.0	0.00011145780979029671	250.0	0.0	1.0	100.0	1	564.0	0.0															
	38	856.9253350636552	32.762801242430385	300.0	-150.0	1.0	100.0	1	865.0	0.0															
	39	554.8117318180975	-99.88558377636643	300.0	-100.0	1.0	100.0	1	1100.0	0.0															
];

%% branch data
%    f_bus    t_bus    r    x    b    rateA    rateB    rateC    ratio    angle    status    angmin    angmax
mpc.branch = [
	1	2	0.0035	0.0411	0.6987	600.0	600.0	600.0	0	0	1	-16.157409822689214	4.365938398896874	-266.7153704136079	-90.6696114566523	269.2819656929669	47.436576761736525				
	1	39	0.001	0.025	0.75	1000.0	1000.0	1000.0	0	0	1	-4.073729923380153	8.250592249883853	188.6353704136079	55.309611456652306	-188.19774679020964	-118.46781074973758				
	2	3	0.0013	0.0151	0.2572	500.0	500.0	500.0	0	0	1	-2.538203032429547	4.938896194027696	345.17471666658565	3.7112125865933834	-343.7479260692656	-14.98128042796059				
	2	25	0.007	0.0086	0.146	500.0	500.0	500.0	0	0	1	-3.409098881028398	3.4320171928336314	-340.81128734295487	73.07351765335238	348.70348798186296	-79.52524279859688				
	2	30	0.0	0.0181	0.0	900.0	900.0	2500.0	1.025	0.0	1	-10.679933301238545	0.1604281826366305	-273.6453950165977	-124.22130700168228	273.6453950165977	140.00005993442485				
	3	4	0.0013	0.0213	0.2214	500.0	500.0	500.0	0	0	1	-4.285724307578558	6.348372370049521	12.958879986899738	37.47693449716758	-12.927391354516446	-60.56223062188655				
	3	18	0.0011	0.0133	0.2138	500.0	500.0	500.0	0	0	1	-3.1512678732195276	3.179915762976069	73.18904608236578	-24.415654069206997	-73.13261525499087	2.059058467438927				
	4	5	0.0008	0.0128	0.1342	600.0	600.0	600.0	0	0	1	-3.311696055856158	2.8705185536054243	-197.93025000642737	-24.962751756037694	198.2296431413432	15.533950086871142				
	4	14	0.0008	0.0129	0.1382	500.0	500.0	500.0	0	0	1	-3.40336930307709	0.26929016371148695	-189.14235863905623	-61.67501762207579	189.4359616194309	51.70212387362951				
	5	6	0.0002	0.0026	0.0434	1200.0	1200.0	1200.0	0	0	1	-1.1172677005051053	0.2578310078088704	-528.2863071933251	-96.06854308467872	528.8275701040333	98.47283023652376				
	5	8	0.0008	0.0112	0.1476	900.0	900.0	900.0	0	0	1	0.6130648407899808	3.7929806037660496	330.05666405198195	80.53459299780758	-329.1786228018188	-83.76772297591381				
	6	7	0.0006	0.0092	0.113	900.0	900.0	900.0	0	0	1	0.6474423084978302	3.477853816444097	429.57242145878706	111.0079779267332	-428.4617540670114	-105.93695804282423				
	6	11	0.0007	0.0082	0.1389	480.0	480.0	480.0	0	0	1	-2.5553917662834715	0.8479775367936184	-319.75998578755855	-81.09779975045087	320.4637936558477	74.34792148675984				
	6	31	0.0	0.025	0.0	1800.0	1800.0	1800.0	1.07	0.0	1	-11.092462913732737	0.29220847551671986	-638.6400057752618	-128.3830084128061	638.6400057752618	241.7996883368223				
	7	8	0.0004	0.0046	0.078	900.0	900.0	900.0	0	0	1	-0.24064227395494575	1.2776958831417358	241.4217540670114	38.73695804282424	-241.19180979553434	-44.22724588013355				
	8	9	0.0023	0.0363	0.3804	900.0	900.0	900.0	0	0	1	-6.600473799907084	11.676879864766178	152.7704325973531	-13.285031143952633	-152.2533673369483	-17.950256513431373				
	9	39	0.001	0.025	1.2	900.0	900.0	900.0	0	0	1	-4.772738433439757	7.849521793292279	147.0533673369483	71.23025651343137	-146.6718741980129	-181.41777302662885				
	10	11	0.0004	0.0043	0.0729	600.0	600.0	600.0	0	0	1	-0.4526366581533504	1.3349916626548182	327.1937002620278	105.66132501377596	-326.76066298309036	-108.98186561857337				
	10	13	0.0004	0.0043	0.0729	600.0	600.0	600.0	0	0	1	-0.6302535746439055	1.5126085791453734	385.16216893358484	81.41788708768183	-384.5961683780726	-83.31520766784006				
	10	32	0.0	0.02	0.0	900.0	900.0	2500.0	1.07	0.0	1	-9.963736057325017	0.13178029288008933	-712.3558691956126	-187.07921210145778	712.3558691956126	299.99923780289487				
	12	11	0.0016	0.0435	0.0	500.0	500.0	500.0	1.006	0.0	1	-0.7677634454753032	0.888084582452776	-6.278648645807633	-34.13856935530012	6.296869327242635	34.63394413181352				
	12	13	0.0016	0.0435	0.0	500.0	500.0	500.0	1.006	0.0	1	-1.008405719430249	0.6531718864491385	-0.545351354192367	-36.261430644699885	0.5652404440928857	36.80216527636952				
	13	14	0.0009	0.0101	0.1723	600.0	600.0	600.0	0	0	1	-1.5126085791453734	3.552338329811104	384.0309279339797	46.513042391470535	-382.7871954396775	-51.18858322400723				
	14	15	0.0018	0.0217	0.366	600.0	600.0	600.0	0	0	1	-4.927437038125079	6.508800552686152	193.35123382024653	-0.5135406496222827	-192.7179766695656	-30.876534877172627				
	15	16	0.0009	0.0094	0.171	600.0	600.0	600.0	0	0	1	-3.684118622691193	1.4495832216809827	-63.28202333043447	-91.52346512282739	63.373813231250594	74.22261537345648				
	16	17	0.0007	0.0089	0.1342	600.0	600.0	600.0	0	0	1	-2.7960340402384176	3.50077212824933	-103.7406520838989	-33.66792445482478	103.81520388465933	20.131608297847407				
	16	19	0.0016	0.0195	0.304	600.0	600.0	2500.0	0	0	1	-6.90414143132642	6.778090716397639	-159.8780741881447	-101.97258060472014	160.36707769459267	74.61347766361375				
	16	21	0.0008	0.0135	0.2548	600.0	600.0	600.0	0	0	1	-4.004974987964455	2.9450030669724314	-94.68758013602144	85.65973765623512	94.82763286496636	-110.3988876835206				
	16	24	0.0003	0.0059	0.068	600.0	600.0	600.0	0	0	1	-0.7219268218648373	1.5412564689019146	31.73249317681442	-50.081847970146654	-31.72367854718751	42.919190819241926				
	17	18	0.0007	0.0082	0.1319	600.0	600.0	600.0	0	0	1	-1.3808282862652839	2.635605857601787	53.28812982232189	12.050479777211923	-53.26738474500916	-26.059058467438923				
	17	27	0.0013	0.0173	0.3216	600.0	600.0	600.0	0	0	1	-3.328884789710083	5.9415723355066365	-157.10333370698123	-32.182088075059326	157.4023832369078	1.1920555820312482				
	19	20	0.0007	0.0138	0.0	900.0	900.0	2500.0	1.06	0.0	1	0.23491269600363754	4.950355349930313	483.62967072268043	36.10972728283628	-481.9718279758266	-3.4265417020043953				
	19	33	0.0007	0.0142	0.0	900.0	900.0	2500.0	1.07	0.0	1	-6.090541362240652	0.17188733853924698	-643.996748417273	-110.72320494645002	647.063496548609	172.93438132497795				
	20	34	0.0009	0.018	0.0	900.0	900.0	2500.0	1.009	0.0	1	-5.9874089591171025	0.1546986046853223	-62.028172024173465	-78.9734582979956	62.122371977513815	80.8574573648014				
	21	22	0.0008	0.014	0.2565	900.0	900.0	900.0	0	0	1	-6.119189251997192	1.105808544602489	-314.0276328649664	18.398887683520606	314.7858631104118	-32.066524250574076				
	22	23	0.0006	0.0096	0.1846	600.0	600.0	600.0	0	0	1	-1.6615776058793874	2.2574537128154435	-164.08903191835213	-61.07760733713643	164.25812163018782	44.30419868617452				
	22	35	0.0	0.0143	0.0	900.0	900.0	2500.0	1.025	0.0	1	-6.543178020394001	0.13750987083139757	-150.6968311920597	93.1441315877105	150.6968311920597	-88.64968580077489				
	23	24	0.0022	0.035	0.361	600.0	600.0	600.0	0	0	1	-3.0710537819012127	9.33348248268111	216.14962062252957	-53.71432481738654	-215.15632145281253	30.840809180758082				
	23	36	0.0005	0.0272	0.0	900.0	900.0	2500.0	0	0	1	-10.284592422598276	0.18334649444186343	-578.4077422527174	-58.26987386878798	579.9999865495345	144.88796361563834				
	25	26	0.0032	0.0323	0.531	600.0	600.0	600.0	0	0	1	-6.193673765364199	11.390400967200765	34.29608334768088	-27.851405345792653	-34.262473024639746	-31.383196488911096				
	25	37	0.0006	0.0232	0.0	900.0	900.0	2500.0	1.025	0.0	1	-8.703228908037206	0.23491269600363754	-562.1995713295439	69.61664814438952	564.0000047675318	0.00011145780979029671				
	26	27	0.0014	0.0147	0.2396	600.0	600.0	600.0	0	0	1	-3.019487580339438	4.910248304271155	384.1048096142547	55.06121482339743	-382.20238323690785	-61.59205558203126				
	26	28	0.0043	0.0474	0.7802	600.0	600.0	600.0	0	0	1	-6.640580845566241	7.041651302157817	-211.09497529984856	-17.72699305660388	212.8313367903609	-50.63798800122922				
	26	29	0.0057	0.0625	1.029	600.0	600.0	600.0	0	0	1	-10.324699468257435	7.64898656499649	-249.94736128976643	-19.55102527788245	253.19995782502662	-60.229908608031735				
	28	29	0.0014	0.0151	0.249	600.0	600.0	600.0	0	0	1	-3.7643327140095084	0.6990085100596044	-377.6313367903609	28.55798800122922	379.4318373170375	-37.107411659453675				
	29	38	0.0008	0.0156	0.0	1200.0	1200.0	2500.0	1.025	0.0	1	-8.94960075994346	0.24064227395494575	-859.4317951420642	75.81732026748541	865.0000065015469	32.762801242430385				
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

