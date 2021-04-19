SRR6703982_frozen_501_1_val_1.Aligned.sortedByCoord.out <- c(0.0114192258415,0.150153356369,0.313055674186,0.476681096401,0.576726715068,0.688264355321,0.756537057887,0.789440734558,0.82726540358,0.826396707691,0.854680281089,0.857999767054,0.882551345265,0.908495748729,0.915697674419,0.888622510386,0.921099118686,0.970784641068,0.984756571029,0.989769771324,0.988736071747,0.995025624102,0.991740109485,0.995020771053,0.9908228831,0.998413052762,0.992647629771,0.990031836006,0.998383934464,0.994292813604,0.986479403657,0.984125674574,0.979816166479,0.989022401677,0.990361843382,1.0,0.989556237139,0.989172846217,0.996093295027,0.991599371045,0.990512287922,0.967523391699,0.973002484761,0.96525216446,0.962932406724,0.940719998447,0.918871568894,0.917638894281,0.892446713515,0.897867569981,0.883575338743,0.882973560585,0.882279574485,0.881580735334,0.871937725667,0.865653026362,0.859644950887,0.853238925341,0.854398804209,0.834578949412,0.831647707419,0.818190200722,0.825765811236,0.817423418876,0.813298326668,0.806581705944,0.799933027915,0.78489342703,0.790886943355,0.766374189541,0.75254299802,0.720575959933,0.722808362775,0.709976899484,0.710685444733,0.696456303141,0.67910179757,0.665309430446,0.653045773964,0.642432154366,0.636259075203,0.622583181271,0.599133245331,0.579298831386,0.53418488178,0.533524867026,0.507813409947,0.497626858718,0.460680591684,0.440365725822,0.423530496564,0.383512249097,0.338378887293,0.298938152735,0.260065224987,0.217644717941,0.169366579959,0.11681290523,0.055664479559,0.0)


pdf("/scratchLocal/ilb4001/RSeqc_2/SRR6703982*_1/rseqc_geneBody_coverage.out.geneBodyCoverage.curves.pdf")
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(1)
plot(x,SRR6703982_frozen_501_1_val_1.Aligned.sortedByCoord.out,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1])
dev.off()