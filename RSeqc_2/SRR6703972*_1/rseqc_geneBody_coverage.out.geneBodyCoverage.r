SRR6703972_time_103FR_1_val_1.Aligned.sortedByCoord.out <- c(0.0,0.083929748029,0.192783950197,0.326789812101,0.407879490151,0.51095130984,0.579510599141,0.62227069275,0.64835389544,0.664428691523,0.726080930634,0.743712085066,0.785404312362,0.808215754436,0.831174880149,0.838763547134,0.864789948424,0.908265739668,0.922522890964,0.919716901824,0.927294208529,0.93960875196,0.944686797083,0.956989980233,0.946845250267,0.956740054075,0.974439370186,0.981823552132,0.990661849908,0.990173357872,0.991661554541,0.981846272692,0.976359257492,0.986628950537,0.981948515211,0.983084543203,0.97800649808,0.977143116806,0.996126144548,1.0,0.994899234317,0.980187671824,0.981891713812,0.976859109808,0.975654920137,0.952479949106,0.938813532365,0.927112444051,0.897121305069,0.902994569786,0.895928475678,0.890804989435,0.893838184173,0.897973326063,0.886260877468,0.879410628678,0.891770613228,0.886703928385,0.874491627474,0.85573580533,0.856076613728,0.851203053643,0.873866812078,0.870333765024,0.87774066753,0.868243473519,0.862199804603,0.846761184196,0.861063776611,0.850703201327,0.836423329471,0.797605252993,0.794503896576,0.782041669507,0.778713107491,0.76928407516,0.754754277145,0.738520437144,0.7224797219,0.71289164565,0.696510122009,0.687274214437,0.663133619612,0.640208574739,0.603162701929,0.599357008157,0.573750937223,0.572467225592,0.543271306205,0.531354372572,0.513336968623,0.478381387317,0.434826074114,0.398677663418,0.354497534819,0.306932042806,0.250164724059,0.189000976984,0.117442573785,0.0430213800468)


pdf("/scratchLocal/ilb4001/RSeqc_2/SRR6703972*_1/rseqc_geneBody_coverage.out.geneBodyCoverage.curves.pdf")
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(1)
plot(x,SRR6703972_time_103FR_1_val_1.Aligned.sortedByCoord.out,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1])
dev.off()
