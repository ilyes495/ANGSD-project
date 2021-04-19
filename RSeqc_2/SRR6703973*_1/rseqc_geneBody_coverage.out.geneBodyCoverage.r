SRR6703973_time_104SR_1_val_1.Aligned.sortedByCoord.out <- c(0.0,0.0740008307598,0.166716331654,0.308902674685,0.399207166206,0.489854797638,0.567088367557,0.619453323942,0.6738455148,0.686794531433,0.728151920681,0.743733181629,0.779482942334,0.810342959311,0.832073648661,0.834936157917,0.859565476513,0.895477777176,0.914969027108,0.918418486211,0.925014899496,0.934415127052,0.941368225244,0.951174802695,0.952746022286,0.96649870871,0.981818099727,0.992089722057,0.995773961099,0.990703617417,0.994455581442,0.988328728035,0.97569124632,0.987105163353,0.990125697567,0.993399071716,0.985258528833,0.987443788265,1.0,0.994405916454,0.993096566795,0.970887287569,0.978368640624,0.97695996099,0.979637355294,0.964187029311,0.95057882285,0.942627909917,0.915023207094,0.922644525112,0.922053060266,0.929633743295,0.925416734392,0.932293077604,0.925899839266,0.923755214824,0.923041845009,0.917185891532,0.918219826263,0.8849984649,0.884754654964,0.875399577396,0.880826605985,0.867263549511,0.873548427877,0.878736161529,0.870392443698,0.853610193061,0.87323237796,0.862838850662,0.846756424843,0.817061277564,0.819124632028,0.816555597696,0.805399035596,0.784869335934,0.764664716187,0.752867024255,0.743069476802,0.729063950443,0.723623376858,0.707531921042,0.684496397031,0.661758862943,0.61834714923,0.61619349479,0.583622293258,0.579703274277,0.544332773473,0.525410413393,0.507738707988,0.465816943888,0.424762059562,0.381792815734,0.340403821495,0.296509002908,0.239782557657,0.175163894457,0.107994256921,0.0411813042929)


pdf("/scratchLocal/ilb4001/RSeqc_2/SRR6703973*_1/rseqc_geneBody_coverage.out.geneBodyCoverage.curves.pdf")
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(1)
plot(x,SRR6703973_time_104SR_1_val_1.Aligned.sortedByCoord.out,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1])
dev.off()