SRR6703989_frozen_510f_1_val_1.Aligned.sortedByCoord.out <- c(0.0,0.103314753125,0.252857943864,0.388251591411,0.472924568752,0.574952999769,0.646518684653,0.685035786141,0.709185659158,0.715980078499,0.752340116758,0.764952669943,0.805640027705,0.82546258122,0.840007915828,0.823575975461,0.85490286619,0.905452026782,0.925030508922,0.933124443418,0.935532174544,0.949061644513,0.950084105676,0.956871928494,0.955216201062,0.961753355981,0.959998680695,0.963976384445,0.980804116231,0.977083676902,0.973759028992,0.972960849632,0.970480556747,0.979399056697,0.980230218675,0.989861143178,0.985995580329,0.993640951219,1.0,0.996761106897,0.997084336555,0.980441307431,0.980454500478,0.977123256044,0.97989379597,0.965625515353,0.948929714041,0.949213364557,0.931897490023,0.935723473729,0.927016062535,0.928876282199,0.927741680135,0.940215706323,0.934232659389,0.93559813978,0.929047791814,0.922504040371,0.927899996702,0.912048550414,0.913018239388,0.899521752037,0.915406180943,0.908420462416,0.904983673604,0.90348626274,0.904099739437,0.890491111184,0.903934826347,0.889910617105,0.872179161582,0.847488373627,0.852303835878,0.839664896599,0.844262673571,0.826023285728,0.808516112009,0.794280814011,0.782684125466,0.773890959464,0.771720703189,0.759048781292,0.737082357598,0.718631880999,0.673795309872,0.673788713348,0.646274613279,0.639664896599,0.597895708961,0.57866024605,0.556838945876,0.513381048188,0.465457304001,0.415178600877,0.369200831162,0.316705696098,0.255120551469,0.184102378047,0.100900425476,0.0229493057159)


pdf("/scratchLocal/ilb4001/RSeqc_2/SRR6703989*_1/rseqc_geneBody_coverage.out.geneBodyCoverage.curves.pdf")
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(1)
plot(x,SRR6703989_frozen_510f_1_val_1.Aligned.sortedByCoord.out,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1])
dev.off()
