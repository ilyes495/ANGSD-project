SRR6703969_time_101SR_1_val_1.Aligned.sortedByCoord.out <- c(0.0,0.0720645059697,0.166251458612,0.30942354048,0.398292869054,0.495378875668,0.571077293301,0.626672104075,0.680663956038,0.690863134774,0.733820008318,0.748129544913,0.786257353624,0.818377095051,0.843648931832,0.845017705225,0.86873503535,0.905869574873,0.928852047305,0.932942216767,0.937193893413,0.94911312367,0.949428062681,0.957196558282,0.954713385312,0.960289420876,0.974481864762,0.986259776232,0.991601626377,0.99244146374,0.998909826501,0.990378209451,0.976823718945,0.98303367021,0.991706606048,0.997125172106,0.98531092152,0.982250360363,1.0,0.990265154421,0.989913876294,0.965651459419,0.972220764171,0.974542429956,0.971974465714,0.95850880416,0.942366161015,0.938566704486,0.906907258537,0.91929889731,0.917643448663,0.925391755866,0.918693245366,0.924826480718,0.916193921677,0.916468483892,0.91778072977,0.912568085373,0.914158931145,0.885507556517,0.884235687435,0.871093040252,0.874880383741,0.858802343469,0.867447015549,0.868823864302,0.860227644377,0.84626131055,0.862084977005,0.85091271748,0.832573576617,0.797768778239,0.801035261056,0.797704175365,0.790315221648,0.771588463542,0.752042056471,0.739440458357,0.727295118042,0.709388008899,0.703206321391,0.688864483359,0.668486314285,0.643775714972,0.598048993205,0.595408350729,0.565275147678,0.55974352659,0.525071971639,0.506724755418,0.488385614555,0.449050539636,0.404103090036,0.361327912075,0.318912087601,0.275442428745,0.22040078008,0.161676767595,0.0956324419483,0.0333431583538)


pdf("/scratchLocal/ilb4001/RSeqc_2/SRR6703969*_1/rseqc_geneBody_coverage.out.geneBodyCoverage.curves.pdf")
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(1)
plot(x,SRR6703969_time_101SR_1_val_1.Aligned.sortedByCoord.out,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1])
dev.off()
