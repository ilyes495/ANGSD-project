SRR6703988_frozen_510_1_val_1.Aligned.sortedByCoord.out <- c(0.0089526006481,0.137474597682,0.304937661339,0.457376533501,0.555815195954,0.665777897612,0.734787319562,0.769004928174,0.798014749572,0.800905745542,0.833136106493,0.84315722724,0.882757381027,0.901786026354,0.909660120733,0.89142038277,0.927405543322,0.970795448304,0.989080124028,0.998017745423,0.991127288705,1.0,0.989125061789,0.993419114526,0.988495933133,0.98858081557,0.986618533331,0.984241825072,0.998037717761,0.995341452089,0.986928104575,0.982698961938,0.979773014375,0.981230995072,0.979678145768,0.988840455968,0.973257039001,0.97273276512,0.975059542534,0.969706955866,0.973132211887,0.951906609346,0.955651422779,0.946144589743,0.945520454171,0.926736469989,0.907832651777,0.907582997548,0.884290257993,0.888274739486,0.87574709028,0.870983687593,0.86575592804,0.864113203214,0.860438292964,0.856983078436,0.847675968783,0.840465954653,0.841164986494,0.819350199973,0.815165995097,0.803851665443,0.812035331066,0.80608357425,0.797450531015,0.795383393999,0.784383628674,0.774037957429,0.778307044743,0.759388247278,0.743100805385,0.713571703191,0.713022463888,0.702012712393,0.700295091299,0.686668963486,0.665837814627,0.651597537411,0.63374226696,0.622707550043,0.619302266361,0.607428711235,0.588300204217,0.564757810433,0.524034212616,0.520509094904,0.497490975,0.486715898481,0.453586782307,0.43320501106,0.415304802848,0.379699116723,0.33527065015,0.295141229397,0.259910024616,0.218182816799,0.169275553359,0.116678400415,0.0578448848345,0.0)


pdf("/scratchLocal/ilb4001/RSeqc_2/SRR6703988*_1/rseqc_geneBody_coverage.out.geneBodyCoverage.curves.pdf")
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(1)
plot(x,SRR6703988_frozen_510_1_val_1.Aligned.sortedByCoord.out,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1])
dev.off()
