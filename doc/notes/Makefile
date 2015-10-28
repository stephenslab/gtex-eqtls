SHELL=/bin/bash
TITLE="GTEx multi-tissue eQTL analysis"
wiki:
	tigernotes dokuwiki {Information,SummaryStats,QuasiPoissonModel,MultiTissueAnalysis}.notes -o today -t $(TITLE) --inline_toc --permission "\@user" -a "Gao Wang"
	scp today.txt gw@tigerwang.org:./tigerwang_data/pages/research
	# google-chrome http://tigerwang.org/research/today

pdf:
	tigernotes doc {Information,SummaryStats,QuasiPoissonModel,MultiTissueAnalysis}.notes -o Notes --footnote --font_size 12 -a "Gao Wang" -t $(TITLE) 

sumstats:
	tigernotes doc SumstatsDB.notes -o SumstatsDB.pdf --footnote --font_size 12 -a "Gao Wang" -t "GTEx summary statistics database"
