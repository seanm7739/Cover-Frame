
#summarize the results $1=current directory
#by seanm 2014-1-16
pysh=/home/hadoop/mgltools_i86Linux2_1.5.6/bin/pythonsh
sum=/home/hadoop/mgltools_i86Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/summarize_results4.py
$pysh $sum -d $1 -o summary.txt -t 1.0