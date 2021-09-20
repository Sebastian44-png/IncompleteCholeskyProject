# Gnuplot script file for plotting data in file "force.dat"
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set terminal png size 1280, 480font "Helvetica,30"
set key font ",15s"
set title "Fehler nach der k-ten Iteration für das vorkonditionierte konjugierte Gradientenverfahren"
set term svg
set output "plot.svg"
set xlabel "k-th Iteration"
set ylabel "Error"
set logscale y
set format y "10^{%T}"
#set format y "%.2t*10^{%T}"
plot    "< head -10 Analyse_Problem_1.dat" using 1:2 title 'Error Jacobi' with linespoints , \
        "< head -10 Analyse_Problem_1.dat" using 1:3 title 'Error Gauss-Seidel' with linespoints , \
        "< head -10 Analyse_Problem_1.dat" using 1:4 title 'Error ICF' with linespoints , \
        "< head -10 Analyse_Problem_1.dat" using 1:5 title 'Error ICNE' with linespoints , \
        "< head -10 Analyse_Problem_1.dat" using 1:6 title 'Error Multigrid' with linespoints
            
