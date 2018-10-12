The present data and the associated source code are freely available under the GNU GPL v3 licence. They correspond to the paper entitled *Convergence rate of individual and global quantities in Direct Numerical Simulations* we have submitted to the journal [Physics of Fluids](https://aip.scitation.org/journal/phf).

The python scripts used to generate the Figures can be used to estimate the convergence rate with similar to

    polyfit(log(abs(array(t[10000:]))),log(abs(array(err_EX0[10000:]))),1)

# Acknowledgements

The author and coworker thank the Slovenian Research Agency for funding the study (research projects P2-0026). We also thank Framasoft and the Institut Jozef Stefan for providing the gitlab services that host the present project at [https://repo.ijs.si/CFLAG/convergence_rate](https://repo.ijs.si/CFLAG/convergence_rate) and [https://framagit.org/CFLAG/convergence_rate](https://framagit.org/CFLAG/convergence_rate).
