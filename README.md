# systematics of hal qcd pot

This is the repository of the published paper
"Systematics of the HAL QCD Potential at Low Energies in Lattice QCD",
Phys. Rev. D99, 014514(2019), DOI:10.1103/PhysRevD.99.014514
[arXiv:1805.02365 [hep-lat]](https://arxiv.org/abs/1805.02365).

In this paper, we discuss the convergence of the derivative expansion
of the non-local HAL QCD potential at the next-to-next-to-leading
order in $\Xi\Xi$ 1S0 channel at the pion mass 510 MeV.


* notes

  + Effective mass wall vs. smeared.ipynb

  plot the effective mass of Xi at L = 64 for wall and smeared source

  + Finite volume method.ipynb

  analyze the finite volume spectra from the fitted parameter (Run "N2LO analysis.ipynb" and "Fitting and Scattering phase shift.ipynb" before using this notebook.), and fit the effective range expansion with the finite volume condition correctly

  + Fitting and Scattering phase shift.ipynb

  fit the potentials and evaluate the scattering phase shifts (Run "N2LO  analysis.ipynb" before using this notebook)

  + LO potential (L = 40, 48 and 64).ipynb

   plot volume dependence of the leading order potential for the wall source

  + LO potential.ipynb

  plot the leading order potentials

  + N2LO analysis.ipynb

   analyze the next-to-next-to-leading order potential by using both wall and smeared sources

  + Non-locality vs Energy dependence.ipynb

  plot energy-dependent local potential in Appendix
  (Run "N2LO analysis.ipynb" before using this notebook)

  + R-correlator.ipynb : plot R-correlator

  plot R-correlator for the wall and the smeared sources

  + luscher_lib.py

  python script for the finite volume formula (only S-wave
    but boosted frame is supported)

  + solve_eigenvalues.py

  sample python script for the finite volume spectra


* data
  + potential

  pickle files of the HAL QCD potential and H0-term 

  + Rcorr

  pickle files of R-correlator (jackknife samples)

  + meff_xi_L64.pkl

  pickle format file of the effective masses of Xi (See "Effective mass wall vs. smeared.ipynb")
